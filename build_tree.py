# tree_main.py - a phylogenetic tree pipeline
# Written by Allison Lopatkin, 2020, Barnard College
#
# Basic Pipeline Steps:
#   * Use BLAST to create a multi-fasta file based on a target accession id
#   * Refine the multi-fasta based on name and quality
#   - Create a phylogenetic tree using ete3
#   - Use Entrez to fetch the Taxonomy ID based on the NCBI accession number
#   - Use ete3's NCBI translator to get species name from Taxonomy ID
#   - Annotate the phylogenetic tree using the species name
#
#   * steps will be skipped if a fasta input is provided
#
# Usage:
# python tree_main.py -o output_dir -a accession -f fasta_file -k api_key -t input_tree (.nw file)
#   - output_dir is the FULL PATH to for the output. A folder will be created here by ete3.
#   - accession is the accession ID to be used for the analysis. This will be the query for BLAST, and
#     will also determine the rooting of the generated phyltogenetic trees.
#   - fasta_file should be the path to a .fasta or .fa file. If this is provided, the BLAST step will be skipped. If it
#     is not provided, the accession ID will be used to generate a .fasta file.
#   - input_tree is a pre-created .nw file that can be input. If it is provided, the ETE pipeline will be skipped. Note
#     that if a tree is provided, the fasta file must also be included.
#   - the -c flag will collapse the tree to the genus level instead of keeping it at the species level.

from ete3 import PhyloTree, TreeStyle, faces, AttrFace, CircleFace, NodeStyle, Tree
import os
import sys
import getopt
from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIWWW
import copy
import csv
import re
global node2labels


# --- SUPPORT FUNCTIONS --- #
def run_blast(target_acc, output_dir):
    """
    Perform a BLAST query using the target accession ID and write results to output_dir.
    :param target_acc: the target accession ID.
    :param output_dir: the output directory where BLAST outputs will be written.
    :return: no output; folder is created by BLAST at output_dir.
    """

    # this is somewhat slow- it's faster via the website
#    result_handle = NCBIWWW.qblast("blastp", "nr", target_acc, hitlist_size=100)

    # get full path
    write_path = os.path.join(output_dir, "all_genuses_fused.xml")

    # write results to XML
#    with open(write_path, 'w') as save_file:
#        blast_results = result_handle.read()
#        save_file.write(blast_results)

    # read the XML file of hits
    blast_qresult = SearchIO.read(write_path, "blast-xml")

    # filter the hits for name (not FLAG-tagged) and quality (%identity > 0.65)
#    filter_func_name = lambda h: "FLAG-tagged" not in h.description
#    filter_func_pcti = lambda h: h.hsps[0].ident_num/h.hsps[0].aln_span > 0.65
    # hit filter iterates over ever hit in blastxml file
#    blast_qresult_filtered = blast_qresult.hit_filter(filter_func_pcti).hit_filter(filter_func_name)

    # write output to TSV; for headers see http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    SearchIO.write(blast_qresult, output_dir+"blast_tab.tsv", "blast-tab")

    # collect the resulting records
    records = []
    for hit in blast_qresult:
        curr_hit = hit[0].hit

        # try to get rid of the extra pipes if they exist (this happens in BlastXML)
        try:
            curr_hit.id = curr_hit.id.split("|")[1]
        except IndexError:
            pass

        # add the current record
        records.append(curr_hit)

    # write the records to a fasta file
    write_path = os.path.join(output_dir, "blast_results.fasta")
    SeqIO.write(records, write_path, "fasta")  # saving FASTA from all records


def get_species_name(acc_id, fasta_file):
    """
    Get a species name given an by searching the fasta input.

    :param acc_id: the accession ID to convert.
    :param fasta_file: the fasta to search.
    :return: the species name- this is a string.
    """
    species_name = ''
    with open(fasta_file) as file:
        for line in file:
            print(line)
            if acc_id in line:
                find = line
                try:
                    species_name = find.split('[')[1].split(']')[0]
                except:
                    pass

    return species_name


# --- LAYOUT & DRAWING FUNCTIONS --- #
def layout(node):
    """Set ETE3 layout properties on a node."""

    if node.is_leaf():
        # Add node name to leaf nodes
        n = AttrFace("name", fsize=20, fgcolor="black")
        faces.add_face_to_node(n, node, 0)
    if "weight" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        c = CircleFace(radius=node.weight, color="RoyalBlue", style="sphere")
        # Let's make the sphere transparent
        c.opacity = 0.3
        # And place as a float face over the tree
        faces.add_face_to_node(c, node, 0, position="float")


def draw_circular_tree(t, output_name="output_annotated_circular.png"):
    """Create a PNG for a circular tree."""

    ts = TreeStyle()

    # Set our custom layout function
    ts.layout_fn = layout

    # Draw a tree
    ts.mode = "c"

    # We will add node names manually
    ts.show_leaf_name = False
    # Show branch data (true or false to print distances on image)
    ts.show_branch_length = False
    ts.show_branch_support = False

    style = NodeStyle()
    style["fgcolor"] = "#000000"  # change hex color code to change line color
    style["size"] = 4
    style["vt_line_color"] = "#000000"
    style["hz_line_color"] = "#000000"
    style["vt_line_width"] = 8
    style["hz_line_width"] = 8
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    t.children[0].img_style = style
    t.children[1].img_style = style
    for n in t.traverse():
        n.img_style = style

    t.render(output_name, tree_style=ts, dpi=300)


# collapse the tree nodes
def collapsed_leaf(node):
    global node2labels

    if len(node2labels[node]) == 1:
        return True
    else:
        return False


# --- MAIN PIPELINE CODE --- #
def main(argv):
    fasta = ''
    output_dir = ''
    accession = ''
    input_tree = ''
    do_collapse = False

    try:
        opts, args = getopt.getopt(argv, "hcf:o:a:t:", ["fasta=", "outputdir=", "accession=", "tree="])
    except getopt.GetoptError:
        print("Invalid arguments. Valid arguments are -h, -f, --fasta, -o, --outputdir, -a, --accession, \
         -t, --tree, -c")
        sys.exit(2)

    # parse the input arguments
    if len(opts) == 0:
        print('Usage: tree_main.py -f <fasta file> -o <output dir> -a <accession id> -t <input_tree> -c')
        sys.exit()

    for opt, arg in opts:
        if opt == '-h':
            print('Usage: tree_main.py -f <fasta file> -o <output dir> -a <accession id> -t <input_tree>')
            sys.exit()
        elif opt in ("-f", "--fasta"):
            fasta = arg
        elif opt in ("-o", "--outputdir"):
            output_dir = arg
        elif opt in ("-a", "--accession"):
            accession = arg
        elif opt in ("-t", "--tree"):
            input_tree = arg
        elif opt in "-c":
            do_collapse = True

    # if we are given neither a fasta or target accession, raise an error
    if fasta == '' and accession == '' and input_tree == '':
        raise ValueError("If fasta is not supplied, accession must be provided. Use -a or --accession to specify.")

    # if we have no output dir, raise an error
    if output_dir == '':
        raise ValueError("No output directory specified. Use -o or --outputdir to specify.")

    # if we are not provided with a fasta but are provided a tree, raise an error
    if fasta == '' and input_tree != '':
        raise ValueError("A fasta file must be provided using -f or --fasta when a tree file is provided.")

    # we can run without a target accession, but the tree will not be rooted
    if accession == '':
        print("No target accession specified; tree will be unrooted. Use -a or --accession to specify.")

    # make output directory
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Could not make directory {}. Using existing directory.".format(output_dir))
    else:
        print("Created working directory at {}".format(output_dir))

    # if we are not provided with a fasta or a tree, generate one from BLAST from target accession
    if fasta == '' and accession != '' and input_tree == '':
        print("Running BLAST with accession ID {}...".format(accession))
        run_blast(accession, output_dir)
        fasta = os.path.join(output_dir, "blast_results_filtered.fasta")

    # check for tree input
    if input_tree == '':
        # run the ete3 build on command line
        print("Running ETE3...")
        cmd_line = "ete3 build -w standard_raxml -a {} -o {} -C 4 --clearall".format(fasta, output_dir)
        os.system(cmd_line)

        # get the output directory from ete3 and pull out nw file
        # note that this is hard-coded to the current workflow (ie., raxml_default)
        generated_dir = output_dir + "/clustalo_default-none-none-raxml_default"
        output_files = os.listdir(generated_dir)
        print(output_files)
        nw_file = [f for f in output_files if ".nw" in f and ".nwx" not in f][0]
        fasta = os.path.join(generated_dir, [f for f in output_files if "used_alg.fa" in f][0])

        # create the PhyloTree using the generated nw file
        t = PhyloTree(os.path.join(generated_dir, nw_file), sp_naming_function=None)
    else:
        t = PhyloTree(input_tree, sp_naming_function=None)

    # traverse the tree and add taxonomy ID and species
    for n in t.traverse():
        if not n.name:
            continue

        try:
            n.species = get_species_name(n.name, fasta)
            if n.species == "":
                n.species = 'NA'
                genus_name = 'NA'
            else:
                genus_name = re.findall('([A-Z][a-z]+)', n.species)
                if not genus_name:
                    genus_name = 'NA'
                else:
                    genus_name = genus_name[0]
            n.add_features(genus=genus_name)

        except ValueError:
            n.add_features(genus='NA')
            n.species = 'NA'

    if do_collapse:
        global node2labels
        # noinspection PyTypeChecker
        node2labels = t.get_cached_content(store_attr="genus")
        t_collapsed = Tree(t.write(features=["species", "genus"], is_leaf_fn=collapsed_leaf, format=2))

        leaves_to_keep = []
        for leaf in t_collapsed.iter_leaves():
            if leaf.species != "Unknown":
                leaves_to_keep.append(leaf)

        t_collapsed.prune(leaves_to_keep, preserve_branch_length=False)
    else:
        t_collapsed = t

    # write a CSV with node names and distances
    with open(output_dir + "/node_dist.csv", 'w') as csvfile:

        # get root node for reference
        root = t_collapsed.get_tree_root()

        # create writer, write header
        csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        header = ["node_name", "genus", "species", "dist_from_root"]
        csvwriter.writerow(header)

        # walk through each node; print name and dist to csv file
        for n in t_collapsed.traverse():
            if n.name == "":
                continue

            node_name = n.name
            node_species = n.species
            try:
                node_genus = n.genus
            except AttributeError:
                node_genus = "NA"
            root_dist = n.get_distance(root)
            row = [node_name, node_genus, node_species, root_dist]
            csvwriter.writerow(row)

    # create an annotated tree by renaming the nodes to include species
    t_annotated = copy.copy(t_collapsed)

    if do_collapse:
        for leaf in t_annotated.iter_leaves():
            leaf.__setattr__("name", leaf.genus)
    else:
        for leaf in t_annotated.iter_leaves():
            leaf.__setattr__("name", leaf.species)

    print("Creating output PNG images...")

    # create a circular tree PNG
    draw_circular_tree(t_collapsed, output_dir + "/output_circular_annotated.png")


if __name__ == "__main__":
    # main execution of the script
    main(sys.argv[1:])
