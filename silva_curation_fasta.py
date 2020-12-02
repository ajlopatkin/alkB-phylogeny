from Bio import Entrez
import re
import csv
import random
import os

input_dir = "/Users/ShoshanaWilliams/Desktop/silva/"
fasta_input = "fused_rds_fastas.fasta"
silva_type_strains = "silva_tstrains_20200601.fasta"
silva_taxID = "silva_taxid.txt"

# set up the entrez parameters
api_key = "194183db70139b8a782130aef1cdfee6cb09"
Entrez.api_key = api_key
Entrez.max_tries = 1
Entrez.email = "scw2145@barnard.edu"

def get_tax_id(acc_id):
    """
    Get the taxonomy ID for a given accession ID using Entrez.
    :param acc_id: the accession ID to be converted to taxonomy ID.
    :param api_key: the api key to use for access to entrez.
    :return: the taxonomy ID for the given accession ID- this is an integer.
    """
    # define and execute the query on Entrez
    handle = Entrez.esummary(db="protein", id=acc_id, retmode="xml")
    records = Entrez.parse(handle)

    # pull out the TaxID from the record and return
    try:
        o = [rec['TaxId'] for rec in records]
        return int(o[0])
    except RuntimeError:
        return 0

def get_species_name(acc_id, fasta_file):
    """
    Get a species name given an by searching the fasta input.
    :param acc_id: the accession ID to convert.
    :param fasta_file: the fasta to search.
    :return: the species name- this is a string.
    """

    with open(fasta_file) as file:
        for line in file:
            if acc_id in line:
                find = line
                break

    return find.split('[')[1].split(']')[0]

# taxID to silva
taxID_filename = input_dir + silva_taxID
tax_id = []
acc_id = []
with open(taxID_filename, 'r') as open_file:
    for line in open_file:
        if ";" in line:
            tax_id.append(line.strip("\n").split(";")[1])
            acc_id.append(line.strip("\n").split(";")[0])
dictionary_accIDtoTaxID = dict(zip(acc_id, tax_id))

# get species name
input_filename = input_dir + fasta_input
species_name_RDS = []
fasta_RDS = []
with open(input_filename, 'r') as open_file:
    for line in open_file:
        if line[0] == ">":
            if "fatty acid desaturase" in line:
                continue
            species_name = ''
            try:
                species_name = line.split('[')[1].split(']')[0]
            except:
                pass
            species_name_RDS.append(species_name)
        #also keep track of the fastas
        else:
            fasta_RDS.append(line.strip("\n"))
dictionary_speciesName_to_fasta = dict(zip(species_name_RDS, fasta_RDS))            

# determine species name and acc ID from all type strains (SILVA)
species_name_Silva = []
accID_Silva = []
silva_filename = input_dir + silva_type_strains
with open(silva_filename, 'r') as open_file:
    for line in open_file:
        if line[0] == ">":
            accID_Silva.append(line.strip("\n").split(" ")[0][1:])
            species_name_Silva.append(re.search("[^;]+$", line.strip("\n"))[0])
dictionary_speciesName_to_accID = dict(zip(species_name_Silva, accID_Silva))

# collect good species and corresponding accession ID's, taxIDs, and FASTAS
curated_species = []
curated_acc = []
curated_taxID = []
curated_fasta = []
for species in species_name_RDS:
    if species in species_name_Silva:
        if species not in curated_species:
            curated_species.append(species)
            curated_acc.append(dictionary_speciesName_to_accID.get(species))
            curated_taxID.append(dictionary_accIDtoTaxID.get(curated_acc[-1].split(".")[0]))
            curated_fasta.append(dictionary_speciesName_to_fasta.get(species))

write_path = open("silva_results.fasta", "w")
for species in curated_species:
    fasta = dictionary_speciesName_to_fasta.get(species)
    acc = dictionary_speciesName_to_accID.get(species)
    write_path.write(">" + acc+ " alk monooxygenase" + " [" + species + "] \n" + fasta + "\n")
write_path.close()
