"""After undersample the sequences, we extract each one in a separated fasta file
to generates the FCGR as .npy files for training purposes"""
from tqdm import tqdm 
from Bio import SeqIO
from pathlib import Path
import pandas as pd
from collections import defaultdict
from parameters import PARAMETERS

print(">> Extracting sequences <<")

PATH_FASTA_GISAID = PARAMETERS["PATH_FASTA_GISAID"]
SPECIE = PARAMETERS["SPECIE"]
FOLDER_FASTA = Path(PARAMETERS["FOLDER_FASTA"]) # here will be saved the selected sequences by clade
FOLDER_FASTA.mkdir(parents=True, exist_ok=True)

# Create folder for each clade
# for clade in PARAMETERS["CLADES"]: 
#     FOLDER_FASTA.joinpath(clade).mkdir(parents=True, exist_ok=True)

# load fasta_id to save
undersample = pd.read_csv("data/train/undersample_by_clade.csv").to_dict("records")
set_fasta_id = set([record.get("fasta_id") for record in undersample])
clades_by_fastaid = {record.get("fasta_id"): record.get("clade") for record in undersample} 

# collect seqrecords by clade
seqrec_by_clade = defaultdict(list)

pbar = tqdm(total=len(set_fasta_id))
# Read fasta with all sequences from GISAID
with open(PATH_FASTA_GISAID) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        
        # save sequence if it was selected
        if record.description in set_fasta_id:
            # save sequence in a fasta file "<accession_id>.fasta"
            clade    = clades_by_fastaid.get(record.description) 

            # append the sequence found into the dictionary in his respective clade list
            seqrec_by_clade[clade].append(record)

            #filename = record.description.replace("/","_") # replace '/' to avoid problems when saving fasta file
            #path_save = FOLDER_FASTA.joinpath(f"{clade}/{filename}.fasta")
            # if not path_save.is_file():
            #     SeqIO.write(record, path_save, "fasta") 
            # remove from the set to be saved   
            set_fasta_id.remove(record.description)
            pbar.update(1)

        # if all sequences has been saved, break the loop
        if not set_fasta_id:
            break

pbar.close()

# save all the sequences per clade in one fasta per clade
for clade, list_seqs in seqrec_by_clade.items():

    path_save = FOLDER_FASTA.joinpath(f"{clade}.fa")
    SeqIO.write(list_seqs, path_save, "fasta") 

if len(set_fasta_id):
    pd.Series(list(set_fasta_id)).to_csv(FOLDER_FASTA.joinpath("not_found.csv"))