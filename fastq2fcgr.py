import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from parameters import PARAMETERS
from src.fcgr import FCGR

KMER = PARAMETERS["KMER"]
SPECIE = PARAMETERS["SPECIE"]
CLADES = PARAMETERS["CLADES"]
#CLADES = [clade+"1" for clade in CLADES]
USE_R1 = PARAMETERS["USE_R1"]
USE_R2 = PARAMETERS["USE_R2"]

FOLDER_FCGR = Path(PARAMETERS["FOLDER_FCGR"])
FOLDER_FASTQ = Path(PARAMETERS["FOLDER_FASTQ"]) # contains all <clade>.fastq files

if all([USE_R1,USE_R2]):
    fastq_files = list(f for f in FOLDER_FASTQ.glob("*.fastq"))
elif USE_R1 is True: 
    fastq_files = list(f for f in FOLDER_FASTQ.glob("*_R1.fastq"))
else:
    fastq_files = list(f for f in FOLDER_FASTQ.glob("*_R2.fastq"))

# initialize FCGR
fcgr = FCGR(KMER)

def preprocess_seq(seq):
    new_seq = []
    for c in seq:
        new_seq.append(c if c in {"A","C","G","T","N"} else "N") 
    return "".join(new_seq)

pbar = tqdm(total=len(fastq_files), desc="Creating FCGR from fastq files")
for fastq_file in fastq_files: 

    filename = str(fastq_file).split("/")[-1]
    # extract clade and "miseq_R1"/"miseq_R2" from the filename
    clade, read_type, extension = filename.split(".") 

    fastq = SeqIO.parse(fastq_file, format="fastq")

    # Read all the reads in the current clade
    current_seq_id = ""
    list_reads = []
    for record in fastq: 

        seq  = str(record.seq)
        name_seq = record.id.split('|')[0]

        # clean list of reads
        if current_seq_id == "":
            current_seq_id = name_seq
        elif name_seq != current_seq_id:
            # create and save FCGR
            sequence_from_reads  = "N".join(list_reads)
            sequence_from_reads = preprocess_seq(sequence_from_reads) 
            array = fcgr(sequence_from_reads)          
            seq_id = current_seq_id.replace("/","_")
            path_save = FOLDER_FCGR.joinpath(f"{SPECIE}/{clade}")
            path_save.mkdir(parents=True, exist_ok=True)
            np.save(path_save.joinpath(f"{seq_id}_{read_type[-2:]}.npy"), array)

            # update name and clean reads
            current_seq_id = name_seq
            list_reads = []

        # append the read 
        list_reads.append(seq)
    ## create and save FCGR the last fatstq for this series of fastq_files
    sequence_from_reads  = "N".join(list_reads)
    sequence_from_reads = preprocess_seq(sequence_from_reads) 
    array = fcgr(sequence_from_reads)          
    seq_id = current_seq_id.replace("/","_")
    path_save = FOLDER_FCGR.joinpath(f"{SPECIE}/{clade}")
    path_save.mkdir(parents=True, exist_ok=True)
    np.save(path_save.joinpath(f"{seq_id}_{read_type[-2:]}.npy"), array)

    pbar.update(1)
pbar.close()
    
