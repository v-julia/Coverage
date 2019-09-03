import argparse
import matplotlib.pyplot as plt
import os
import subprocess
import sys
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from pathlib import Path


def fetch_from_GB(query):
    '''
    Downloads sequences found using query from GenBank
    Writes them to file in fasta-format
    Input:
        query - str - query for search in GenBank
    Output:
        f_name - str - name of ouput file with sequences in fasta-format
    '''
    Entrez.email = "vjulia94@gmail.com"
    # list with ids obtained by query
    id_list = Entrez.read(Entrez.esearch(db="nucleotide", term=query, idtype="acc"))['IdList']

    # output file with sequences will be saved in current working dir
    f_name = Path(os.getcwd(),"file.fasta")

    file_out = open(f_name, 'w')
    # Saving sequences 
    i = 0
    # fetch each sequence from GenBank by its ID and writes to ouput file
    for id in id_list:
        i+=1
        if i % 1000 ==0:
            print("Downloaded {} sequences".format(i))
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        for line in handle:
            file_out.write(line)
        handle.close()
    file_out.close()
    return(f_name)

def find_coverage(input_file, reference, path_to_blast):
    '''
    Performs standalone blast of sequences from input_file against reference sequence
    Check whether sequence from input_file overlaps with reference sequence
    Input:
        input_file - name of file with sequences in fasta format
        reference - name of file with reference sequence, fasta format
        output_dir - output directory to save fasta file with sequences that overlap with reference
    Output:
        pos_coverage - list - list with coverage counts in each position of reference sequence

    '''
    output_dir = os.path.split(input_file)[0]
    # sequence objects from input file
    records = list(SeqIO.parse(input_file, "fasta"))

    # id of reference sequence
    ref_id = os.path.splitext(os.path.split(reference)[1])[0]

    seq_ordered = [] # copy of records
    for rec in records:
        seq_ordered.append(rec.id)

    # name of output blast file
    blast_table = Path(output_dir, "blast_" + ref_id + ".out")
    # name of local database
    local_db = Path(output_dir, "local_db_"+ref_id)

    if sys.platform == 'win32' or sys.platform == 'cygwin':
        makeblast_command = "{}makeblastdb.exe -in {} -dbtype nucl -out {}".format(path_to_blast,\
           reference, local_db)
        blastn_command = "{blast_path}blastn.exe -db {db} -query {input} -outfmt 6 -out {out_blast} \
                             -strand plus -evalue 1e-20 -word_size 7".format(blast_path = path_to_blast, \
                            input = input_file, db = local_db, out_blast = blast_table)
    else:
        makeblast_command = '{}makeblastdb -in {} -dbtype nucl -out {}'.format(path_to_blast, reference, local_db)
        blastn_command = '{blast_path}blastn -db {db} -query {input} -outfmt 6 -out \
                            {out_blast} -strand plus -evalue 1e-20 -word_size 7'.format(blast_path = path_to_blast, \
                            input = input_file, db = local_db,  out_blast = blast_table)


    #creating local database using reference sequence
    subprocess.call(makeblast_command, shell=True)

    #blast against reference sequences
    subprocess.call(blastn_command, shell=True)
    print(output_blast)


    #dataframe with blast results
    blast_out_df = pd.read_csv(blast_table, sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore'])

    # length of reference sequence
    reference_length = len(list(SeqIO.parse(reference, "fasta"))[0].seq) 
    print(list(SeqIO.parse(reference, "fasta"))[0])
    #list for sequence counts in each position
    pos_coverage = [0]*reference_length

    def add_cov(row):
        for i in range(row["sstart"]-1, row["send"], 1):
            pos_coverage[i] +=1
    blast_out_df.apply(add_cov, axis=1)

    plt.hist(range(1,len(pos_coverage)+1,1), bins= range(1,len(pos_coverage)+1,1), weights = pos_coverage)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(pos_coverage))
    plt.savefig(os.path.splitext(args.input_file)[0]+"_"+ref_id+"_cov.svg")
    plt.savefig(os.path.splitext(args.input_file)[0]+"_"+ref_id+"cov_.png")
    plt.show()
    return pos_coverage

#path_to_blast = 'C:/Program Files/NCBI/blast-2.9.0+/bin'

#reference = 'D:/MY_FILES/DATA/Lukashev/Astroviruses/reference_strains/HAstV-1_ref_complete.fasta'
#input_file = 'D:/MY_FILES/DATA/Lukashev/Astroviruses/MAstV-1/HAstV-1/Ann/HAst-1_allgb_wref.fasta'


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str,
                        help="File with sequences in fasta or genbank format for calculating coverage")
    parser.add_argument("-q", "--query", type=str,
                        help="Query for GenBank Nucleotide database")
    parser.add_argument("-ref", "--reference", type=str,
                        help="Path to file with reference sequence in fasta-format or GenBank \
                       accession number of reference sequence", required= True)
    parser.add_argument("-path_blast", "--path_blast", type=str,
                        help="Path to blast program",  required = True)
    args = parser.parse_args()

    if not (args.input_file or args.args.query):
        print("File with sequences or query for Nucleotide database are not defined. \
        Please, use genbank_coverage.py --help")
    else:
        if os.path.exists(args.input_file):
            out_directory = os.path.split(args.input_file)[0]
            file_ext = os.path.splitext(args.input_file)[1]
            # converts genbank to fasta
            if file_ext == ".gb" or file_ext == ".genbank":
                with open(args.input_file) as input_handle:
                    input_f_name = os.path.splitext(args.input_file)[0]+".fasta"
                    with open(input_f_name, "w") as output_handle:
                        sequences = SeqIO.parse(input_handle, "genbank")
                        count = SeqIO.write(sequences, output_handle, "fasta")
                input_handle.close()
                output_handle.close()
                args.input_file = input_f_name
        else:
            out_directory = os.getcwd()
            # downloades sequences from GenBank and saves them to fasta_file
            args.input_file = fetch_from_GB(args.query)

        # searches reference sequence in file with sequences downloaded from GenBank
        if not os.path.exists(args.reference):
            ref_seq = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))[args.reference]
            
            ref_f_name = os.path.join(out_directory,args.reference+'.fasta')
            SeqIO.write(ref_seq, ref_f_name, "fasta")
            args.reference = ref_f_name

        #founds coverage of reference sequence by sequences downloaded from GenBank
        # list with 
        pos_coverage = find_coverage(args.input_file, args.reference, args.path_blast)

        fout_cov = open(os.path.splitext(args.input_file)[0]+"_cov.txt", "w")
        fout_cov.write(",".join(str(x) for x in pos_coverage))
        fout_cov.close()





