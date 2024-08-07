import argparse
import matplotlib.pyplot as plt
import os
import subprocess
import sys
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from pathlib import Path
from time import time


def fetch_from_GB(query, filename="file.fasta"):
    '''
    Downloads sequences found using query from GenBank
    Writes them to file in fasta-format
    Input:
        query - str - query for search in GenBank
    Output:
        f_name - str - name of ouput file with sequences in fasta-format
    '''
    Entrez.email = "A.N.Other@example.com"
    # list with ids obtained by query
    record = Entrez.read(Entrez.esearch(db="nucleotide", term=query, idtype="acc", RetMax=1000000))
    id_list = record['IdList']
    print("Query to GenBank Nucleotide database:\"{}\"".format(query))
    print("Number of records found: {}".format(record["Count"]))
    print("{} records will be downloaded".format(len(id_list)))
    # output file with sequences will be saved in current working dir
    f_name = Path(os.getcwd(),filename)

    # output file wuth sequences in fasta format
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
    print('Finished')
    return(f_name)

def make_cov_list(blast_out_df, reference_length, merged=0):
    '''
    Creates list with coverage of each sequence position for one blast output table
    
    Input:
        blast_out_df - pandas dataframe - dataframe with blast results
        reference_length - int - length of reference sequence blast was run against
        merged -  int - 1 for making coverage from merged blast outputs, 0 for single blast output table
    Output:
        pos_coverage - list - list with coverage values for each position in the sequence
        blast was run against.
    '''

    # list for sequence counts in each position of reference sequence
    pos_coverage = [0]*reference_length

    # adds counts to pos_coverage list according to blast hits in blast_output
    if merged == 0:
        def add_cov(row):
            for i in range(row["sstart"]-1, row["send"], 1):
                #pos_coverage[rel_pos_list.index(i)] +=1
                pos_coverage[i] +=1
    else:
        def add_cov(row):
            for i in range(row["sstart_al"]-1, row["send_al"], 1):
                #pos_coverage[rel_pos_list.index(i)] +=1
                pos_coverage[i] +=1

    blast_out_df.apply(add_cov, axis=1)

    return pos_coverage

def plot_cov(pos_cov_list, out_dir, title):
    '''
    Plots coverage of reference sequence
    Input:
        pos_cov_list - list - list with coverage count for each position in
                                reference seq
        out_dir - str - path to output directory
        title - str - title of plot
    Output:
        saves coberage plot in .png and .svg formats in output directory
    '''
    #plt.hist(range(1,len(pos_coverage)+1,1), bins= range(1,len(pos_coverage)+1,1), weights = pos_coverage)
    plt.plot(pos_cov_list)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(pos_cov_list))
    plt.savefig(str(Path(out_dir, title + "_cov.svg")))
    print("Output figure saved as {}".format(str(Path(out_dir, title + "_cov.svg"))))
    plt.savefig(str(Path(out_dir, title + "_cov.png")))
    plt.clf()


def find_coverage(input_file, reference, path_to_blast):
    '''
    Performs standalone blast of sequences from input_file against reference sequence,
    saves blast output table in 'blast_out' folder.
    Checks whether sequence from input_file overlaps with reference sequence, 
    calculates coverage for each position of reference, plots coverage.
    
    Input:
        input_file - name of file with sequences in fasta format
        reference - name of file with reference sequence, fasta format
        path_to_blast - Path to blast program
    Output:
        pos_coverage - list - list with coverage counts in each position of reference sequence
    '''
    
    # path to ouput dorectory
    output_dir = os.path.split(input_file)[0]

    # id of reference sequence
    ref_id = os.path.splitext(os.path.split(reference)[1])[0]

    # creates "blast_out" directory if it doesn't exist
    if not os.path.exists(Path(output_dir, "blast_out")):
        try:
            print(Path(output_dir, "blast_out"))
            os.system('mkdir ' + str(Path(output_dir, "blast_out")))
        except:
            print('Can\'t create output folder')

    # name of output blast file
    blast_table = str(Path(output_dir, "blast_out", "blast_" + ref_id + ".out"))
    # name of local database
    local_db = str(Path(output_dir, "local_db_"+ref_id))

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


    print("Performing blast search")

    #creating local database using reference sequence
    subprocess.call(makeblast_command, shell=True)

    #blast against reference sequences
    subprocess.call(blastn_command, shell=True)
    #print(blast_table)

    print('Finished')
    
    if sys.platform == "win32" or sys.platform == "cygwin":
        remove_db_command = "del {}*".format(local_db)
    else:
        remove_db_command = "rm {}*".format(local_db)
    
    try:
        subprocess.call(remove_db_command, shell=True)
    except:
        print("Could not remove blast db")
    
    #dataframe with blast results
    blast_out_df = pd.read_csv(blast_table, sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore'])

    # length of reference sequence
    reference_length = len(list(SeqIO.parse(reference, "fasta"))[0].seq) 
    print("Reference sequence ID: {}".format(list(SeqIO.parse(reference, "fasta"))[0].id))

    #list for sequence counts in each position
    print('Calculating coverage using blast results')
    t1 = time()
    pos_coverage = make_cov_list(blast_out_df,reference_length, 0)
    '''
    global pos_coverage
    pos_coverage = [0]*reference_length

    print('Calculating coverage using blast results')
    t1 = time()
    def add_cov(row):
        for i in range(row["sstart"]-1, row["send"], 1):
            pos_coverage[i] +=1


    blast_out_df.apply(add_cov, axis=1)
    '''
    t2 = time()
    print('Finished\n Time: {:.4}'.format(t2-t1))

    print('Drawing coverage plot')
    t1 = time()
    # title for pictures
    pic_title = os.path.splitext(os.path.split(input_file)[1])[0] + '_' + ref_id
    plot_cov(pos_coverage, os.path.split(input_file)[0], pic_title)
    t2 = time()
    print('Finished\n Time {:.4}'.format(t2-t1))

#    plt.show()
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
        # reads input file
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
            records = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))
            if args.reference in records.keys():
                ref_seq = records[args.reference]
                ref_f_name = os.path.join(out_directory,args.reference+'.fasta')
                SeqIO.write(ref_seq, ref_f_name, "fasta")
            else:
                ref_f_name = fetch_from_GB(args.reference)
            args.reference = ref_f_name
        #finds coverage of reference sequence by sequences downloaded from GenBank
        # returns list with sequence counts for each position in sequence
        pos_coverage = find_coverage(args.input_file, args.reference, args.path_blast)

        # saves list with coverage to a text file
        fout_cov = open(os.path.splitext(args.input_file)[0]+"_"+args.reference.strip('fasta')+"_cov.txt", "w")
        fout_cov.write(",".join(str(x) for x in pos_coverage))
        fout_cov.close()
