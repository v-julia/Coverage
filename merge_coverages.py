import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathlib import Path
from time import time
from Bio import SeqIO

def merges_coverage(input_dir, out_dir, path_alignment, title, first_name):
    '''
    Input:
        input_dir - str - directory with blast results
        out_dir - str - directory to save output files
        path_alignment - str - path to the alignment of sequences which were
                        used as references for blast searches
        title - 
    '''

    #input_dir = os.getcwd()
    blast_outputs_names = os.listdir(input_dir)

    # dictionary for complete first blast output dataframe and the other dataframes
    # with rows absent in the first dataframes
    blast_out_dict = {}

    # reading the first df
    blast_out_dict[first_name] = pd.read_csv(Path(input_dir,first_name), sep='\t', header = None, \
                                    names=['qseqid','sseqid','pident','length','mismatch',\
                                    'gapopen','qstart','qend','sstart','send','evalue','bitscore']) #.transpose()



    # getting the other blast output dataframes with the rows distinct from the first df
    for name in blast_outputs_names:
        if name == first_name:
            continue
        else:
            print(name)
            blast_out_dict[name] = compare_blast_out(Path(input_dir,first_name),Path(input_dir,name))
            print(len(blast_out_dict[name]))
    # print(blast_out_dict["blast_AY593810.1.out"])

    # alignment of reference sequences
    records_temp = SeqIO.to_dict(SeqIO.parse(open(path_alignment), "fasta"))

    # dict_relative[seq_id] = [] list with relative positions of nucleotides in sequences 
    # of seq_id in the alignment of reference sequences
    # seq_id - id of reference sequence in records_temp
    dict_relative = {}
    # Dictionary with coverages according to the relative positions of sequences in the alignments
    dict_coverage = {}

    reference_length = len(records_temp[list(records_temp)[0]].seq)
    final_coverage = [0]*reference_length
    print("Total length of alignment {}".format(str(reference_length)))

    for proto_id in blast_outputs_names:
        print(proto_id)
        # reference sequence with gaps
        seq = str(records_temp[proto_id.strip(".out").strip("blast_")].seq)
        seq_ng_length = len(seq.replace('-',''))
        rel_pos_list = make_pos_list(seq)
        #print(rel_pos_list)
        coverage_list = make_cov_list(blast_out_dict[proto_id], seq_ng_length)
        #print(len(coverage_list))
        #print(coverage_list)
        temp_rel_cov = [0]*len(records_temp[list(records_temp)[0]].seq)
        for i in range(1,seq_ng_length+1):
            temp_rel_cov[rel_pos_list.index(i)] = coverage_list[i-1]
        final_coverage = np.add(final_coverage, temp_rel_cov)

    final_coverage_nogap = []
    for i in range(len(final_coverage)):
        if i > 100 and i < 8000:
            if final_coverage[i] < 200:
                continue
            else:
                final_coverage_nogap.append(final_coverage[i])
        '''
        else:
            if i < 3100 and i > 2900:
                if final_coverage[i]<200:
                    continue
                else:
                    final_coverage_nogap.append(final_coverage[i])
            else:
                    final_coverage_nogap.append(final_coverage[i])
        '''
        '''
        if i > 100 and i < 3100:
            if final_coverage[i] < 50:
                continue
            else:
                final_coverage_nogap.append(final_coverage[i])
        else:
            if i < 3950 and i > 3400:
                if final_coverage[i]<1300:
                    continue
                else:
                    final_coverage_nogap.append(final_coverage[i])
            else:
                    final_coverage_nogap.append(final_coverage[i])
        '''
    print(len(final_coverage))
    #plt.hist(range(1,len(final_coverage_nogap)+1,1), bins= range(1,len(final_coverage_nogap)+1,1), weights = final_coverage_nogap)
    plt.hist(range(1,len(final_coverage)+1,1), bins= range(1,len(final_coverage)+1,1), weights = final_coverage)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(final_coverage))
    plt.savefig(out_dir+title+".png")
    plt.savefig(out_dir+title+".svg")
    plt.show()
    plt.clf()
    plt.hist(range(1,len(final_coverage_nogap)+1,1), bins= range(1,len(final_coverage_nogap)+1,1), weights = final_coverage_nogap)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(final_coverage_nogap))
    plt.savefig(out_dir+title+"_nogap.png")
    plt.savefig(out_dir+title+"_nogap.svg")
    #plt.show()

    final_coverage_s = [str(x) for x in final_coverage]
    final_coverage_nogap_s = [str(x) for x in final_coverage_nogap]
    with open(os.path.join(out_dir, title+"_cov.txt"),"w") as out_file:
        out_file.write(','.join(final_coverage_s))
        out_file.write(','.join(final_coverage_nogap_s))
    out_file.close()

def compare_blast_out(blast_output1, blast_output2):
    '''
    Compares two dataframes, returns the second dataframe with the rows which 
    first column values (qseqid) are absent in the first column of the first dataframe

    Input:
        blast_output1 - str - path to file with the first dataframe
        blast_output2 - str - path to file with the second dataframe

    Output:
        blast_output2df - pandas dataframe - blast_output2 with rows which first columns
        are absent in blast_output1 dataframe
    '''

    blast_output1df = pd.read_csv(blast_output1, sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore'])
    blast_output2df = pd.read_csv(blast_output2, sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore'])

    #get the rows of blast_table2 which query names absent in blast_table1
    blast_output2df = blast_output2df[:][blast_output2df['qseqid'].isin(blast_output1df['qseqid']) != True]

    return blast_output2df #.transpose()


def make_pos_list(record_seq):
    '''
    Creates a list with length equal to reference alignment length
    for sequence record_seq. Values in positions of alignment where record_seq
    contains gap are zeros. Other values correspond to position of nucleotide in
    record_seq.
    Input:
        record_seq - str - nucleotide sequences from alignment
    Output:
        pos_list - list - list looks like like [0,0,1,2,0,0,3,4,5], 
    where zero-elements correspond to positions with gaps, and the other
    elements are serial numbers of nucleotides0
    '''
    pos_list = []
    k=1
    for i in range(len(record_seq)):
        if record_seq[i] != '-':
            pos_list.append(k)
            k+=1
        else:
            pos_list.append(0)
    return pos_list

def make_cov_list(blast_out_df, reference_length):
    '''
    Creates list with coverage of each sequence position for one blast output table
    
    Input:
        record_seq - str - nucleotide sequence
        reference_length - int - length of reference sequence blast was run against
    Output:
        pos_coverage - list - list with coverage values for each position in the sequence
        blast was run against.
    '''

    #ID of reference seq blast was performed against
    ref_id = blast_out_df['sseqid'][list(blast_out_df.index)[0]]

    # list for sequence counts in each position of reference sequence
    pos_coverage = [0]*reference_length

    # adds counts to pos_coverage list according to blast hits in blast_output

    def add_cov(row):
        for i in range(row["sstart"]-1, row["send"], 1):
            #pos_coverage[rel_pos_list.index(i)] +=1
            pos_coverage[i] +=1

    blast_out_df.apply(add_cov, axis=1)

    return pos_coverage



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str,
                        help="Input directory with blast output tables", required=True)
    parser.add_argument("-o", "--output_dir", type=str,
                        help="Output directory", required=True)
    parser.add_argument("-al", "--alignment", type=str,
                        help="Path to file with alignment of reference sequences", required=True)
    parser.add_argument("-t", "--title", type=str,
                        help="Title of output file figure", required=True)
    parser.add_argument("-f", "--first", type=str,
                        help="The first blast output file", required=True)
    args = parser.parse_args()

    pos_coverage = merges_coverage(args.input_dir, args.output_dir, args.alignment, args.title, args.first)