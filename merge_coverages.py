from time import time

t1 = time()
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathlib import Path

from Bio import SeqIO

t2 = time()

print(t2-t1)

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
            blast_out_dict[name] = pd.read_csv(Path(input_dir,name), sep='\t', header = None, \
                                    names=['qseqid','sseqid','pident','length','mismatch',\
                                    'gapopen','qstart','qend','sstart','send','evalue','bitscore']) #.transpose()
            #blast_out_dict[name] = compare_blast_out(Path(input_dir,first_name),Path(input_dir,name))
            print(len(blast_out_dict[name]))
    # print(blast_out_dict["blast_AY593810.1.out"])


    # alignment of reference sequences
    records_temp = SeqIO.to_dict(SeqIO.parse(open(path_alignment), "fasta"))




    # list with positions of sequence's nucleotides in alignment 
    # of reference sequences
    global rel_pos_list

    # list with ids of all blast hits from all files
    seq_ids_all = []

    # Adds new columns to each dataframe - 'sstart_al' and 'send_al'
    # with positions in sseq according to the alignment of all reference sequence
    for name in blast_outputs_names:
        '''
        # sequence blast was run against
        seq = str(records_temp[name.strip(".out").strip("blast_")].seq)
        # length of sequence without gaps
        seq_ng_length = len(seq.replace('-',''))

        rel_pos_list = make_pos_list(seq)

        blast_out_dict[name]['sstart_al'] = blast_out_dict[name].apply(get_sstart_pos_in_al, axis=1)
        blast_out_dict[name]['send_al'] = blast_out_dict[name].apply(get_send_pos_in_al, axis=1)
        print(blast_out_dict[name].head())
        '''
        # adds all qseqid ids to list
        [seq_ids_all.append(id) for id in list(blast_out_dict[name]['qseqid'])]

    # unique ids in all files with blast output
    seq_ids_all = np.unique((np.array(seq_ids_all)))
    print(len(seq_ids_all))

    # dictionary with length of coverage for sequences from blast output
    # length_dict[seq_id][blast_out_name] = length
    # seq_id - qseq, blast_out_name - name of blast table,
    # length - sum of 'length' columns in table blast_out_name for sequence with id seq_id
    length_dict = {}
    t1 = time()
    for seq_id in seq_ids_all:
        length_dict[seq_id] = {}
        for name, df in blast_out_dict.items():
            length_dict[seq_id][name] = sum(df[df['qseqid'] == seq_id]['length'])

    # dataframe
    # column names - ids of sequences
    # indices - names of blast output tabels
    # values - length of blast hit seq_id in blast output table
    length_df = pd.DataFrame.from_dict(length_dict)
    t2 = time()
    print(t2-t1)
    print(length_df.head())

    # series where value is the name of blast ouput table where length of vlast hit is maximal
    # if there are two tables with maximal length the first one is chosen
    length_df = length_df.apply(lambda col: col[col==max(col)].index[0], axis=0)
    print(length_df.head())
    # new table with results from several blast runs
    # for each sequence has string from the blast table where query length was the highest
    blast_table_new = length_df.apply(lambda x: blast_out_dict[str(x[0])][blast_out_dict[str(x[0])]==x.name],axis=0)
    print(blast_table_new.head())
    # sorted by sseqid
    blast_table_new = blast_table_new.sort_values(by='sseqid')

    # dictionary with lists of relative positions of ref sequences in the alignment
    rel_pos_l_dict = {}
    print(blast_outputs_names)
    for name in blast_outputs_names:
        # sequence blast was run against
        seq_name = name.strip(".out").strip("blast_")
        seq = str(records_temp[seq_name].seq)
        rel_pos_l_dict[seq_name] = make_pos_list(seq)

    def get_sstart_pos_in_al(row):
        r = rel_pos_l_dict[row['sseqid']].index(row['sstart'])
        return rel_pos_l_dict[row['sseqid']].index(row['sstart'])

    def get_send_pos_in_al(row):
        r = rel_pos_l_dict[row['sseqid']].index(row['send'])
        return rel_pos_l_dict[row['sseqid']].index(row['send'])
    blast_table_new['sstart_al'].apply(get_sstart_pos_in_al, axis=0)
    blast_table_new['send_al'].apply(get_send_pos_in_al, axis=0)

    print(blast_table_new.head())



    reference_length = len(records_temp[list(records_temp)[0]].seq)

    final_coverage = make_cov_list(blast_out_dict[proto_id], reference_length)

    #plt.hist(range(1,len(final_coverage_nogap)+1,1), bins= range(1,len(final_coverage_nogap)+1,1), weights = final_coverage_nogap)
    plt.hist(range(1,len(final_coverage)+1,1), bins= range(1,len(final_coverage)+1,1), weights = final_coverage)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(final_coverage))
    plt.savefig(out_dir+title+".png")
    plt.savefig(out_dir+title+".svg")
    plt.show()
    plt.clf()
    '''
    plt.hist(range(1,len(final_coverage_nogap)+1,1), bins= range(1,len(final_coverage_nogap)+1,1), weights = final_coverage_nogap)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(final_coverage_nogap))
    plt.savefig(out_dir+title+"_nogap.png")
    plt.savefig(out_dir+title+"_nogap.svg")
    #plt.show()
    '''
    final_coverage_s = [str(x) for x in final_coverage]
    #final_coverage_nogap_s = [str(x) for x in final_coverage_nogap]
    with open(os.path.join(out_dir, title+"_cov.txt"),"w") as out_file:
        out_file.write(','.join(final_coverage_s))
        #out_file.write(','.join(final_coverage_nogap_s))
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