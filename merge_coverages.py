from time import time
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from genbank_coverage import make_cov_list, plot_cov



def merges_coverage(input_dir, out_dir, path_alignment, title):
    '''
    Reads files with blast results for different reference sequences from input_dir,
    for each query sequence finds the blast table with maximal hit length,
    merges strings with maximal hit lengths from different blast tables, 
    calculates final coverage, draws a coverage plot, saves coverage values to output dir.
    
    Input:
        input_dir - str - directory with blast results
        out_dir - str - directory to save output files
        path_alignment - str - path to the alignment of sequences which were
                        used as references for blast searches
        title - title for coverage plot
    Output
        Plots coverage and saves file with coverage values in output directory
    '''

    #input_dir = os.getcwd()
    blast_outputs_names = os.listdir(input_dir)

    # dictionary for  blast output dataframes

    blast_out_dict = {}

    # getting the other blast output dataframes with the rows distinct from the first df
    for name in blast_outputs_names:
        blast_out_dict[name] = pd.read_csv(Path(input_dir,name), sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore']) #.transpose()
        print("File name:{}, number of rows:{}".format(name,len(blast_out_dict[name])))

    # list with ids of all blast hits from all files
    seq_ids_all = []

    # Adds new columns to each dataframe - 'sstart_al' and 'send_al'
    # with positions in sseq according to the alignment of all reference sequence
    t1 = time()
    for name in blast_outputs_names:
        # adds all qseqid ids to list
        [seq_ids_all.append(id) for id in list(blast_out_dict[name]['qseqid'])]
    t2 = time()
    print('Retrieving ids from table {:.4}'.format(t2-t1))
    # unique ids in all files with blast output
    t1 = time()
    seq_ids_all = np.unique((np.array(seq_ids_all)))
    t2 = time()
    print('Finished. Time: {:.4}'.format(t2-t1))
    
    print('Unique ids {}'.format(len(seq_ids_all)))


    print('Creating tabe with hit lengths')
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
    print('Finished. Time: {:.4}'.format(t2-t1))
    print(length_df.head())

    t1 = time()
    # series where value is the name of blast ouput table where length of vlast hit is maximal
    # if there are two tables with maximal length the first one is chosen

    def max_length(col):
        return col[col==max(col)].index[0]

    length_df1 = length_df.apply(lambda col: col[col==max(col)].index[0], axis=0)
    #length_df = length_df.apply(max_length, axis=0)
    #length_df.to_csv(Path(input_dir,'lengths.txt'))
    t2 = time()
    print('Time for lambda creating Series with tables\' names {:.4}'.format(t2-t1))

    t1 = time()
    length_df = length_df.apply(max_length, axis=0)
    t2 = time()
    print('Time for creating series with tables\' names {:.4} '.format(t2-t1))

    print(length_df)
    
    t1 = time()
    # new table with results from several blast runs
    # for each sequence has string from the blast table where query length was the highest
    print("Merging dataframes")
    t1 = time()
    blast_table_new = pd.DataFrame(columns = blast_out_dict[name].columns)
    for i in range(len(length_df)):
        seq_id = length_df.index[i]
        table = blast_out_dict[length_df[i]]
        blast_table_new = blast_table_new._append(table[table['qseqid']==seq_id], ignore_index=True)
    t2 = time()
    print('Finished. Time: {:.4} '.format(t2-t1))
    print(blast_table_new)
    
    # alignment of reference sequences
    records_temp = SeqIO.to_dict(SeqIO.parse(open(path_alignment), "fasta"))

    # dictionary with lists of relative positions of ref sequences in the alignment
    global rel_pos_l_dict
    rel_pos_l_dict = {}

    # calculating relative postitions of reference sequences
    for name in blast_outputs_names:
        # sequence blast was run against
        seq_name = name.strip(".out").strip("blast_")
        seq = str(records_temp[seq_name].seq)
        rel_pos_l_dict[seq_name] = make_pos_list(seq)

    # calculates hit's start position in alignment of references
    def get_sstart_pos_in_al(row):
        r = rel_pos_l_dict[row['sseqid']].index(row['sstart'])
        return rel_pos_l_dict[row['sseqid']].index(row['sstart'])
    # calculates hit's end position in alignment of references
    def get_send_pos_in_al(row):
        r = rel_pos_l_dict[row['sseqid']].index(row['send'])
        return rel_pos_l_dict[row['sseqid']].index(row['send'])

    print("Adding the rows with new positions")
    t1 = time()
    blast_table_new['sstart_al'] = blast_table_new.apply(get_sstart_pos_in_al, axis=1)
    blast_table_new['send_al'] = blast_table_new.apply(get_send_pos_in_al, axis=1)
    t2 = time()
    print("Finished. Time: {:.4}".format(t2-t1))

    print(blast_table_new.head())
    blast_table_new.to_csv(Path(out_dir,'blast_new.txt'))


    reference_length = len(records_temp[list(records_temp)[0]].seq)

    print('Making final coverage')
    t1 = time()
    final_coverage = make_cov_list(blast_table_new, reference_length, 1)
    t2 = time()
    print('Finished. TIme: {:.4}'.format(t2-t1))
    plot_cov(final_coverage, out_dir, title+'_plot')

    '''
    print('Drawing final histogram')
    t1 = time()
    #plt.hist(range(1,len(final_coverage_nogap)+1,1), bins= range(1,len(final_coverage_nogap)+1,1), weights = final_coverage_nogap)
    plt.hist(range(1,len(final_coverage)+1,1), bins= range(1,len(final_coverage)+1,1), weights = final_coverage)
    plt.xlabel("Position in genome, nt")
    plt.ylabel("Number of sequences in GenBank")
    plt.xlim(0,len(final_coverage))
    plt.savefig(out_dir+title+".png")
    plt.savefig(out_dir+title+".svg")
    t2= time()
    print('Finished {}'.format(str(t2-t1)))
    plt.show()
    plt.clf()
    '''
    # saves coverage to a text file
    final_coverage_s = [str(x) for x in final_coverage]
    #final_coverage_nogap_s = [str(x) for x in final_coverage_nogap]
    with open(os.path.join(out_dir, title+"_cov.txt"),"w") as out_file:
        out_file.write(','.join(final_coverage_s))
        #out_file.write(','.join(final_coverage_nogap_s))
    out_file.close()

# This function is not used any more
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
        record_seq - str - nucleotide sequence from alignment
    Output:
        pos_list - list - list looks like like [0,0,1,2,0,0,3,4,5], 
    where zero-elements correspond to positions with gaps, and the other
    elements are serial numbers of nucleotides in record_seq
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
    args = parser.parse_args()

    merges_coverage(args.input_dir, args.output_dir, args.alignment, args.title)