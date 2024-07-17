import argparse
import os
import sys
import time

from Bio import SeqIO
from pathlib import Path

from genbank_coverage import fetch_from_GB, find_coverage
from merge_coverages import merges_coverage

"""
For each sequences in input_file (or sequences loaded from GenBank nucleotide using query) 
performs standalone blast against reference sequences from reference, saves blast output table in 'blast_out' folder.
Then merges results using merge_coverage.py. Draws a coverage plot and saves coverage values to the directory of input file.
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str,
                        help="File with sequences in fasta or genbank format for calculating coverage")
    parser.add_argument("-q", "--query", type=str,
                        help="Query for GenBank Nucleotide database")
    parser.add_argument("-ref", "--reference", type=str,
                        help="List with reference sequences'ids (with version!) separated by comma. \
                        Example: AB084913.1,AB0834732.2", required= True)
    parser.add_argument("-f", "--figname", type=str,
                        help="Name for final figure with coverage (without extension)", required= True)
    parser.add_argument("-path_blast", "--path_blast", type=str,
                        help="Path to blast program",  required = True)
    parser.add_argument("-path_mafft", "--path_mafft", type=str,
                        help="Path to mafft program",  required = True)
    args = parser.parse_args()

    if not (args.input_file or args.query):
        print("File with sequences or query for Nucleotide database are not defined. \
        Please, use genbank_coverage.py --help")
    else:
        # reads sequences from file if it exists
        if args.input_file != None:
            # ADD EXCEPTION ERROR IF INPUT FILE DOES NOT EXIST
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

        # list with names of files with sequences
        ref_files_list = []
        # list to reference records
        ref_records_list = []
        
        #list with reference ids
        args.reference = args.reference.strip('\'').split(",")
        args.reference = [x.strip() for x in args.reference]


        # searches reference sequence in file with sequences downloaded from GenBank
        print("Checking the presence of reference sequences in loaded entries/file")
        
        
        for ref_seq_id in args.reference:
            print(ref_seq_id)
            time.sleep(2)
            ref_f_name = str(Path(out_directory,ref_seq_id+".fasta"))
            # First check whether file with reference in present in current directory
            if os.path.exists(ref_f_name):
                ref_seq = list(SeqIO.parse(ref_f_name, "fasta"))[0]
                ref_records_list.append(ref_seq)
            # Then check whether reference sequences are present in input file
            else:
                print("Searching reference sequences in input file")
                time.sleep(2) 
                # all records loaded from GenBank
                records = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))
                if ref_seq_id in records.keys():
                     ref_seq = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))[ref_seq_id]
                     ref_records_list.append(ref_seq)
                     ref_f_name = os.path.join(out_directory,ref_seq_id+'.fasta')
                     SeqIO.write(ref_seq, ref_f_name, "fasta")
                # Otherwise download from GenBank Nucleotide
                else:
                    print("Loading reference sequences")
                    time.sleep(2) 
                    ref_f_name = fetch_from_GB(ref_seq_id, ref_f_name)
                    ref_seq = list(SeqIO.parse(ref_f_name, "fasta"))[0]
                    ref_records_list.append(ref_seq)
                    
            ref_files_list.append(ref_f_name)

        # writes reference sequence to one file
        reference_file_n = str(Path(os.path.split(args.input_file)[0],"reference.fasta"))
        SeqIO.write(ref_records_list, reference_file_n, "fasta")
        print("Finished. Reference sequences were saved in {}".format(reference_file_n))
        
        t0 = time.time()
        # aligns reference sequences
        print("Aligning reference sequences")
        reference_aln_file_n = str(Path(os.path.split(args.input_file)[0],"reference_aln.fasta"))
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            os.system(('{} --op 15 --ep 3 --retree 1 ' + reference_file_n+ ' > ' + \
                   reference_aln_file_n).format(str(Path(args.path_mafft,'mafft.bat'))))
        else:
            os.system(('{} --op 15 --ep 3 --retree 1 ' + reference_file_n+ ' > ' + \
                   reference_aln_file_n).format(str(Path(args.path_mafft,'mafft'))))

        print('Calculating coverage using reference sequences')

        input_dir = os.path.split(args.input_file)[0]
        for ref_seq_id, ref_seq_file in zip(args.reference,ref_files_list):
            # founds coverage of reference sequence by sequences downloaded from GenBank
            pos_coverage = find_coverage(args.input_file, ref_seq_file, args.path_blast)
            # writes coverage to the text file
            fout_cov = open(input_dir+ref_seq_id+"_cov.txt", "w")
            fout_cov.write(",".join(str(x) for x in pos_coverage))
            fout_cov.close()

        out_dir = str(Path(input_dir, "blast_out"))
        print('Merging blast results from different reference sequences')
        t1 = time.time()
        merges_coverage(out_dir, input_dir, reference_aln_file_n, args.figname)
        t2 = time.time()
        print("Finished merging {:.4}".format(t2-t1))

        tn = time.time()
        print('Total time {:.4}'.format(tn-t0))