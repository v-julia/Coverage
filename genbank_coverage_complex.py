import argparse
import os
import sys
from Bio import SeqIO
from pathlib import Path
from genbank_coverage import fetch_from_GB, find_coverage
from merge_coverages import merges_coverage


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str,
                        help="File with sequences in fasta or genbank format for calculating coverage")
    parser.add_argument("-q", "--query", type=str,
                        help="Query for GenBank Nucleotide database")
    parser.add_argument("-ref", "--reference", type=str,
                        help="List with reference sequences'ids separated by comma", required= True)
    parser.add_argument("-t", "--title", type=str,
                        help="Title for coverage plot", required= True)
    parser.add_argument("-path_blast", "--path_blast", type=str,
                        help="Path to blast program",  required = True)
    parser.add_argument("-path_mafft", "--path_mafft", type=str,
                        help="Path to mafft program",  required = True)
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

        # list with names of files with sequences
        ref_files_list = []
        # list to reference records
        ref_records_list = []
        #list with reference ids
        args.reference = args.reference.strip('\'').split(",")
        args.reference = [x.strip() for x in args.reference]


        # searches reference sequence in file with sequences downloaded from GenBank
        print("Checking the presence of reference sequence")
        # all records loaded from GenBank
        records = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))
        for ref_seq_id in args.reference:
            print(ref_seq_id)
            if os.path.exists(Path(os.path.split(args.input_file)[0],ref_seq_id+".fasta")):
                ref_f_name = str(Path(os.path.split(args.input_file)[0],ref_seq_id+".fasta"))
                ref_seq = list(SeqIO.parse(ref_f_name, "fasta"))[0]
                ref_records_list.append(ref_seq)
            else:
                if ref_seq_id in records.keys():
                     ref_seq = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))[ref_seq_id]
                     ref_records_list.append(ref_seq)
                     ref_f_name = os.path.join(out_directory,ref_seq_id+'.fasta')
                     SeqIO.write(ref_seq, ref_f_name, "fasta")
                else: 
                    ref_f_name = fetch_from_GB(ref_seq_id)
            ref_files_list.append(ref_f_name)

        # writes reference sequence to one file
        reference_file_n = str(Path(os.path.split(args.input_file)[0],"reference.fasta"))
        SeqIO.write(ref_records_list, reference_file_n, "fasta")

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

        input_dir = os.path.splitext(args.input_file)[0]
        for ref_seq_id, ref_seq_file in zip(args.reference,ref_files_list):
            # founds coverage of reference sequence by sequences downloaded from GenBank
            pos_coverage = find_coverage(args.input_file, ref_seq_file, args.path_blast)
            # writes coverage to the text file
            fout_cov = open(input_dir+"_"+ref_seq_id+"_cov.txt", "w")
            fout_cov.write(",".join(str(x) for x in pos_coverage))
            fout_cov.close()


        out_dir = str(Path(input_dir, "blast_out"))
        print('Merging blast results from different reference sequences')
        merges_coverage(out_dir, input_dir, reference_aln_file_n, args.title)
