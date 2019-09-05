import os
from pathlib import Path




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str,
                        help="File with sequences in fasta or genbank format for calculating coverage")
    parser.add_argument("-q", "--query", type=str,
                        help="Query for GenBank Nucleotide database")
    parser.add_argument("-ref", "--reference", type=str,
                        help="List with reference sequence", required= True)
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
