import argparse
from .TeloVision import *

def main():
    argParser = argparse.ArgumentParser(description="TeloVision - Determine \
                                and visualise the presence of telomeres. \
                                For more information see: \
                                https://github.com/TimVerschuren/TeloVision")
    argParser.add_argument("-i", 
                           "--input", 
                           type=str, 
                           help="Location of fasta file", 
                           required=True)
    argParser.add_argument("-o", 
                           "--output", 
                           type=str, 
                           help="Name of output file", 
                           required=True)
    argParser.add_argument("-k", 
                           "--kmer_size", 
                           type=int, 
                           help="Size of k-mer sliding window. \
                            (default: 5000)")
    argParser.add_argument("-r", 
                           "--min_repeat_length", 
                           type=int, 
                           help="Minimum length for a repetitive \
                           region to be counted as a telomere. (default: 30)")
    argParser.add_argument("-s", 
                           "--sequence_size", 
                           type=int, 
                           help="Size of sequence taken from the 5' and 3' \
                           ends of the sequence for analysis. (default: 200)")
    argParser.add_argument("-g",
                           "--gap_size",
                           type=int,
                           help="Size of gaps allowed within telomeric \
                            repeats. (default: length of detected repeat)")
    argParser.add_argument("--sorted",
                           action="store_true",
                           help="Sorts the scaffolds by size if selected.")
    
    args = argParser.parse_args()
    if not args.kmer_size:
        args.kmer_size = 5000
    if not args.min_repeat_length:
        args.min_repeat_length = 30
    if not args.sequence_size:
        args.sequence_size = 200
    if not args.gap_size:
        args.gap_size= False
    if not args.sorted:
        args.sorted = False
    
    visualiseGC(args.input, 
            findTelomeres(args.input, args.output, args.sorted)\
                .telomere_position(args.min_repeat_length, 
                                   args.sequence_size,
                                   args.gap_size),
                args.output, args.sorted
                ).gc_content(args.kmer_size)
    
    print(f"Figure saved as {args.output}.html. Repeat data saved as \
{args.output}_info.tsv. Thank you for using TeloVision!\n")

if __name__ == "__main__":
    main()
