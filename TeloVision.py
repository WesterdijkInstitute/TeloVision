#!/usr/bin/env python

"""
TeloVision is a Python package which determines the presence of telomeres 
and visualises scaffolds in genome assemblies.

usage: telovision [-h] -i INPUT -o OUTPUT [-k KMER_SIZE] [-r MIN_REPEAT_LENGTH]
 [-s SEQUENCE_SIZE] [-g GAP_SIZE] [--sorted]

"""

"""Import Statements"""
import pandas as pd
from Bio import SeqIO
import plotly.express as px
import Levenshtein

"""Authorship Information"""
__author__ = "Tim Verschuren"
__credits__ = ["Tim Verschuren", "Jérôme Collemare"]

__licence__ = "MIT"
__date__ = "20-09-2023"
__version__ = "0.3.2"
__maintainer__ = "Tim Verschuren"
__email__ = "t.verschuren@wi.knaw.nl"
__status__ = "Development"


class findTelomeres:
    """Search for repeating sequences and determine presence of telomeres.

    Attributes:
        fasta_file (str): Path to fasta file.
        output (str): Name of output file.
        sort (bool): True if scaffolds should be sorted by size.
    """
    def __init__(self, fasta_file: str, output: str, sort=False):
        if sort == True:
            self.fasta = dict(sorted(read_fasta(fasta_file).items(), 
                                    key=lambda item: len(item[1]), 
                                    reverse=True))
        else:
            self.fasta = read_fasta(fasta_file)
        self.output = output

    def identify_repeats(self, seq_slice: str, gap: int) -> int:
        """Iterates over nucleotide sequence and uses
        k-mers of various sizes to determine the presence
        of a repeating sequence.

        Attributes:
            seq_slice (str): A sub portion of a nucleotide sequence.
            gap (int): Size of gap allowed within repeating sequence to count
            as one repeat.

        Returns:
            repeat_length (int): The length of the repeating
            sequence, if present.
        """
        kmer_dict = {}
        repeat_coord = {}
        rep_qc_dict = {}
        for k in range(5,25):
            if gap == False:
                gap = k
            for i in range(0, len(seq_slice) - k + 1):
                # If a match is found, initialize k-mer repeat search.
                if seq_slice[i:i+k] == seq_slice[i+k:i+2*k] or \
                    seq_slice[i:i+k] == seq_slice[i+k+1:i+2*k+1]:
                    j = i
                    while True:
                        # If a match is found, move over to the next k-mer.
                        if self.nuc_diff(seq_slice[i:i+k], 
                                    seq_slice[j:j+k]) < 2 or \
                            self.nuc_diff(seq_slice[i:i+k], 
                                     seq_slice[j+1:j+k+1]) < 2:
                            j += k
                        else:
                            # If the gap between current repeat and previous
                            # repeat is small enough, connect the two
                            if seq_slice[i:i+k] in kmer_dict:
                                if abs(i - repeat_coord[seq_slice[i:i+k]])\
                                    <= gap:
                                    kmer_dict[seq_slice[i:i+k]] += \
                                        seq_slice[
                                            repeat_coord[seq_slice[i:i+k]]:j
                                            ]
                                    repeat_coord[seq_slice[i:i+k]] = j
                                # If repeats are too distant, save the largest
                                # one.
                                elif len(kmer_dict[seq_slice[i:i+k]]) \
                                    - len(seq_slice[i:j]) < 0:
                                    kmer_dict[seq_slice[i:i+k]] \
                                        = seq_slice[i:j]
                                    repeat_coord[seq_slice[i:i+k]] = j
                            # If repeat is new, add to repeat dict and save
                            # its final index.
                            else:
                                kmer_dict[seq_slice[i:i+k]] = \
                                    seq_slice[i:j]
                                repeat_coord[seq_slice[i:i+k]] = j
                            break

        # Retrieve length of found repeats from dictionary
        seq_len = [len(seq) for seq in list(kmer_dict.values())]
        if len(seq_len) == 0:
            repeat_sequence = "NA"
            repeat = "NA"
        # Retrieve full sequence containing repeat and the repeating sequence.
        else:
            # Extract 2 largest repeats
            sorted_kmer_dict = dict(list(sorted(kmer_dict.items(), 
                                            key=lambda item: len(item[1]), 
                                            reverse=True))[:3])

            # If all keys have the same length, select the longest repeat.
            top_keys = list(sorted_kmer_dict.keys())
            if all(len(key) == len(top_keys[0]) for key in top_keys):
                for key, value in sorted_kmer_dict.items():
                    rep_qc_dict[key] = len(value)

                # Select position of largest repeat in dictionary.
                qc_score = [score for score in list(rep_qc_dict.values())]
                position = qc_score.index(max(qc_score))
                repeat_sequence = list(sorted_kmer_dict.values())[position]
                repeat = list(sorted_kmer_dict.keys())[position]                
            
            # If not all keys have the same length, determine repeat quality
            else:
                diff_keys = {}
                # If no key selected yet
                for key, value in sorted_kmer_dict.items():
                    rep_key = self.repeat_pattern(key)
                    if len(diff_keys) == 0:
                        diff_keys[rep_key] = value
                    else:
                        # Determine if key of similar length is present
                        diff_key_len = list(len(key_) \
                                            for key_ in diff_keys.keys())
                        if len(rep_key) not in diff_key_len:
                            diff_keys[rep_key] = value
                            
                for key, value in diff_keys.items():
                    rep_qc_dict[key] = self.repeat_qc(key, value)
            
                # Select position of largest repeat in dictionary.
                qc_score = [score for score in list(rep_qc_dict.values())]
                position = qc_score.index(max(qc_score))
                repeat_sequence = list(diff_keys.values())[position]
                repeat = list(diff_keys.keys())[position]

        return repeat_sequence, repeat
    
    def telomere_position(self, rep_len=30, 
                          seq_size=200, gap=False) -> pd.DataFrame:
        """Loop over the scaffolds of a fasta file and determine
        whether a telomere is present at the beginning and end of
        each scaffold.

        Attributes:
            rep_len (int): Minimum length of a repetitive sequence
            to qualify as a telomeric repeat.
            seq_size (int): Size of sequence taken from the top
            and bottom of the scaffolds. 
            gap (int): Size of gap allowed within repeating sequence to count
            as one repeat. 

        Returns:
            df (pd.Dataframe): Dataframe containing the name, length
            and presence of telomeres for each scaffold.

        Yields:
            A tsv file containing information about the repeats.
        """
        telomeres = {}
        telo_pos = []
        scaffolds = []
        lengths = []
        scaf_data = []
        len_data = []
        gc_data = []
        repeat = []
        repeat_sequence = []
        telo_class = []

        for key, value in self.fasta.items():
            # Initialize identify repeats for 5' and 3' regions.
            telo_bin = []
            scaffolds.append(key)
            lengths.append(len(value))
            five_prime = self.identify_repeats(value[:seq_size], gap)
            three_prime = self.identify_repeats(value[-seq_size:], gap)

            # Gather repeat length and GC-content data for repeat info file.
            scaf_data.append(f"{key}_5'")
            scaf_data.append(f"{key}_3'")
            repeat.append(five_prime[1])
            repeat.append(three_prime[1])

            if five_prime[0] == "NA":
                len_data.append("NA")
            else:
                len_data.append(len(five_prime[0]))
            if three_prime[0] == "NA":
                len_data.append("NA")
            else:
                len_data.append(len(three_prime[0]))

            if five_prime[0] == "NA":
                gc_data.append("NA")
            else:
                gc_data.append(self.calculate_gc(five_prime[0]))
            if three_prime[0] == "NA":
                gc_data.append("NA")
            else:
                gc_data.append(self.calculate_gc(three_prime[0]))

            # Determine if a repeat is telomeric based on the GC-content of
            # the repeat.
            if len(five_prime[0]) >= rep_len and \
                self.calculate_gc(five_prime[0]) > 0.2:
                telo_bin.append(1)
                telomeres[key] = telo_bin
                telo_class.append("Y")
            if len(five_prime[0]) < rep_len or \
                self.calculate_gc(five_prime[0]) <= 0.2:
                telo_bin.append(0)
                telomeres[key] = telo_bin
                telo_class.append("N")
            if len(three_prime[0]) >= rep_len and \
                self.calculate_gc(three_prime[0]) > 0.2:
                telo_bin.append(1)
                telomeres[key] = telo_bin
                telo_class.append("Y")
            if len(three_prime[0]) < rep_len or \
                self.calculate_gc(three_prime[0]) <= 0.2:
                telo_bin.append(0)
                telomeres[key] = telo_bin
                telo_class.append("N")

            # Save repeating sequence.
            repeat_sequence.append(five_prime[0])
            repeat_sequence.append(three_prime[0])

        # Save telomere position data
        for value in telomeres.values():
            telo_pos.append(value)

        # Create a DataFrame containing information about the repeats and
        # output as a .tsv file.
        telo_data = pd.DataFrame(data={"Repeat": repeat, 
                                "Length": len_data,
                                "GC%": gc_data, 
                                "Telomere": telo_class,
                                "Repetitive Sequence": repeat_sequence
                                }, 
                                index=scaf_data)
        telo_data.to_csv(f"{self.output}_info.tsv", sep="\t")

        # Save scaffold data as DataFrame for visualization.
        df = pd.DataFrame(data={'Scaffolds': scaffolds, 
                                'Lengths': lengths, 
                                'Telomeres': telo_pos})
        return df

    def calculate_gc(self, sequence: str) -> int:
        """Calculates the GC content of a given sequence.

        Attributes:
        sequence (str): Nucleotide sequence.

        Returns:
        Integer between 0 and 1.
        """
        return (sequence.count("G") + sequence.count("C"))/len(sequence)

    def nuc_diff(self, seq1: str, seq2: str) -> int:
        """Calculate number of different nucleotides between two sequences.

        Attributes:
        seq1 (str): Nucleotide sequence
        seq2 (str): Nucleotide sequence

        Returns:
        nuc_diff (int): Number of different nucleotides. 
        """
        if len(seq1) != len(seq2):
            return len(seq1)
        else:
            return sum(seq1[nuc] != seq2[nuc] for nuc in range(len(seq1)))

    def repeat_qc(self, rep: str, seq: str) -> float:
        """Determine how well a repeat fits the repeating sequence.
        
        Attributes:
            rep (str): Short repeating nucleotide sequence
            seq (str): Full repeating nucleotide sequence
        
        Returns:
            int: Quality score of repeat, higher is better.
        """
        distance = Levenshtein.distance(str(rep)*round(len(seq)/len(rep)), 
                                        str(seq)) + 1

        return len(seq)/(distance*len(rep))

    def repeat_pattern(self, rep: str) -> str:
        """Determine if the detected repeat is miss reported as a larger repeat.

        Attributes:
            rep(str): Detected repeating nucleotide sequence.

        Returns:
            str: The pattern within the detected repeat.
        """

        pattern = rep
        # Check for different sized repeats.
        for i in range(5, len(rep) // 2 + 1):
            # If string is found to repeat it self, save the repeat.
            if(not len(rep) % len(rep[:i]) and \
            rep[:i] * (len(rep)//len(rep[:i])) == rep):
                pattern = rep[:i]

        return pattern


class visualiseGC:
    """Calculation and visualisation of GC content and telomere position.

    Attributes:
        fasta_file (str): Path to fasta file.
        telo_df (pd.Dataframe): Dataframe containing the name, length
        and presence of telomeres for each scaffold.
        output (str): Name of output file.
        sort (bool): True if scaffolds should be sorted by size.
    """
    def __init__(self, fasta_file: str, telo_df: pd.DataFrame, 
                 output: str, sort=False):
        if sort == True:
            self.fasta = dict(sorted(read_fasta(fasta_file).items(), 
                                    key=lambda item: len(item[1]), 
                                    reverse=True))
        else:
            self.fasta = read_fasta(fasta_file)
        self.scaffolds = []
        self.GC_cont = []
        self.length_list = []
        self.telo_df = telo_df
        self.output = output

    def gc_content(self, k=5000) -> px:
        """Calculates the GC content of each scaffold using 
        a k-mer sliding window. GC content and present telomeres
        are subsequently visualised using a plotly bar plot. 
        Telomeres will be visualised as black bars at either 
        the top or bottom of a scaffold.

        Attributes:
            k (int): Size of k-mer sliding window.

        Yields:
            px.bar: html file of plotly bar plot.
        """
        for key, value in self.fasta.items():
            # Check if a telomeres are present on the scaffolds 5' end.
            if list(self.telo_df.loc[self.telo_df["Scaffolds"] \
                                     == key]["Telomeres"])[0][0] == 1:
                self.GC_cont.append(0)
                self.scaffolds.append(key)
                # Create a black bar for the telomere
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)
            else:
                self.GC_cont.append(100)
                self.scaffolds.append(key)
                # Create a white bar if no telomere is present.
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)

            # Iterate over scaffold k-mer wise and calculate the GC-content for
            # each k-mer. 
            for i in range(0, len(value)-k+1, k):
                self.GC_cont.append((value[i:i+k].count('C') + \
                                     value[i:i+k].count('G'))/k*100)
                self.scaffolds.append(key)
                self.length_list.append(k/1e6)
            # Check if a telomeres are present on the scaffolds 3' end.
            if list(self.telo_df.loc[self.telo_df["Scaffolds"] == \
                                     key]["Telomeres"])[0][1] == 1:
                self.GC_cont.append(0)
                self.scaffolds.append(key)
                # Create a black bar for the telomere
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)
            else:
                self.GC_cont.append(100)
                self.scaffolds.append(key)
                # Create a white bar if no telomere is present.
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)

        fig_df = pd.DataFrame(data={"Scaffolds": self.scaffolds, 
                                    "Length": self.length_list, 
                                    "GC%": self.GC_cont})
        # Create figure
        fig = px.bar(fig_df, 
                     x="Scaffolds", 
                     y = "Length", 
                     color="GC%", 
                     color_continuous_scale=px.colors.sequential.Hot, 
                     range_color=(0,fig_df["GC%"].mean()))
        
        fig.update_traces(marker_line_width = 0,
                  selector=dict(type="bar"))
        fig.update_xaxes(showgrid=False, title="Scaffolds", tickangle=45)
        fig.update_yaxes(showgrid=False, title="Physical position (Megabases)")
        fig.write_html(f"{self.output}.html")


def read_fasta(fasta_file) -> dict:
    """Read content of fasta file and store data
    in a dictonary with the scaffold names as the 
    key and the nucleotide sequence as the values.

    Attributes:
        fasta_file (str): Path to the fasta file.

    Returns:
        fasta_dict (dict): Dictionary containing 
        scaffold names as keys and nucleotide sequences
        as the values.
    """
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = (record.seq).upper()

    return fasta_dict
