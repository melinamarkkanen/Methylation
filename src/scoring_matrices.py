"""Analyzing Pac-Bio methylation data using GFF and FASTA files.
1. Parsing the GFF files to extract sequence information, methylation counts, context data, and context sequences.
2. Counting nucleotide counts for each contig from the FASTA files.
3. Calculating the occurrence n sized n-mers from the context sequences which are from the GFF files.
4. Writing the log transformed and normalized count matrices (PWM) into a tsv file. One file for each modification type.
5. Flatten the matrices and create one file for each modification type.
"""
import os
import sys
from collections import defaultdict, Counter
import logging
import time
from itertools import product
from datetime import datetime
import pandas as pd
import numpy as np
import argparse


def set_up_logger():
    """Function to set up the logger for scoring_matrices script."""
    logger = logging.getLogger('scoring_matrices')
    logger.setLevel(logging.INFO)

    script_directory = os.path.dirname(os.path.abspath(__file__))
    log_directory = os.path.join(script_directory, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    log_file = os.path.join(log_directory, f'scoring_matrices_{datetime.now().strftime("%Y%m%d_%H%M")}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    return logger

def parse_gff(folder_path: str) -> tuple:
    """Parse a GFF file and extract the context sequences for each modification type
    
    Parameter:
    - folder_path: absolute path to the folder of GFF and FASTA files for each contig
    
    Return:
    - context_sequences: defaultdict with modification type as key and list of context sequences as values
    - mod_types: list of modification types
    """
    context_sequences = defaultdict(list)
    mod_types = []

    with open(folder_path, 'r') as file:
        for line in file:
            if line.startswith('##source-commandline'):
                if '--identify' in line:
                    mod_types = line.split('--identify ')[1].split(' ')[0].split(',')
                    mod_types.append('modified_base')
            elif not line.startswith('#'):
                parts = line.strip().split('\t')
                methylation_type = parts[2]
                context = parts[8] # The ninth column contains the context sequence
                context = next((item.split('context=')[1] for item in context.split(';') if 'context=' in item), None) # Extract the context sequence
                IPDRatio = next((item.split('IPDRatio=')[1] for item in context.split(';') if 'IPDRatio=' in item), None)
                context_sequences[methylation_type].append((context))

    return context_sequences, mod_types



def count_bases_in_fasta(fasta_file_path:str) -> Counter:
    """Count the number of each base (A, C, G, T) in a fasta file.

    Parameter:
    - fasta_file_path (str): absolute path to the folder of fasta files
    
    Return
    - base_counts (Counter): Counter object with the counts of each base in the fasta file.
    """

    base_counts = Counter()
    with open(fasta_file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                base_counts.update(line.strip())

    return base_counts

def create_PWM(context_seqs:list, fasta_counts:dict, n:int, logger, window_size=None) -> pd.DataFrame:
    """
    Create a DataFrame with log transformed and normalized counts of n-tuples at each position.
    If window_size is specified, counts are calculated within a window around the middle of each sequence.

    Formula for log transformation and normalization:
    log((q + p0) / (p0 + 1))

    p0 is the probability of a nucleotide in the FASTA file.
    q is the probability of a nucleotide in the context sequences.
    

    Parameters:
    - context_seqs: List of sequences.
    - fasta_counts: Dictionary with counts of each nucleotide in the FASTA file.
    - n: Size of the n-tuple.
    - window_size: Size of the window around the middle (modified cite) of each sequence. If None, consider all positions.
    
    Returns:
    - DataFrame with log transformed and normalized counts of n-tuples.
    """
    possible_nucleotides = [''.join(p) for p in product('ATCG', repeat=n)]
    if window_size is not None:
        positions = range(-window_size, window_size + 1)
    else:
        if not context_seqs:
            sequence_length = 41
        else:
            sequence_length = len(context_seqs[0])
        start_position = -(sequence_length // 2)
        end_position = (sequence_length // 2)
        positions = range(start_position, end_position + 1)

    PWM = pd.DataFrame(index=possible_nucleotides, columns=positions, dtype=float).fillna(0.0)

    if not context_seqs or len(context_seqs) < 20:
        return PWM

    counts = defaultdict(int)
    total_counts_by_position = defaultdict(int)

    for seq in context_seqs:
        seq_length = len(seq)
        if window_size is not None:
            middle = seq_length // 2
            start_position = max(middle - window_size, 0)
            end_position = min(middle + window_size + 1, seq_length - n + 1)
            for i in range(start_position, end_position):
                n_tuple = seq[i:i+n]
                if 'N' not in n_tuple:
                    pos = i - middle
                    counts[(n_tuple, pos)] += 1
        else:
            for i in range(seq_length - n + 1):
                n_tuple = seq[i:i+n]
                if 'N' not in n_tuple:
                    pos = i - (seq_length // 2)
                    counts[(n_tuple, pos)] += 1

    for (nucleotide, pos), count in counts.items():
        total_counts_by_position[pos] += count

    total_nucleotides_in_fasta = sum(fasta_counts.values())

    for (nucleotide, pos), count in counts.items():
        p0 = fasta_counts.get(nucleotide, 0) / total_nucleotides_in_fasta
        q = (count + p0) / (total_counts_by_position[pos] + 1)
        normalized_count = q / p0
        PWM.at[nucleotide, pos] = np.log(normalized_count)

    return PWM

def main():
    start_time = time.time()
    logger = set_up_logger()

    parser = argparse.ArgumentParser(description="Create scoring matrices from GFF and FASTA files.")
    parser.add_argument("input_dir", help="Path to the directory containing GFF and FASTA files.")
    parser.add_argument("output_dir", help="Path to the output directory.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_directory = args.output_dir

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    logger.info("Output from scoring matrices will be saved to: %s", output_directory)

    flattened_matrices = {}

    for file in os.listdir(input_dir):
        if file.endswith('.gff'):
            contig = file.replace('_basemods.gff', '')
            logger.info("Analyzing contig %s", contig)
            gff_file = os.path.join(input_dir, file)
            fasta_file = os.path.join(input_dir, f"{contig}.fasta")
            if os.path.exists(fasta_file):
                fasta_counts = count_bases_in_fasta(fasta_file)
                context_sequences, mod_types= parse_gff(gff_file)
                
                logger.info("Found modification types: %s", mod_types)
                for mod_type in mod_types:
                    context_seqs = context_sequences.get(mod_type, [])
                    PWM_df = create_PWM(context_seqs, fasta_counts, 1, logger)
                    flattened_matrix = PWM_df.values.flatten()
                    if mod_type not in flattened_matrices:
                        flattened_matrices[mod_type] = []
                    flattened_matrices[mod_type].append((contig, flattened_matrix))
                    output_file_path = os.path.join(output_directory, f"{file.replace('_basemods.gff', '')}_{mod_type}.tsv")
                    logger.info(f"Checking if the file exists at: {output_file_path}")

                    if not os.path.exists(output_file_path):
                        PWM_df.to_csv(output_file_path, sep='\t')
                    else:
                        logger.info(f"File '{output_file_path}' already exists. Skipping save.")
            else:
                logger.warning("No corresponding FASTA file found for %s . Skipping.", file)

    sub_folder = os.path.join(output_directory, 'flattened')
    if not os.path.exists(sub_folder):
        os.makedirs(sub_folder)
    for modification_type, data in flattened_matrices.items():
        output_file = os.path.join(sub_folder, f"{modification_type}.tsv")
        with open(output_file, 'w') as f:
            for contig_name, matrix in data:
                flattened_array_str = '\t'.join(map(str, matrix))
                f.write(f"{contig_name}\t{flattened_array_str}\n")
        logger.info("Flattened matrices for modification type: %s in folder %s", modification_type, output_file)
    elapsed_time = time.time() - start_time
    logger.info("Script executed in %.2f seconds", elapsed_time)

if __name__ == "__main__":
    main()
