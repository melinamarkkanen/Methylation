"""Python script to calculate information content from methylation site context sequences 
and plotting the sequence logos using logomaker library."""
from collections import defaultdict
import os
import sys
from datetime import datetime
import logging
import time
import logomaker
import pandas as pd
import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
from scoring_matrices import parse_gff
import argparse

def set_up_logger():
    """Function to set up the logger for create_logos script.
    return: the logger instance"""
    logger = logging.getLogger('create_logos')
    logger.setLevel(logging.INFO)

    script_directory = os.path.dirname(os.path.abspath(__file__))
    log_directory = os.path.join(script_directory, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    log_file = os.path.join(log_directory, f'create_logos_{datetime.now().strftime("%Y%m%d_%H%M")}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    return logger

def calculate_information_content(sequences: list) -> tuple:
    """Function to calculate information content of given sequences.
    Converting to heights in bits.
    Parameters:
    - sequences (list): list of sequences to calculate information content.
    Returns:
    - frequency_matrix (pd.DataFrame): matrix with the frequency of each base in each position.
    - R_i_list (list): list of information content values for each position.
    - R_i_df (pd.DataFrame): DataFrame with the information content values.
    - heights (dict): dictionary with the height of each nucleotide in each position.

    """
    n = len(sequences)
    s = 4
    sequence_df = pd.DataFrame([list(seq) for seq in sequences])
    sequence_df = sequence_df.apply(lambda col: col.map(lambda x: x if x in ['A', 'C', 'G', 'T'] else np.nan))
    base_counts = sequence_df.apply(pd.Series.value_counts).fillna(0).astype(int)

    frequency_matrix = base_counts.divide(n)
    R_i_list = []
    e_n = (1 / np.log(2)) * ((s - 1) / (2 * n))
    heights = {base: [] for base in ['A', 'C', 'G', 'T']}

    for i in range(frequency_matrix.shape[1]):
        freq = frequency_matrix.iloc[:, i]
        if freq.max() == n:
            H_i = 0
        else:
            H_i = entropy(freq[freq > 0], base=2)

        R_i = max(0, np.log2(s) - (H_i + e_n))


        for base in ['A', 'C', 'G', 'T']:
            heights[base].append(freq.get(base, 0) * R_i)

        R_i_list.append(R_i)

    R_i_df = pd.DataFrame(R_i_list, columns=['R_i'])

    return frequency_matrix, R_i_list, R_i_df.T, heights

def plot_seq_logos(heights: dict, contig_id: str, mod_type: str, output_dir: str, logger: logging.Logger):
    """Plot sequence logos using logomaker and save the figure as a png file.
    
    Parameters:
    - heights (dict): dictionary with the height of each nucleotide in each position.
    - contig_id (str): ID of the contig.
    - mod_type (str): type of modification.
    - output_dir (str): path to the output directory.
    - logger (logging.Logger): logger object to log info. 
    
    Returns:
    - None
    """

    figure_path = os.path.join(output_dir, f"{contig_id}_{mod_type}_logo.png")
    if os.path.exists(figure_path):
        logger.info("File %s already exists", figure_path)
        return
    nucleotide_df = pd.DataFrame(heights)
    fig, ax = plt.subplots(figsize=(10, 3))
    logomaker.Logo(nucleotide_df, ax=ax, flip_below=True)
    ax.set_title(f'Sequence Logo for: {contig_id} modification type: {mod_type}')
    ax.set_ylim(0, 2)
    sequence_length = nucleotide_df.shape[0]
    half_length = sequence_length // 2

    x_ticks = list(range(-half_length, half_length + 1, 5))
    ax.set_xticks([i + half_length for i in x_ticks])
    ax.set_xticklabels(x_ticks)

    plt.savefig(figure_path)
    plt.close(fig)

def main():
    """Main function to parse input arguments and directly calculate and plot sequence logos."""
    start_time = time.time()
    logger = set_up_logger()

    parser = argparse.ArgumentParser(description="Create sequence logos from input directory of GFF and FASTA files.")
    parser.add_argument("input_directory", help="Path to the directory containing input GFF and FASTA files.")
    parser.add_argument("output_directory", help="Path to the directory where output files will be saved.")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_directory)
    output_dir = os.path.abspath(args.output_directory)

    if not os.path.isdir(input_dir):
        logger.error("The input path is not a directory: %s", input_dir)
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    gff_files = [
        os.path.join(input_dir, file) for file in os.listdir(input_dir)
        if file.endswith('.gff')]

    if not gff_files:
        logger.warning("No GFF files found in the input directory")
        sys.exit(1)

    for file in gff_files:
        logger.info("Processing file: %s", file)
        context_sequences, mod_types = parse_gff(file)
        base_name = os.path.basename(file)
        contig_id = base_name.split('_basemods.gff')[0]
        for mod_type in mod_types:
            context_seqs = context_sequences.get(mod_type, [])
            if len(context_seqs) < 20:
                logger.warning("Skipping %s for %s due to insufficient data", contig_id, mod_type)
                continue
            _, _, _, heights = calculate_information_content(context_seqs)
            
            plot_seq_logos(heights, contig_id, mod_type, output_dir, logger)
    elapsed_time = time.time() - start_time
    logger.info("Script executed in %.2f seconds", elapsed_time)
if __name__ == "__main__":
    main()