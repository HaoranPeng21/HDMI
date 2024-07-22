#!/usr/bin/env python
import importlib.util
import subprocess

def check_and_install_module(module_name):
    if importlib.util.find_spec(module_name) is None:
        print(f"{module_name} is not installed. Installing...")
        subprocess.check_call(["pip", "install", module_name])
    else:
        print(f"{module_name} is already installed.")

modules = ["numpy", "pandas", "Bio", "pysam"]

for module in modules:
    check_and_install_module(module)

import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import pysam



def run_subprocess_command_with_pipes(cmd):
    process = subprocess.run(cmd, shell=True, check=True)
    return process

def generate_coverage_data(sample_id, read1, read2, output, threads):
    output_bam = os.path.join(output, f'{sample_id}_simi.bam')
    output_index = os.path.join(output, 'simi_sequences_index')
    # Execute bowtie2
    cmd1 = (f"bowtie2 -a --very-sensitive -x {output_index} -1 {read1} -2 {read2} -p {threads} --no-unal | "
            f"samtools view -@ {threads} -bS - | "
            f"samtools sort -@ {threads} -o {output_bam}")
    run_subprocess_command_with_pipes(cmd1)
    
    return output_bam



def get_pair_reads(bam_file,contig,position):
    
    reads = list(bam_file.fetch(contig, position, position + 1))
    forward_reads = [read for read in reads if not read.is_reverse and read.pos < position-10 and read.pos + read.alen > position+10]
    reverse_reads = [read for read in reads if read.is_reverse and read.pos < position-10 and read.pos + read.alen > position+10]
        
    return forward_reads, reverse_reads
    
def get_spanning_reads(bam_file, contig_name,position):
    """Return a tuple containing lists of forward and reverse reads crossing the given position in the BAM file."""
        
    forward_reads, reverse_reads = get_pair_reads(bam_file, contig_name, position)

    reads_number = len(forward_reads)+len(reverse_reads)
        
    return reads_number

def calculate_reads_counts(bam_file, elements_info_path, output_path):
    elements_info = pd.read_csv(elements_info_path)
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []

    for _, row in elements_info.iterrows():
        hgt_id = row['HGT_ID']
        contig_name = row['contig_name']
        contig_start = row['contig_start']
        contig_end = row['contig_end']
        hgt_id_sim = row['HGT_ID_sim']
        hgt_id_sim_position = row['HGT_ID_sim_position']

        start_counts = get_spanning_reads(bam, contig_name, contig_start)
        end_counts = get_spanning_reads(bam, contig_name, contig_end)
        sim_counts = get_spanning_reads(bam, hgt_id_sim, hgt_id_sim_position)

        result = f"{hgt_id}\t{start_counts},{end_counts},{sim_counts}"
        results.append(result)

    bam.close()

    with open(output_path, 'w') as f:
        f.writelines("\n".join(results))

def main():
    parser = argparse.ArgumentParser(description='HGT detection pipeline.')

    parser.add_argument('-r1', '--read1', required=True, help='Path to read1 FASTA.')
    parser.add_argument('-r2', '--read2', required=True, help='Path to read2 FASTA.')
    parser.add_argument('-i', '--sample_id', default="sample", help='Sample ID.')
    parser.add_argument('-o', '--output', default="./", help='Working directory.')
    parser.add_argument('-table_dir', '--HGT_table_path', default="./", help='Path to the HGT table.')
    parser.add_argument('-threads', '--threads', default=1, type=int, help='Number of threads.')

    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    bam_file = os.path.join(args.output, f'{args.sample_id}_simi.bam')
    
    if not os.path.exists(bam_file):
        bam_file = generate_coverage_data(args.sample_id, args.read1, args.read2, args.output, args.threads)
    else:
        print(f"{bam_file} already exists. Skipping alignment step.")
    
    bai_file = bam_file + ".bai"
    if not os.path.exists(bai_file):
        run_subprocess_command_with_pipes(f"samtools index {bam_file}")
    else:
        print(f"{bai_file} already exists. Skipping indexing step.")
    output_path = os.path.join(args.output,"temp", f"{args.sample_id}_reads_counts_very_all.tsv")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    calculate_reads_counts(bam_file, args.HGT_table_path, output_path)

    #if os.path.exists(bam_file):
    #    os.remove(bam_file)

    #if os.path.exists(bai_file):
    #    os.remove(bai_file)

if __name__ == '__main__':
    main()

