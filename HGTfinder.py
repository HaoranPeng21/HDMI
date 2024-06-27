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



def build_genome_contig_index(genome_dir):
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    genome_to_contigs = {}
    for genome in genome_files:
        genome_name = os.path.splitext(genome)[0]
        genome_to_contigs[genome_name] = {}
        for record in SeqIO.parse(os.path.join(genome_dir, genome), 'fasta'):
            genome_to_contigs[genome_name][record.id] = len(record.seq)
    return genome_to_contigs


def run_subprocess_command_with_pipes(cmd):
    process = subprocess.run(cmd, shell=True, check=True)
    return process

def combine_fasta_files(genome_path, combined_output_path):
    with open(combined_output_path, 'w') as combined:
        for root, _, files in os.walk(genome_path):
            for file in files:
                if (file.endswith('.fa') or file.endswith('.fasta')) and file != os.path.basename(combined_output_path):
                    with open(os.path.join(root, file), 'r') as f:
                        combined.write(f.read())


def generate_coverage_data(sample_id, read1, read2, output, combined_genome_index, threads):
    output_bam = os.path.join(output, f'{sample_id}_combined_genome.bam')

    # Execute bowtie2
    cmd1 = (f"bowtie2 -a --very-sensitive -x {combined_genome_index} -1 {read1} -2 {read2} -p {threads} --no-unal | "
            f"samtools view -q30 -@ {threads} -bS - | "
            f"samtools sort -@ {threads} -o {output_bam}")
    run_subprocess_command_with_pipes(cmd1)
    return output_bam



def index_exists(index_prefix):
    return any(os.path.exists(index_prefix + ext) for ext in ['.1.bt2', '.1.bt2l'])



def get_pair_reads(bam,contig,position):

    reads = list(bam.fetch(contig, position, position + 1))
    forward_reads = [read for read in reads if not read.is_reverse and read.pos < position-10 and read.pos + read.alen > position+10]
    reverse_reads = [read for read in reads if read.is_reverse and read.pos < position-10 and read.pos + read.alen > position+10]

    return forward_reads, reverse_reads

def get_spanning_reads(bam, q_contig_id, q_start, q_end,s_contig_id,s_start,s_end,sth):
    """Return a tuple containing lists of forward and reverse reads crossing the given position in the BAM file."""

    q_start_forward_reads, q_start_reverse_reads = get_pair_reads(bam, q_contig_id, q_start)
    q_end_forward_reads, q_end_reverse_reads = get_pair_reads(bam, q_contig_id, q_end)
    s_start_forward_reads, s_start_reverse_reads = get_pair_reads(bam, s_contig_id, s_start)
    s_end_forward_reads, s_end_reverse_reads = get_pair_reads(bam, s_contig_id, s_end)
    q_HGT_site_loose = len(q_start_forward_reads) + len(q_start_reverse_reads) >= sth and len(q_end_forward_reads) + len(q_end_reverse_reads) >= sth
    q_HGT_site_strict = len(q_start_forward_reads) > 0 and len(q_start_reverse_reads) > 0 and len(q_end_forward_reads) > 0 and len(q_end_reverse_reads) > 0
    s_HGT_site_loose = len(s_start_forward_reads) + len(s_start_reverse_reads) >= sth and len(s_end_forward_reads) + len(s_end_reverse_reads) >= sth
    s_HGT_site_strict = len(s_start_forward_reads) > 0 and len(s_start_reverse_reads) > 0 and len(s_end_forward_reads) > 0 and len(s_end_reverse_reads) > 0 

    return q_HGT_site_loose, q_HGT_site_strict,s_HGT_site_loose,s_HGT_site_strict

def calculate_output_value(bam,details,sth):
    q_contig_id, q_start, q_end, hgt_length, genome_name_mag1,s_contig_id,s_start, s_end, genome_name_mag2 = details
    # For q_contig_id
    q_HGT_site_loose, q_HGT_site_strict,s_HGT_site_loose,s_HGT_site_strict = get_spanning_reads(bam, q_contig_id, q_start, q_end,s_contig_id,s_start,s_end,sth)

    mag1_value = 2 if q_HGT_site_loose else 0
    mag2_value = 2 if s_HGT_site_loose else 0

    loose_mode = f"{mag1_value:01}{mag2_value:01}"

    mag1_value1 = 2 if q_HGT_site_strict and q_HGT_site_loose else 0
    mag2_value1 = 2 if s_HGT_site_strict and s_HGT_site_loose else 0

    strict_mode = f"{mag1_value1:01}{mag2_value1:01}"

    return loose_mode, strict_mode

def load_coverage_data(bam):
    coverage_data = {}

    contig_lengths = {ref['SN']: ref['LN'] for ref in bam.header['SQ']}

    for pileup_column in bam.pileup(min_mapping_quality=30, stepper="nofilter"):
        contig = bam.get_reference_name(pileup_column.reference_id)
        position = pileup_column.pos
        coverage = pileup_column.n

        if contig not in coverage_data:
            coverage_data[contig] = np.zeros(contig_lengths[contig], dtype=np.int16)
        coverage_data[contig][position] = np.int16(coverage)

    return coverage_data


def compute_HGT_coverage(coverage_data, contig_id, start, end):
    """Computes the percentage of positions with coverage greater than zero for a given range."""
    if contig_id not in coverage_data:
        return 0

    if start > end:
        start, end = end, start

    covered_positions = np.count_nonzero(coverage_data[contig_id][start:end + 1])
    total_positions = end - start + 1
    coverage_percentage = (covered_positions / total_positions)

    return coverage_percentage


def extract_HGT_coverage(table_file):
    df = pd.read_csv(table_file, sep=",")
    df["HGT_ID"] = ["HGT_seq_" + str(i + 1) for i in range(len(df))]
    df.set_index("HGT_ID", inplace=True)
    df['details'] = df[['query_congtig_id', 'q_start', 'q_end', 'Length', 'MAG 1', "subject_contig_id", "s_start", "s_end", "MAG 2"]].values.tolist()
    HGT = dict(zip(df.index, df['details']))
    output_path = os.path.join(os.path.dirname(table_file), 'HGT_events.csv')
    df.to_csv(output_path, sep=",")
    return HGT

def calculate_coverage_value(mag1_ratio, mag2_ratio):
    #mag1_value = 2 if mag1_ratio >= 0.90  else 0
    #mag2_value = 2 if mag2_ratio >= 0.90  else 0
    mag1_value = round(mag1_ratio*100)
    mag2_value = round(mag2_ratio*100)

    return f"{mag1_value:01}_{mag2_value:01}"

def process_HGT(coverage_data, genome_index, details, num_samples=1000):
    q_contig_id, q_start, q_end, hgt_length, genome_name_mag1,s_contig_id,s_start, s_end, genome_name_mag2 = details
    ### Compute mag1 coverage
    hgt_coverage_mag1 = compute_HGT_coverage(coverage_data, q_contig_id, q_start, q_end)

    ### compute mag2 covergage
    hgt_coverage_mag2 = compute_HGT_coverage(coverage_data, s_contig_id, s_start, s_end)

    return hgt_coverage_mag1,hgt_coverage_mag2


def main():
    parser = argparse.ArgumentParser(description='HGT detection pipeline step2')

    parser.add_argument('-r1', '--read1', required=True, help='Path to read1 FASTA.')
    parser.add_argument('-r2', '--read2', required=True, help='Path to read2 FASTA.')
    parser.add_argument('-i', '--sample_id', default="sample", help='Sample ID.')
    parser.add_argument('-o', '--output', default="./", help='Working directory.')
    parser.add_argument('-mag_dir', '--genome_path', default="./", help='Directory containing genomes.')
    parser.add_argument('-table_dir', '--HGT_table_path', default="./", help='Path to the HGT table.')
    parser.add_argument('-threads', '--threads', default=1, type=int, help='Number of threads.')
    parser.add_argument('-sth', '--sth', default=3, type=int, help='reads span sites number')

    args = parser.parse_args()

    # Building the combined genome index if necessary
    combined_genome = os.path.join(args.genome_path, 'combined_genome.fa')
    os.makedirs(f"{args.output}", exist_ok=True)
    combined_genome_index = os.path.join(args.output, 'combined_genome_index')

    if not index_exists(combined_genome_index):
        print("Generating combined genome index...")
        combine_fasta_files(args.genome_path, combined_genome)
        run_subprocess_command_with_pipes(f"bowtie2-build {combined_genome} {combined_genome_index}")
        #os.remove(combined_genome)
    else:
        print("Index files already exist. Skipping the bowtie2-build step.")
        print("Processing HGT data...")

    genome_index = build_genome_contig_index(args.genome_path)

    # Align reads to combined genome
    bam_file = os.path.join(args.output, f'{args.sample_id}_combined_genome.bam')

    if not os.path.exists(bam_file):
        bam_file = generate_coverage_data(args.sample_id, args.read1, args.read2, args.output, combined_genome_index, args.threads)
    else:
        print(f"{bam_file} already exists. Skipping alignment step.")

    bai_file = bam_file + ".bai"
    if not os.path.exists(bai_file):
        run_subprocess_command_with_pipes(f"samtools index {bam_file}")
    else:
        print(f"{bai_file} already exists. Skipping indexing step.")

    bam = pysam.AlignmentFile(bam_file, "rb")
    # Extract coverage data
    coverage_data = load_coverage_data(bam)

    # Process each HGT event
    hgt_data = extract_HGT_coverage(args.HGT_table_path)

    # Save readsplit results
    loose_output = {}
    strict_output = {}

    for hgt_id, details in hgt_data.items():
        loose_value,strict_value = calculate_output_value(bam, details,args.sth)

        loose_output[hgt_id] = loose_value
        strict_output[hgt_id] = strict_value

    bam.close()
    os.makedirs(f"{args.output}/read_split", exist_ok=True)
    loose_file = os.path.join(args.output, 'read_split' ,f'{args.sample_id}_loose.csv')
    loose_df = pd.DataFrame(list(loose_output.items()), columns=['HGT_ID', f'LLD_{args.sample_id}'])
    loose_df.to_csv(loose_file, sep=",", index=False)

    strict_file = os.path.join(args.output,'read_split', f'{args.sample_id}_strict.csv')
    strict_df = pd.DataFrame(list(strict_output.items()), columns=['HGT_ID', f'LLD_{args.sample_id}'])
    strict_df.to_csv(strict_file, sep=",", index=False)

    # Save coverage results
    output_values = {}
    for hgt_id, details in hgt_data.items():
        ratio_mag1, ratio_mag2 = process_HGT(coverage_data, genome_index, details)
        #print(f"HGT ID: {hgt_id}, Value1: {ratio_mag1} Value2: {ratio_mag2}")
        value = calculate_coverage_value(ratio_mag1, ratio_mag2)
        output_values[hgt_id] = value
    os.makedirs(f"{args.output}/fraction", exist_ok=True)
    output_file = os.path.join(args.output,'fraction', f'{args.sample_id}_fraction.csv')
    output_df = pd.DataFrame(list(output_values.items()), columns=['HGT_ID', f'LLD_{args.sample_id}'])
    output_df.to_csv(output_file, sep=",", index=False)

    if os.path.exists(bam_file):
        os.remove(bam_file)

    if os.path.exists(bai_file):
        os.remove(bai_file)

if __name__ == '__main__':
    main()
