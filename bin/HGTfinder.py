#!/usr/bin/env python3
"""
HGT Validation Script

Validate HGT events using read coverage analysis.

Author: Haoran Peng (penghr21@gmail.com)
GitHub: https://github.com/HaoranPeng21/HDMI
"""

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
import random



def read_group_info_and_select_representatives(group_info_file, seed=42):
    """
    Read group info file and select one representative genome per group using a fixed seed.
    This ensures reproducible results across runs.
    
    Args:
        group_info_file: Path to the group info file
        seed: Random seed for reproducible selection (default: 42)
    
    Returns:
        Dictionary mapping group numbers to selected genome names.
    """
    # Set random seed for reproducibility
    random.seed(seed)
    
    group_genomes = {}
    with open(group_info_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                genome_file, group_id = line.split()
                genome_name = os.path.splitext(genome_file)[0]  # Remove .fa/.fasta extension
                group_id = int(group_id)
                
                if group_id not in group_genomes:
                    group_genomes[group_id] = []
                group_genomes[group_id].append(genome_name)
    
    # Sort genomes within each group for consistent ordering before selection
    for group_id in group_genomes:
        group_genomes[group_id].sort()
    
    # Select one representative from each group (reproducible due to fixed seed)
    selected_representatives = {}
    for group_id, genomes in group_genomes.items():
        selected_representatives[group_id] = random.choice(genomes)
    
    print(f"Selected representative genomes (seed={seed}):")
    for group_id, genome in sorted(selected_representatives.items()):
        print(f"  Group {group_id}: {genome}")
    
    return selected_representatives

def build_genome_contig_index(genome_dir, selected_genomes=None):
    """
    Build genome contig index. If selected_genomes is provided, only index those genomes.
    """
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    genome_to_contigs = {}
    
    for genome in genome_files:
        genome_name = os.path.splitext(genome)[0]
        
        # If selected_genomes is provided, only process selected genomes
        if selected_genomes is not None:
            if genome_name not in selected_genomes.values():
                continue
        
        genome_to_contigs[genome_name] = {}
        for record in SeqIO.parse(os.path.join(genome_dir, genome), 'fasta'):
            # Add MAG prefix to contig ID to match BAM file format
            contig_id_with_prefix = f"{genome_name}_{record.id}"
            genome_to_contigs[genome_name][contig_id_with_prefix] = len(record.seq)
    
    return genome_to_contigs


def run_subprocess_command_with_pipes(cmd):
    process = subprocess.run(cmd, shell=True, check=True)
    return process

def combine_fasta_files(genome_path, combined_output_path, selected_genomes=None):
    """
    Combine FASTA files. If selected_genomes is provided, only combine those genomes.
    """
    with open(combined_output_path, 'w') as combined:
        for root, _, files in os.walk(genome_path):
            for file in files:
                if (file.endswith('.fa') or file.endswith('.fasta')) and file != os.path.basename(combined_output_path):
                    genome_name = os.path.splitext(file)[0]
                    
                    # If selected_genomes is provided, only process selected genomes
                    if selected_genomes is not None:
                        if genome_name not in selected_genomes.values():
                            continue
                    
                    with open(os.path.join(root, file), 'r') as f:
                        combined.write(f.read())


def generate_coverage_data(sample_id, read1, read2, output, combined_genome_index, threads, suffix=""):
    """Generate coverage data by aligning reads to genome index."""
    output_bam = os.path.join(output, 'intermediate', f'{sample_id}_combined_genome{suffix}.bam')

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
    
    # Add MAG prefixes to contig IDs to match BAM file format
    df['query_congtig_id_with_prefix'] = df['MAG 1'].str.replace('.fa', '', regex=False) + '_' + df['query_congtig_id']
    df['subject_contig_id_with_prefix'] = df['MAG 2'].str.replace('.fa', '', regex=False) + '_' + df['subject_contig_id']
    
    df['details'] = df[['query_congtig_id_with_prefix', 'q_start', 'q_end', 'Length', 'MAG 1', "subject_contig_id_with_prefix", "s_start", "s_end", "MAG 2"]].values.tolist()
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


# Functions from HGTfinder_abun.py for abundance calculation
def get_genome_total_length(genome_to_contigs, genome_name):
    return sum(genome_to_contigs[genome_name].values())

def get_read_count_from_fastq(fastq_file):
    if fastq_file.endswith('.gz'):
        cmd = f"zcat {fastq_file} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            line_count = int(result.stdout.strip())
            return line_count // 4  # Each read has 4 lines in FASTQ format
        else:
            raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    else:
        cmd = ['awk', '{s++}END{print s/4}', fastq_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return int(result.stdout.strip())
        else:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")

def get_average_read_length(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_length = 0
        count = 0
        for read in bam.fetch():
            if read.mapping_quality >= 30:
                total_length += len(read.query_sequence)
                count += 1
    return total_length / count

def load_coverage_data_for_abundance(bam_file, fastq1, fastq2):
    coverage_data = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        contig_lengths = {ref['SN']: ref['LN'] for ref in bam.header['SQ']}
        coverage_data = {contig: np.zeros(length, dtype=np.int16) for contig, length in contig_lengths.items()}
        for pileup_column in bam.pileup(min_base_quality=30, min_mapping_quality=40, stepper="nofilter"):
            contig = bam.get_reference_name(pileup_column.reference_id)
            position = pileup_column.pos
            coverage = sum(1 for read in pileup_column.pileups if not read.is_del and not read.is_refskip)
            
            if contig not in coverage_data:
                coverage_data[contig] = np.zeros(contig_lengths[contig], dtype=np.int16)
            coverage_data[contig][position] = np.int16(coverage)

    reads1 = get_read_count_from_fastq(fastq1)
    reads2 = get_read_count_from_fastq(fastq2)
    total_reads = reads1 + reads2
    
    return coverage_data, total_reads

def get_avg_and_median_coverage_for_genome(coverage_data, genome_to_contigs, genome_name, genome_length):
    all_coverages = []

    for contig, lengths in genome_to_contigs[genome_name].items():
        if contig in coverage_data:
            all_coverages.extend(coverage_data[contig])

    if not all_coverages:
        print(f"Warning: No coverage data found for genome {genome_name}")
        return 0.0, 0.0

    avg_coverage = np.sum(all_coverages) / genome_length
    median_coverage = np.median(all_coverages)

    return avg_coverage, median_coverage

def calculate_relative_abundance(avg_coverage, genome_length, avg_read_length, total_reads):
    K = avg_coverage
    L = avg_read_length
    S = genome_length
    T = total_reads

    return (K * S / L) / T


def main():
    parser = argparse.ArgumentParser(description='HGT detection pipeline step2')

    parser.add_argument('-r1', '--read1', required=True, help='Path to read1 FASTA.')
    parser.add_argument('-r2', '--read2', required=True, help='Path to read2 FASTA.')
    parser.add_argument('-i', '--sample_id', default="sample", help='Sample ID.')
    parser.add_argument('-o', '--output', default="./", help='Working directory.')
    parser.add_argument('-mag_dir', '--genome_path', default="./", help='Directory containing genomes.')
    parser.add_argument('-table_dir', '--HGT_table_path', default="./", help='Path to the HGT table.')
    parser.add_argument('-group_info', '--group_info', required=True, help='Path to the group info file.')
    parser.add_argument('-seed', '--seed', default=42, type=int, help='Random seed for reproducible representative genome selection (default: 42).')
    parser.add_argument('-threads', '--threads', default=1, type=int, help='Number of threads.')
    parser.add_argument('-sth', '--sth', default=3, type=int, help='reads span sites number')

    args = parser.parse_args()

    # Select representative genomes from each group
    print("Selecting representative genomes from each group...")
    selected_representatives = read_group_info_and_select_representatives(args.group_info, args.seed)
    
    # Create output directories
    os.makedirs(f"{args.output}", exist_ok=True)
    os.makedirs(f"{args.output}/intermediate", exist_ok=True)
    
    # Building the combined genome index for selected representatives only (shared across samples)
    # Use the main output directory for index files (not the sample-specific output)
    main_output_dir = args.output.split('intermediate')[0].rstrip('/')
    index_dir = os.path.join(main_output_dir, 'index')
    os.makedirs(index_dir, exist_ok=True)
    
    combined_genome = os.path.join(index_dir, 'combined_genome_representatives.fa')
    combined_genome_index = os.path.join(index_dir, 'combined_genome_representatives_index')

    if not index_exists(combined_genome_index):
        print("Generating combined genome index for representative genomes...")
        combine_fasta_files(args.genome_path, combined_genome, selected_representatives)
        run_subprocess_command_with_pipes(f"bowtie2-build {combined_genome} {combined_genome_index}")
        os.remove(combined_genome)
    else:
        print("Representative genome index files already exist. Skipping the bowtie2-build step.")
        print("Processing HGT data...")

    # Build genome index for all genomes (for HGT processing)
    genome_index = build_genome_contig_index(args.genome_path)
    
    # Build genome index for selected representatives only (for abundance calculation) 
    representatives_genome_index = build_genome_contig_index(args.genome_path, selected_representatives)

    # Align reads to representative genomes for abundance calculation
    bam_file_representatives = os.path.join(args.output, 'intermediate', f'{args.sample_id}_combined_genome_representatives.bam')

    if not os.path.exists(bam_file_representatives):
        bam_file_representatives = generate_coverage_data(args.sample_id, args.read1, args.read2, args.output, combined_genome_index, args.threads, "_representatives")
    else:
        print(f"{bam_file_representatives} already exists. Skipping alignment step.")

    # For HGT processing, we still need alignment to all genomes
    # Check if we need to create a full genome alignment (only if doing HGT analysis)
    bam_file_all = os.path.join(args.output, 'intermediate', f'{args.sample_id}_combined_genome_all.bam')
    combined_genome_all_index = os.path.join(index_dir, 'combined_genome_all_index')
    
    # Only create full genome index and alignment if it doesn't exist
    if not os.path.exists(bam_file_all):
        if not index_exists(combined_genome_all_index):
            print("Generating combined genome index for all genomes...")
            combined_genome_all = os.path.join(index_dir, 'combined_genome_all.fa')
            combine_fasta_files(args.genome_path, combined_genome_all)
            run_subprocess_command_with_pipes(f"bowtie2-build {combined_genome_all} {combined_genome_all_index}")
            os.remove(combined_genome_all)
        
        print("Aligning reads to all genomes for HGT analysis...")
        bam_file_all = generate_coverage_data(args.sample_id, args.read1, args.read2, args.output, combined_genome_all_index, args.threads, "_all")
    else:
        print(f"{bam_file_all} already exists. Skipping full genome alignment step.")

    # Index both BAM files
    bai_file_representatives = bam_file_representatives + ".bai"
    if not os.path.exists(bai_file_representatives):
        run_subprocess_command_with_pipes(f"samtools index {bam_file_representatives}")
    else:
        print(f"{bai_file_representatives} already exists. Skipping indexing step.")

    bai_file_all = bam_file_all + ".bai"
    if not os.path.exists(bai_file_all):
        run_subprocess_command_with_pipes(f"samtools index {bam_file_all}")
    else:
        print(f"{bai_file_all} already exists. Skipping indexing step.")

    # Use the full genome BAM for HGT analysis
    bam = pysam.AlignmentFile(bam_file_all, "rb")
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

    # Calculate abundance for representative genomes
    print("Calculating abundance for representative genomes...")
    avg_read_length = get_average_read_length(bam_file_representatives)
    coverage_data_abundance, total_reads = load_coverage_data_for_abundance(bam_file_representatives, args.read1, args.read2)
    
    # Calculate median coverage for each representative genome (representing its group)
    median_output = {}
    for group_id, genome_name in selected_representatives.items():
        if genome_name in representatives_genome_index:
            genome_length = get_genome_total_length(representatives_genome_index, genome_name)
            avg_coverage, median_coverage = get_avg_and_median_coverage_for_genome(
                coverage_data_abundance, representatives_genome_index, genome_name, genome_length)
            # Use group_id as key instead of genome_name to represent the entire group
            median_output[f"Group_{group_id}"] = median_coverage
    
    # Save abundance results in abundance folder
    os.makedirs(f"{args.output}/abundance", exist_ok=True)
    median_output_file = os.path.join(args.output, 'abundance', f'{args.sample_id}_median.csv')
    median_df = pd.DataFrame(list(median_output.items()), columns=['Group_id', f'LLD_{args.sample_id}'])
    median_df.to_csv(median_output_file, sep=",", index=False)
    print(f"Median coverage results saved to: {median_output_file}")

    # Keep intermediate files for potential reuse across samples
    print(f"Intermediate files saved in: {os.path.join(args.output, 'intermediate')}")
    print(f"BAM files:")
    print(f"  - Representative genomes: {bam_file_representatives}")
    print(f"  - All genomes: {bam_file_all}")
    print(f"Analysis complete!")

if __name__ == '__main__':
    main()
