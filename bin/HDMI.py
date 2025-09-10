#!/usr/bin/env python3
"""
HDMI: HGT Detection from MAGs in Individual

A comprehensive command-line tool for detecting and validating HGT events.
Author: Haoran Peng (penghr21@gmail.com)
GitHub: https://github.com/HaoranPeng21/HDMI

Usage:
    HDMI detect     - Detect HGT candidates from MAGs
    HDMI index      - Build optimized genome indices and extract HGT sequences (integrates connect)
    HDMI profile    - Validate HGT events and analyze read coverage (integrates validate)
    HDMI summary    - Merge results and generate final element table (integrates merge)
    HDMI --help     - Show this help message

Author: Haoran Peng (penghr21@gmail.com)
Version: 1.0
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO


def get_script_dir():
    """Get the directory containing HDMI scripts."""
    return os.path.dirname(os.path.abspath(__file__))


def merge_fasta_files_with_duplicates(file1, file2, output_file):
    """
    Merge two FASTA files, handling duplicate sequence IDs by keeping the first occurrence.
    
    Args:
        file1: Path to first FASTA file
        file2: Path to second FASTA file  
        output_file: Path to output merged FASTA file
    """
    sequences = {}
    
    # Read sequences from first file
    if os.path.exists(file1):
        for record in SeqIO.parse(file1, "fasta"):
            if record.id not in sequences:
                sequences[record.id] = record
    
    # Read sequences from second file (only add if not already present)
    if os.path.exists(file2):
        for record in SeqIO.parse(file2, "fasta"):
            if record.id not in sequences:
                sequences[record.id] = record
    
    # Write merged sequences to output file
    with open(output_file, 'w') as out_handle:
        SeqIO.write(sequences.values(), out_handle, "fasta")


def extract_sample_prefix(read1_path):
    """Extract sample prefix from read1 filename."""
    basename = os.path.basename(read1_path)
    
    # Remove common suffixes
    suffixes = [
        '_R1.fastq.gz', '_R1.fq.gz', '_R1.fastq', '_R1.fq',
        '_1.fastq.gz', '_1.fq.gz', '_1.fastq', '_1.fq',
        '.fastq.gz', '.fq.gz', '.fastq', '.fq'
    ]
    
    for suffix in suffixes:
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
    
    # If no suffix found, remove extension
    return os.path.splitext(basename)[0]


def find_file_in_output(output_dir, filename):
    """Find a file in the output directory structure."""
    possible_paths = [
        os.path.join(output_dir, filename),
        os.path.join(output_dir, 'intermediate', filename),
        os.path.join(output_dir, 'intermediate', '01_detection', filename),
        os.path.join(output_dir, 'intermediate', '03_final', filename),
        os.path.join(output_dir, 'intermediate', '03_final', 'merged', filename),
        os.path.join(output_dir, 'intermediate', '04_connect', filename)
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # If not found in output, check current directory
    if os.path.exists(filename):
        return filename
    
    return None


def auto_create_dirs(base_output, *subdirs):
    """Automatically create directory structure."""
    full_path = os.path.join(base_output, *subdirs)
    os.makedirs(full_path, exist_ok=True)
    return full_path


def cleanup_bam_files(output_dir):
    """Clean up BAM and BAI files to save disk space."""
    print("\n=== Cleaning up BAM files ===")
    
    # Find all BAM and BAI files in the output directory
    bam_files = []
    bai_files = []
    
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith('.bam'):
                bam_files.append(os.path.join(root, file))
            elif file.endswith('.bai'):
                bai_files.append(os.path.join(root, file))
    
    # Remove BAM files
    for bam_file in bam_files:
        try:
            os.remove(bam_file)
            print(f"Removed: {bam_file}")
        except OSError as e:
            print(f"Warning: Could not remove {bam_file}: {e}")
    
    # Remove BAI files
    for bai_file in bai_files:
        try:
            os.remove(bai_file)
            print(f"Removed: {bai_file}")
        except OSError as e:
            print(f"Warning: Could not remove {bai_file}: {e}")
    
    total_removed = len(bam_files) + len(bai_files)
    if total_removed > 0:
        print(f"Cleaned up {total_removed} BAM/BAI files")
    else:
        print("No BAM/BAI files found to clean up")


def run_command(cmd, description="Running command"):
    """Execute a command and handle errors."""
    print(f"\n=== {description} ===")
    print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    try:
        result = subprocess.run(cmd, check=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with return code {e.returncode}")
        raise


def cmd_detect(args):
    """HDMI detect: HGT candidate detection."""
    
    # Handle count-only mode
    if hasattr(args, 'count_only') and args.count_only:
        print("=== HDMI Detect: Genome Pair Count and Performance Estimation ===")
        
        # Count genomes
        genome_files = []
        for root, _, files in os.walk(args.genome_path):
            for file in files:
                if file.endswith(('.fa', '.fasta')) and not file.startswith('combined_genome'):
                    genome_files.append(os.path.join(root, file))
        
        genome_count = len(genome_files)
        total_pairs = genome_count * (genome_count - 1) // 2
        
        print(f"Genome count: {genome_count}")
        print(f"Total genome pairs: {total_pairs}")
        
        # Performance estimation (350k pairs ~ 8h/10CPU)
        estimated_hours = (total_pairs / 350000) * 8
        estimated_cpu_hours = estimated_hours * 10
        
        print(f"Estimated processing time: {estimated_hours:.1f} hours")
        print(f"Estimated CPU hours: {estimated_cpu_hours:.1f} CPU-hours")
        print(f"Recommended: Use {max(1, min(10, int(estimated_hours)))} CPU cores")
        
        return
    
    script_path = os.path.join(get_script_dir(), 'HGTdetect.py')
    
    # Set default batch parameters (always use batch mode)
    # Parameters now have defaults in argparse, so no need to set them here
    
    # Create intermediate/01_detection directory
    output_dir = auto_create_dirs(args.output, 'intermediate', '01_detection')
    
    cmd = [
        'python', script_path,
        '-i', args.genome_path,
        '-o', output_dir,
        '-m', args.group_info
    ]
    
    # Always add batch parameters since they have defaults
    task_number = args.task_number if args.task_number is not None else 1
    total_tasks = args.total_tasks if args.total_tasks is not None else 1
    cmd.extend(['-number', str(task_number)])
    cmd.extend(['-total', str(total_tasks)])
    
    # Add threads parameter
    if hasattr(args, 'threads') and args.threads:
        cmd.extend(['--threads', str(args.threads)])
    
    run_command(cmd, "HDMI Detect: Finding HGT candidates")
    
    # Check for HGTdetect.py output
    hgt_events_file = os.path.join(output_dir, f'HGT_events_raw_batch_{task_number}.csv')
    
    if os.path.exists(hgt_events_file):
        print(f"\nSuccess! HGT candidates detected: {hgt_events_file}")
    else:
        print(f"\nWarning: Expected output file not found: {hgt_events_file}")
        return
    # Define final output directory (parent of batch directory)
    # When output is intermediate/01_detection, final_output should be the main output directory
    if os.path.basename(output_dir) == '01_detection':
        final_output = os.path.dirname(os.path.dirname(output_dir))
    else:
        final_output = args.output
    

    
    if total_tasks > 1:
        # Multi-task mode: merge batch results
        # Define combined file paths (in final output directory)
        combined_files = {
            'hgt_events': os.path.join(final_output, 'HGT_events_raw.csv'),
            'contig_q': os.path.join(final_output, 'sequences_contig_q.fa'),
            'contig_s': os.path.join(final_output, 'sequences_contig_s.fa'),
            'matched_q': os.path.join(final_output, 'sequences_matched_seq_q.fa'),
            'matched_s': os.path.join(final_output, 'sequences_matched_seq_s.fa')
        }
        
        # Create empty combined files if this is the first batch
        if task_number == 1:
            for file_path in combined_files.values():
                with open(file_path, 'w') as f:
                    pass  # Create empty file
        
        # Define batch file paths (actual files generated by HGTdetect.py)
        batch_files = {
            'hgt_events': os.path.join(output_dir, f'HGT_events_raw_batch_{task_number}.csv'),
            'contig_q': os.path.join(output_dir, f'sequences_contig_q{task_number}.fa'),
            'contig_s': os.path.join(output_dir, f'sequences_contig_s{task_number}.fa'),
            'matched_q': os.path.join(output_dir, f'sequences_matched_seq_q{task_number}.fa'),
            'matched_s': os.path.join(output_dir, f'sequences_matched_seq_s{task_number}.fa')
        }
        
        # Append batch results to combined files
        for file_type, batch_file in batch_files.items():
            if os.path.exists(batch_file):
                with open(batch_file, 'r') as batch_f:
                    content = batch_f.read()
                    if content.strip():  # Only append if file is not empty
                        with open(combined_files[file_type], 'a') as combined_f:
                            if file_type == 'hgt_events':
                                # For CSV files, only append header if it's the first batch
                                if task_number == 1:
                                    combined_f.write(content)
                                else:
                                    # Skip header for subsequent batches
                                    lines = content.split('\n')
                                    if len(lines) > 1:
                                        combined_f.write('\n'.join(lines[1:]) + '\n')
                            else:
                                # For FASTA files, append all content
                                combined_f.write(content)
                print(f"  ✓ Merged {os.path.basename(batch_file)} to combined file")
        
        # Create combined contig file on last batch
        combined_contig_file = os.path.join(final_output, 'sequences_contig_combined.fa')
        if task_number == total_tasks:  # Only on last batch
            with open(combined_contig_file, 'w') as combined:
                # Add sequences from q file
                if os.path.exists(combined_files['contig_q']):
                    with open(combined_files['contig_q'], 'r') as q:
                        combined.write(q.read())
                # Add sequences from s file
                if os.path.exists(combined_files['contig_s']):
                    with open(combined_files['contig_s'], 'r') as s:
                        combined.write(s.read())
            print(f"  ✓ Created combined contig file: {combined_contig_file}")
    else:
        # Single task mode: copy files to final output (remove batch number from filename)
        source_files = {
            'hgt_events': os.path.join(output_dir, f'HGT_events_raw_batch_{task_number}.csv'),
            'contig_q': os.path.join(output_dir, f'sequences_contig_q{task_number}.fa'),
            'contig_s': os.path.join(output_dir, f'sequences_contig_s{task_number}.fa'),
            'matched_q': os.path.join(output_dir, f'sequences_matched_seq_q{task_number}.fa'),
            'matched_s': os.path.join(output_dir, f'sequences_matched_seq_s{task_number}.fa')
        }
        
        target_files = {
            'hgt_events': os.path.join(final_output, 'HGT_events_raw.csv'),
            'contig_q': os.path.join(final_output, 'sequences_contig_q.fa'),
            'contig_s': os.path.join(final_output, 'sequences_contig_s.fa'),
            'matched_q': os.path.join(final_output, 'sequences_matched_seq_q.fa'),
            'matched_s': os.path.join(final_output, 'sequences_matched_seq_s.fa')
        }
        
        # Copy files to final output (removing batch number from filename)
        for file_type, source_file in source_files.items():
            if os.path.exists(source_file):
                shutil.copy2(source_file, target_files[file_type])
                print(f"  ✓ Copied {os.path.basename(source_file)} to {os.path.basename(target_files[file_type])}")
        
        # Create combined contig file with duplicate handling
        combined_contig_file = os.path.join(final_output, 'sequences_contig_combined.fa')
        merge_fasta_files_with_duplicates(
            target_files['contig_q'], 
            target_files['contig_s'], 
            combined_contig_file
        )
        print(f"  ✓ Created combined contig file: {combined_contig_file}")


# Removed cmd_validate function as it's now integrated into cmd_profile


# Removed cmd_merge and cmd_connect functions as they are now integrated into cmd_summary

############ Removed legacy cmd_profile (summary-like) to avoid confusion


def cmd_index(args):
    """HDMI index: Build optimized genome indices based on HGT detection results."""
    import subprocess
    import random
    import pandas as pd
    from Bio import SeqIO
    script_dir = get_script_dir()
    
    print("=== HDMI INDEX: Building optimized genome indices ===")
    print("This builds indices only for contigs involved in HGT events")
    
    # Check required files
    genome_path = args.genome_path
    group_info = args.group_info
    
    if not os.path.exists(genome_path):
        print(f"ERROR: Genome directory not found: {genome_path}")
        sys.exit(1)
    
    if not os.path.exists(group_info):
        print(f"ERROR: Group info file not found: {group_info}")
        sys.exit(1)
    
    # Check if HGT events file exists (from HDMI detect)
    hgt_events_file = None
    if hasattr(args, 'output') and args.output:
        hgt_events_file = os.path.join(args.output, 'HGT_events_raw.csv')
        if not os.path.exists(hgt_events_file):
            # Try to find in intermediate directory
            hgt_events_file = os.path.join(args.output, 'intermediate', '01_detection', 'HGT_events_raw.csv')
    
    if not hgt_events_file or not os.path.exists(hgt_events_file):
        print("ERROR: HGT events file not found. Please run HDMI detect first.")
        print("Expected location: output/HGT_events_raw.csv or output/intermediate/01_detection/HGT_events_raw.csv")
        sys.exit(1)
    
    print(f"Found HGT events file: {hgt_events_file}")
    
    # Read HGT events to extract involved contigs
    print("Extracting contigs involved in HGT events...")
    hgt_df = pd.read_csv(hgt_events_file)
    
    # Add HGT_ID column if it doesn't exist
    if 'HGT_ID' not in hgt_df.columns:
        hgt_df['HGT_ID'] = [f"HGT{i+1}" for i in range(len(hgt_df))]
        # Save updated HGT_events_raw.csv
        hgt_df.to_csv(hgt_events_file, index=False)
        print(f"Added HGT_ID column to HGT_events_raw.csv")
    
    # Extract unique contigs from HGT events
    involved_contigs = set()
    for _, row in hgt_df.iterrows():
        # Add query contig
        query_contig = f"{row['MAG 1'].replace('.fa', '')}_{row['query_congtig_id']}"
        involved_contigs.add(query_contig)
        # Add subject contig  
        subject_contig = f"{row['MAG 2'].replace('.fa', '')}_{row['subject_contig_id']}"
        involved_contigs.add(subject_contig)
    
    print(f"Found {len(involved_contigs)} unique contigs involved in HGT events")
    
    # Read group info and select representatives (using fixed seed for consistency)
    print("Selecting representative genomes...")
    random.seed(42)
    
    group_genomes = {}
    with open(group_info, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) != 2:
                    print(f"ERROR: Line {line_num} has {len(parts)} parts: {repr(line)}")
                    sys.exit(1)
                genome_file, group_id = parts
                genome_name = os.path.splitext(genome_file)[0]
                group_id = int(group_id)
                
                if group_id not in group_genomes:
                    group_genomes[group_id] = []
                group_genomes[group_id].append(genome_name)
    
    # Sort and select representatives
    for group_id in group_genomes:
        group_genomes[group_id].sort()
    
    selected_representatives = {}
    for group_id, genomes in group_genomes.items():
        selected_representatives[group_id] = random.choice(genomes)
    
    print("Selected representative genomes:")
    for group_id, genome in sorted(selected_representatives.items()):
        print(f"  Group {group_id}: {genome}")
    
    # Create index directory in output folder
    if hasattr(args, 'output') and args.output:
        index_dir = auto_create_dirs(args.output, 'index')
    else:
        index_dir = auto_create_dirs(os.path.dirname(genome_path), 'intermediate', 'index')
    print(f"Index files will be stored in: {index_dir}")
    
    # Build representative genome index (for abundance calculation)
    print("\n--- Building Representative Genome Index ---")
    combined_genome_rep = os.path.join(index_dir, 'combined_genome_representatives.fa')
    combined_genome_rep_index = os.path.join(index_dir, 'combined_genome_representatives_index')
    
    # Check if representative index already exists
    if any(os.path.exists(combined_genome_rep_index + ext) for ext in ['.1.bt2', '.1.bt2l']):
        print("Representative genome index already exists. Skipping...")
    else:
        # Combine representative genomes
        with open(combined_genome_rep, 'w') as combined:
            for group_id, genome_name in selected_representatives.items():
                genome_file = genome_name + '.fa'
                genome_file_path = os.path.join(genome_path, genome_file)
                if os.path.exists(genome_file_path):
                    with open(genome_file_path, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                # Add MAG prefix to contig names to avoid duplicates
                                combined.write(f">{genome_name}_{line[1:]}")
                            else:
                                combined.write(line)
        
        # Build index
        cmd_rep = ['bowtie2-build', '--large-index', combined_genome_rep, combined_genome_rep_index]
        if hasattr(args, 'threads') and args.threads:
            cmd_rep.extend(['--threads', str(args.threads)])
        run_command(cmd_rep, "Building representative genome index")
    
    # Build optimized HGT contigs index (only contigs involved in HGT events)
    print("\n--- Building Optimized HGT Contigs Index ---")
    combined_genome_hgt = os.path.join(index_dir, 'combined_genome_hgt_contigs.fa')
    combined_genome_hgt_index = os.path.join(index_dir, 'combined_genome_hgt_contigs_index')
    
    # Check if HGT contigs index already exists
    if any(os.path.exists(combined_genome_hgt_index + ext) for ext in ['.1.bt2', '.1.bt2l']):
        print("HGT contigs index already exists. Skipping...")
    else:
        # Extract only the contigs involved in HGT events
        print(f"Extracting {len(involved_contigs)} contigs from genome files...")
        with open(combined_genome_hgt, 'w') as combined:
            contigs_written = 0
            for genome_file in os.listdir(genome_path):
                if genome_file.endswith(('.fa', '.fasta')):
                    genome_name = os.path.splitext(genome_file)[0]
                    genome_file_path = os.path.join(genome_path, genome_file)
                    
                    # Parse FASTA and extract only involved contigs
                    for record in SeqIO.parse(genome_file_path, "fasta"):
                        contig_id = f"{genome_name}_{record.id}"
                        if contig_id in involved_contigs:
                            combined.write(f">{contig_id}\n{str(record.seq)}\n")
                            contigs_written += 1
        
        print(f"Extracted {contigs_written} contigs to optimized index")
        
        # Build index
        cmd_hgt = ['bowtie2-build', '--large-index', combined_genome_hgt, combined_genome_hgt_index]
        if hasattr(args, 'threads') and args.threads:
            cmd_hgt.extend(['--threads', str(args.threads)])
        run_command(cmd_hgt, "Building optimized HGT contigs index")
    
    # === INTEGRATED CONNECT FUNCTIONALITY ===
    print("\n--- Integrated Connect: Generating Simulated Sequences ---")
    
    # Check if simulated sequences already exist
    simi_sequences_fasta = os.path.join(index_dir, 'simi_sequences.fasta')
    simi_sequences_index = os.path.join(index_dir, 'simi_sequences_index')
    elements_info_raw = os.path.join(index_dir, 'elements_info_raw.csv')
    
    if (os.path.exists(simi_sequences_fasta) and 
        any(os.path.exists(simi_sequences_index + ext) for ext in ['.1.bt2', '.1.bt2l']) and
        os.path.exists(elements_info_raw)):
        print("Simulated sequences and index already exist. Skipping connect step...")
    else:
        # Add script directory to Python path
        script_dir = get_script_dir()
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        
        # Import connect_seq functionality
        from connect_seq import extract_and_simulate_sequences
        
        print("Generating simulated sequences (removing HGT fragments)...")
        
        # Create elements_info by splitting HGT events into query and subject records
        elements_info_rows = []
        
        for idx, row in hgt_df.iterrows():
            # Get HGT event ID from the dataframe
            hgt_event_id = row['HGT_ID']
            
            # Query contig record
            query_contig_id = f"{row['MAG 1'].replace('.fa', '')}_{row['query_congtig_id']}"
            # Build robust element ID without relying on Sequence_* name format
            query_element_id = f"{query_contig_id}_{row['q_start']}_{row['q_end']}"
            elements_info_rows.append({
                'ID': query_element_id,  # Renamed from HGT_ID to ID
                'HGT_ID': hgt_event_id,  # New column for HGT event ID
                'Element_Type': 'query',  # New column for element type
                'contig_id': query_contig_id,
                'start': row['q_start'],
                'end': row['q_end']
            })
            
            # Subject contig record  
            subject_contig_id = f"{row['MAG 2'].replace('.fa', '')}_{row['subject_contig_id']}"
            # Build robust element ID without relying on Sequence_* name format
            subject_element_id = f"{subject_contig_id}_{row['s_start']}_{row['s_end']}"
            elements_info_rows.append({
                'ID': subject_element_id,  # Renamed from HGT_ID to ID
                'HGT_ID': hgt_event_id,  # New column for HGT event ID
                'Element_Type': 'subject',  # New column for element type
                'contig_id': subject_contig_id,
                'start': row['s_start'],
                'end': row['s_end']
            })
        
        # Create DataFrame from the rows
        elements_info = pd.DataFrame(elements_info_rows)
        print(f"Generated {len(elements_info)} elements_info rows")
        if len(elements_info) > 0:
            print(f"First few rows:\n{elements_info.head()}")
        
        # Generate simulated sequences
        extract_and_simulate_sequences(elements_info, combined_genome_hgt, index_dir)
        
        # Build index for simulated sequences
        if os.path.exists(simi_sequences_fasta):
            print("Building bowtie2 index for simulated sequences...")
            cmd_simi = ['bowtie2-build', '--large-index', simi_sequences_fasta, simi_sequences_index]
            if hasattr(args, 'threads') and args.threads:
                cmd_simi.extend(['--threads', str(args.threads)])
            run_command(cmd_simi, "Building simulated sequences index")
            print("✓ Simulated sequences index built successfully")
        else:
            print("✗ Simulated sequences file not found")
    
    print(f"\n✓ Integrated index+connect building completed!")
    print(f"Representative genome index: {combined_genome_rep_index}")
    print(f"Optimized HGT contigs index: {combined_genome_hgt_index}")
    print(f"Simulated sequences index: {simi_sequences_index}")
    print(f"Contigs in HGT index: {len(involved_contigs)}")
    print(f"Elements info raw: {elements_info_raw}")
    print(f"\nUse: HDMI profile [options] for each sample")


# Removed cmd_test function as it's no longer needed

def cmd_connect(args):
    """Connect sequences: Extract HGT sequences and generate simulated sequences"""
    
    # Auto-find HGT events file if not provided
    if not hasattr(args, 'hgt_info') or not args.hgt_info:
        hgt_events = find_file_in_output(args.output, 'HGT_events.csv')
        if hgt_events:
            args.hgt_info = hgt_events
            print(f"Auto-found HGT events file: {args.hgt_info}")
        else:
            print("ERROR: HGT events file not found. Please run HDMI merge first or specify with -i")
            sys.exit(1)
    
    # Auto-find contig sequences file if not provided
    if not hasattr(args, 'contig_seq') or not args.contig_seq:
        contig_seq = find_file_in_output(args.output, 'sequences_contig_combined.fa')
        if contig_seq:
            args.contig_seq = contig_seq
            print(f"Auto-found contig sequences file: {args.contig_seq}")
        else:
            print("ERROR: Contig sequences file not found. Please run HDMI detect first or specify with -s")
            sys.exit(1)
    
    # Create output directory
    connect_output = auto_create_dirs(args.output, 'intermediate', '04_connect')
    
    script_path = os.path.join(get_script_dir(), 'connect_seq.py')
    
    cmd = [
        'python', script_path,
        '-i', args.hgt_info,
        '-s', args.contig_seq,
        '-o', connect_output
    ]
    
    print(f"=== HDMI Connect: Extracting HGT sequences ===")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        print(f"✓ Connect sequences completed successfully")
        
        # Build bowtie2 index for simi_sequences.fasta
        simi_sequences_fasta = os.path.join(connect_output, 'simi_sequences.fasta')
        simi_sequences_index = os.path.join(connect_output, 'simi_sequences_index')
        
        if os.path.exists(simi_sequences_fasta):
            print(f"Building bowtie2 index for simi_sequences.fasta...")
            index_cmd = f"bowtie2-build --threads {args.threads} {simi_sequences_fasta} {simi_sequences_index}"
            try:
                subprocess.run(index_cmd, shell=True, check=True)
                print("✓ Bowtie2 index built successfully")
            except subprocess.CalledProcessError as e:
                print(f"✗ Bowtie2 index building failed: {e}")
                raise
        else:
            print(f"✗ simi_sequences.fasta not found: {simi_sequences_fasta}")
            raise FileNotFoundError(f"simi_sequences.fasta not found: {simi_sequences_fasta}")
            
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Connect sequences failed with return code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        raise


def cmd_profile(args):
    """HGT Profile Analysis: Analyze read coverage for simulated sequences"""
    
    # Auto-extract sample prefix if not provided
    if not hasattr(args, 'prefix') or not args.prefix:
        args.prefix = extract_sample_prefix(args.read1)
        print(f"Auto-extracted sample prefix: {args.prefix}")
    
    # Auto-find HGT table if not provided
    if not hasattr(args, 'hgt_table') or not args.hgt_table:
        hgt_table = find_file_in_output(args.output, 'HGT_events_raw.csv')
        if hgt_table:
            args.hgt_table = hgt_table
            print(f"Auto-found HGT table: {args.hgt_table}")
        else:
            print("ERROR: HGT table not found. Please run HDMI detect first or specify with -t")
            sys.exit(1)
    
    # Set sample_id for compatibility
    args.sample_id = args.prefix
    
    # Auto-create output directory structure
    sample_output = auto_create_dirs(args.output, 'intermediate', '02_validation', args.prefix)
    expected_dirs = ['read_split', 'fraction', 'abundance', 'intermediate', 'temp']
    for dirname in expected_dirs:
        os.makedirs(os.path.join(sample_output, dirname), exist_ok=True)
    
    script_path = os.path.join(get_script_dir(), 'HGTfinder.py')
    
    cmd = [
        'python', script_path,
        '-r1', args.read1,
        '-r2', args.read2,
        '-i', args.prefix,
        '-o', sample_output,
        '-mag_dir', args.genome_path,
        '-table_dir', args.hgt_table,
        '-group_info', args.group_info
    ]
    
    if args.threads:
        cmd.extend(['-t', str(args.threads)])
    if hasattr(args, 'seed') and args.seed:
        cmd.extend(['-seed', str(args.seed)])
    if hasattr(args, 'sth') and args.sth:
        cmd.extend(['-sth', str(args.sth)])
    
    run_command(cmd, f"HDMI Profile: Processing sample {args.sample_id} (integrates validate functionality)")
    
    # Check for expected outputs
    success = True
    for dirname in expected_dirs:
        dirpath = os.path.join(sample_output, dirname)
        if not os.path.exists(dirpath) or not os.listdir(dirpath):
            print(f"WARNING: {dirname} directory is missing or empty")
            success = False
    
    if success:
        print(f"✓ Profile analysis completed successfully for sample {args.sample_id}")
    else:
        print(f"WARNING: Some expected outputs are missing for sample {args.sample_id}")

    
    # Clean up BAM files to save disk space
    cleanup_bam_files(sample_output)

def cmd_summary(args):
    """HDMI summary: Merge results and generate final element table (integrates merge functionality)"""
    
    # Auto-find validation directory if not provided
    if not hasattr(args, 'samples_dir') or not args.samples_dir:
        validation_dir = os.path.join(args.output, 'intermediate', '02_validation')
        if os.path.exists(validation_dir):
            args.samples_dir = validation_dir
            print(f"Auto-found validation directory: {args.samples_dir}")
        else:
            print("ERROR: Validation directory not found. Please run HDMI profile first or specify with -i")
            sys.exit(1)
    
    # Auto-find HGT events file if not provided
    if not hasattr(args, 'hgt_events') or not args.hgt_events:
        hgt_events = find_file_in_output(args.output, 'HGT_events_raw.csv')
        if hgt_events:
            args.hgt_events = hgt_events
            print(f"Auto-found HGT events file: {args.hgt_events}")
        else:
            print("ERROR: HGT events file not found. Please run HDMI detect first or specify with -hgt")
            sys.exit(1)
    
    # Group info must be provided by user
    if not hasattr(args, 'group_info') or not args.group_info:
        print("ERROR: Group info file must be specified with -group")
        sys.exit(1)
    
    # Create output directory and subdirectories
    merge_output = auto_create_dirs(args.output, 'intermediate', '03_final')
    os.makedirs(os.path.join(merge_output, 'merged'), exist_ok=True)
    
    script_path = os.path.join(get_script_dir(), 'merge_and_filter.py')
    
    cmd = [
        'python', script_path,
        '-i', args.samples_dir,
        '-hgt', args.hgt_events,
        '-group', args.group_info,
        '-o', merge_output
    ]
    
    # Use default sample pattern to match all directories
    cmd.extend(['--sample_pattern', '*'])
    
    if args.threshold:
        cmd.extend(['-t', str(args.threshold)])
    if args.temp_dir:
        cmd.extend(['--temp_dir', args.temp_dir])
    
    run_command(cmd, "HDMI Summary: Merging and filtering results (integrates merge functionality)")


def main():
    # Create main parser
    parser = argparse.ArgumentParser(
        prog='HDMI',
        description='HGT Detection from MAGs in Individual',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  HDMI detect -i genomes/ -o output -m groups.txt             # Detect HGT candidates
  HDMI index -g genomes/ -m groups.txt -o output              # Build optimized indices and extract sequences
  HDMI profile -r1 reads_R1.fq -r2 reads_R2.fq --prefix sample1 -o output -g genomes/ -m groups.txt
  HDMI summary -o output -group groups.txt                    # Merge results and generate final summary
        """
    )
    
    # Add subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # HDMI detect
    detect_parser = subparsers.add_parser('detect', help='Detect HGT candidates from MAGs')
    detect_parser.add_argument('-i', '--genome_path', required=True,
                             help='Directory containing genome FASTA files')
    detect_parser.add_argument('-o', '--output', default='result',
                             help='Output directory for detection results (default: result)')
    detect_parser.add_argument('-m', '--group_info', required=True,
                             help='Group info file (Group_info_test.txt format)')
    detect_parser.add_argument('-number', '--task_number', type=int, default=1,
                             help='Task number for parallel processing (1-indexed, default: 1)')
    detect_parser.add_argument('-total', '--total_tasks', type=int, default=1,
                             help='Total number of parallel tasks (default: 1)')
    detect_parser.add_argument('-t', '--threads', type=int, default=1,
                             help='Number of threads for parallel processing (default: 1)')
    detect_parser.add_argument('--count-only', action='store_true',
                             help='Only count genome pairs and estimate performance without running detection')
    
    
    # HDMI index
    index_parser = subparsers.add_parser('index', help='Build optimized genome indices and extract HGT sequences (integrates connect functionality)')
    index_parser.add_argument('-g', '--genome_path', required=True,
                            help='Directory containing genome FASTA files')
    index_parser.add_argument('-m', '--group_info', required=True,
                            help='Group info file')
    index_parser.add_argument('-o', '--output',
                            help='Output directory (default: genome_folder parent/intermediate/index)')
    index_parser.add_argument('-t', '--threads', type=int, default=1,
                            help='Number of threads for bowtie2-build (default: 1)')
    
    # Removed test command as it's no longer needed
    
    
    # HDMI profile
    profile_parser = subparsers.add_parser('profile', help='Validate HGT events and analyze read coverage (integrates validate functionality)')
    profile_parser.add_argument('-r1', '--read1', required=True,
                              help='Path to read1 FASTQ file')
    profile_parser.add_argument('-r2', '--read2', required=True,
                              help='Path to read2 FASTQ file')
    profile_parser.add_argument('--prefix', '--sample_id',
                              help='Sample prefix (auto-extracted from read1 filename if not provided)')
    profile_parser.add_argument('-o', '--output', default='result',
                              help='Output directory (default: result)')
    profile_parser.add_argument('-g', '--genome_path', required=True,
                              help='Directory containing genome FASTA files')
    profile_parser.add_argument('--hgt_table',
                              help='HGT events table (auto-found in output directory if not provided)')
    profile_parser.add_argument('-m', '--group_info', required=True,
                              help='Group info file')
    profile_parser.add_argument('-t', '--threads', type=int, default=1,
                              help='Number of threads (default: 1)')
    profile_parser.add_argument('--seed', type=int, default=42,
                              help='Random seed for reproducibility (default: 42)')
    profile_parser.add_argument('--sth', type=int, default=2,
                              help='Read span threshold (default: 2)')
    
    # HDMI summary
    summary_parser = subparsers.add_parser('summary', help='Merge results and generate final element table (integrates merge functionality)')
    summary_parser.add_argument('-i', '--samples_dir',
                                help='Directory containing validation results (auto-found in output/intermediate/02_validation if not provided)')
    summary_parser.add_argument('-hgt', '--hgt_events',
                                help='HGT events file (auto-found in output directory if not provided)')
    summary_parser.add_argument('-group', '--group_info', required=True,
                                help='Group info file')
    summary_parser.add_argument('-o', '--output', default='result',
                                help='Output directory (default: result)')
    summary_parser.add_argument('--threshold', type=float, default=1.0,
                                help='Abundance threshold (default: 1.0)')
    summary_parser.add_argument('--temp_dir',
                                help='Temporary directory for merged files')
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        print("\nError: No command specified")
        sys.exit(1)
    
    # Execute commands
    try:
        if args.command == 'detect':
            cmd_detect(args)
        elif args.command == 'index':
            cmd_index(args)
        elif args.command == 'profile':
            cmd_profile(args)
        elif args.command == 'summary':
            cmd_summary(args)
        else:
            print(f"ERROR: Unknown command '{args.command}'")
            print("Available commands: detect, index, profile, summary")
            sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {args.command} command failed - {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
