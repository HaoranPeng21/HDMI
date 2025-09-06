#!/usr/bin/env python3
"""
HDMI: HGT Detection from MAGs in Individual

A comprehensive command-line tool for detecting and validating HGT events.
Author: Haoran Peng (penghr21@gmail.com)
GitHub: https://github.com/HaoranPeng21/HDMI

Usage:
    HDMI index      - Build genome indices for faster processing
    HDMI detect     - Detect HGT candidates from MAGs
    HDMI validate   - Validate HGT events for a single sample
    HDMI merge      - Merge and filter results from multiple samples
    HDMI connect    - Extract HGT sequences and generate simulated sequences
    HDMI profile    - Analyze read coverage for simulated sequences
    HDMI summary    - Generate final element table with metagenomic evidence
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


def cmd_validate(args):
    """HDMI validate: HGT validation for a single sample."""
    
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
    expected_dirs = ['read_split', 'fraction', 'abundance', 'intermediate']
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
    if args.seed:
        cmd.extend(['-seed', str(args.seed)])
    if args.sth:
        cmd.extend(['-sth', str(args.sth)])
    
    run_command(cmd, f"HDMI Validate: Processing sample {args.sample_id}")
    
    # Check for expected outputs
    success = True
    for dirname in expected_dirs:
        dirpath = os.path.join(sample_output, dirname)
        if not os.path.exists(dirpath):
            print(f"Warning: Expected directory not found: {dirpath}")
            success = False
    
    if success:
        print(f"\nSuccess! Sample validation complete: {sample_output}")
        print("Output directories created:")
        for dirname in expected_dirs:
            print(f"  - {dirname}/")
    else:
        print(f"\nWarning: Some expected outputs may be missing")
    
    # Clean up BAM files to save disk space
    cleanup_bam_files(sample_output)


def cmd_merge(args):
    """HDMI merge: Merge and filter results from multiple samples."""
    
    # Auto-find validation directory if not provided
    if not hasattr(args, 'samples_dir') or not args.samples_dir:
        validation_dir = os.path.join(args.output, 'intermediate', '02_validation')
        if os.path.exists(validation_dir):
            args.samples_dir = validation_dir
            print(f"Auto-found validation directory: {args.samples_dir}")
        else:
            print("ERROR: Validation directory not found. Please run HDMI validate first or specify with -i")
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
    
    run_command(cmd, "HDMI Merge: Merging and filtering results")
    
    final_file = os.path.join(merge_output, 'HGT_events.csv')
    if os.path.exists(final_file):
        print(f"\nSuccess! Final HGT events: {final_file}")
        
        # Copy final results to output root directory
        import shutil
        files_to_copy = [
            ('HGT_events.csv', 'HGT_events.csv'),
            ('HGT_events_raw.csv', 'HGT_events_raw.csv')
        ]
        
        for src, dst in files_to_copy:
            src_path = os.path.join(merge_output, src)
            dst_path = os.path.join(args.output, dst)
            if os.path.exists(src_path):
                shutil.copy2(src_path, dst_path)
                print(f"  ✓ Copied {src} to output root")
    else:
        print(f"\nWarning: Expected final output not found: {final_file}")


def cmd_index(args):
    """HDMI index: Pre-build genome indices for faster sample processing."""
    import subprocess
    import random
    script_dir = get_script_dir()
    
    print("=== HDMI INDEX: Pre-building genome indices ===")
    print("This needs to be done only once and will speed up sample processing")
    
    # Check required files
    genome_path = args.genome_path
    group_info = args.group_info
    
    if not os.path.exists(genome_path):
        print(f"ERROR: Genome directory not found: {genome_path}")
        sys.exit(1)
    
    if not os.path.exists(group_info):
        print(f"ERROR: Group info file not found: {group_info}")
        sys.exit(1)
    
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
    
    # Build representative genome index
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
        cmd_rep = ['bowtie2-build', combined_genome_rep, combined_genome_rep_index]
        if hasattr(args, 'threads') and args.threads:
            cmd_rep.extend(['--threads', str(args.threads)])
        run_command(cmd_rep, "Building representative genome index")
            # os.remove(combined_genome_rep)  # Keep temporarily for inspection
    
    # Build all genome index
    print("\n--- Building All Genome Index ---")
    combined_genome_all = os.path.join(index_dir, 'combined_genome_all.fa')
    combined_genome_all_index = os.path.join(index_dir, 'combined_genome_all_index')
    
    # Check if all genome index already exists
    if any(os.path.exists(combined_genome_all_index + ext) for ext in ['.1.bt2', '.1.bt2l']):
        print("All genome index already exists. Skipping...")
    else:
        # Combine all genomes
        with open(combined_genome_all, 'w') as combined:
            for root, _, files in os.walk(genome_path):
                for file in sorted(files):  # Sort for consistency
                    if (file.endswith('.fa') or file.endswith('.fasta')) and file != os.path.basename(combined_genome_all):
                        mag_name = os.path.splitext(file)[0]  # Get mag name without extension
                        with open(os.path.join(root, file), 'r') as f:
                            for line in f:
                                if line.startswith('>'):
                                    # Add MAG prefix to contig names to avoid duplicates
                                    combined.write(f">{mag_name}_{line[1:]}")
                                else:
                                    combined.write(line)
        
        # Build index
        cmd_all = ['bowtie2-build', combined_genome_all, combined_genome_all_index]
        if hasattr(args, 'threads') and args.threads:
            cmd_all.extend(['--threads', str(args.threads)])
        run_command(cmd_all, "Building all genome index")
            # os.remove(combined_genome_all)  # Keep temporarily for inspection
    
    print("\n=== INDEX BUILDING COMPLETE ===")
    print(f"Representative genome index: {combined_genome_rep_index}.*")
    print(f"All genome index: {combined_genome_all_index}.*")
    print("You can now run sample validation much faster!")
    print("Use: HDMI validate [options] for each sample")


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
            index_cmd = f"bowtie2-build {simi_sequences_fasta} {simi_sequences_index}"
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
    if not hasattr(args, 'HGT_table_path') or not args.HGT_table_path:
        hgt_table = find_file_in_output(args.output, 'new_elements_info.csv')
        if hgt_table:
            args.HGT_table_path = hgt_table
            print(f"Auto-found HGT table: {args.HGT_table_path}")
        else:
            print("ERROR: HGT table not found. Please run HDMI connect first or specify with -table_dir")
            sys.exit(1)
    
    # Set sample_id for compatibility
    args.sample_id = args.prefix
    
    # Create output directory
    profile_output = auto_create_dirs(args.output, 'intermediate', '05_profile', args.prefix)
    
    # Copy simi_sequences_index files to output directory if they exist in parent
    table_dir = os.path.dirname(args.HGT_table_path)
    source_index_files = [
        os.path.join(table_dir, 'simi_sequences_index.1.bt2'),
        os.path.join(table_dir, 'simi_sequences_index.2.bt2'),
        os.path.join(table_dir, 'simi_sequences_index.3.bt2'),
        os.path.join(table_dir, 'simi_sequences_index.4.bt2'),
        os.path.join(table_dir, 'simi_sequences_index.rev.1.bt2'),
        os.path.join(table_dir, 'simi_sequences_index.rev.2.bt2')
    ]
    
    # Check if all index files exist
    if all(os.path.exists(f) for f in source_index_files):
        print(f"Copying simi_sequences_index files to {profile_output}...")
        import shutil
        try:
            for source_file in source_index_files:
                target_file = os.path.join(profile_output, os.path.basename(source_file))
                shutil.copy2(source_file, target_file)
            print(f"✓ Index files copied to {profile_output}")
        except Exception as e:
            print(f"✗ Index copying failed: {e}")
            raise
    
    script_path = os.path.join(get_script_dir(), 'HGTprofile.py')
    
    cmd = [
        'python', script_path,
        '-r1', args.read1,
        '-r2', args.read2,
        '-i', args.prefix,
        '-o', profile_output,
        '-table_dir', args.HGT_table_path,
        '-t', str(args.threads)
    ]
    
    print(f"=== HDMI Profile: Analyzing read coverage ===")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        print(f"✓ Profile analysis completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Profile analysis failed with return code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        raise
    
    # Clean up BAM files to save disk space
    cleanup_bam_files(profile_output)

def cmd_summary(args):
    """Summary Analysis: Calculate metagenomic evidence from profile results"""
    
    # Auto-find profile directory if not provided
    if not hasattr(args, 'profile_dir') or not args.profile_dir:
        profile_dir = os.path.join(args.output, 'intermediate', '05_profile')
        if os.path.exists(profile_dir):
            args.profile_dir = profile_dir
            print(f"Auto-found profile directory: {args.profile_dir}")
        else:
            print("ERROR: Profile directory not found. Please run HDMI profile first or specify with -p")
            sys.exit(1)
    
    # Auto-find species median file if not provided
    if not hasattr(args, 'species_median') or not args.species_median:
        species_median = find_file_in_output(args.output, 'species_median.csv')
        if species_median:
            args.species_median = species_median
            print(f"Auto-found species median file: {args.species_median}")
        else:
            print("ERROR: Species median file not found. Please run HDMI merge first or specify with -s")
            sys.exit(1)
    
    script_path = os.path.join(get_script_dir(), 'calculate_ME_strict.py')
    
    cmd = [
        'python', script_path,
        '-p', args.profile_dir,
        '-s', args.species_median,
        '-o', args.output
    ]
    
    print(f"=== HDMI ME Strict: Calculating metagenomic evidence ===")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        print(f"✓ ME strict analysis completed successfully")
        
        # Copy final result to output directory and rename
        import shutil
        me_strict_file = os.path.join(args.output, 'ME_connect_Process_stricter.csv')
        element_table_file = os.path.join(args.output, 'element_table.csv')
        if os.path.exists(me_strict_file):
            shutil.copy2(me_strict_file, element_table_file)
            print(f"  ✓ Copied ME_connect_Process_stricter.csv to element_table.csv")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: ME strict analysis failed with return code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        raise


def main():
    # Create main parser
    parser = argparse.ArgumentParser(
        prog='HDMI',
        description='HGT Detection from MAGs in Individual',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  HDMI index -g genomes/ -m groups.txt -o output              # Pre-build indices (run once)
  HDMI detect -i genomes/ -o output -m groups.txt             # Detect HGT candidates
  HDMI validate -r1 reads_R1.fq -r2 reads_R2.fq --prefix sample1 -o output -g genomes/ -m groups.txt
  HDMI merge -o output -group groups.txt                      # Merge and filter results
  HDMI connect -o output                                       # Extract HGT sequences
  HDMI profile -r1 reads_R1.fq -r2 reads_R2.fq --prefix sample1 -o output
  HDMI summary -o output/element_table.csv                    # Generate final summary
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
    
    # HDMI validate
    validate_parser = subparsers.add_parser('validate', help='Validate HGT events for a single sample')
    validate_parser.add_argument('-r1', '--read1', required=True,
                               help='Path to read1 FASTQ file')
    validate_parser.add_argument('-r2', '--read2', required=True,
                               help='Path to read2 FASTQ file')
    validate_parser.add_argument('--prefix', '--sample_id',
                               help='Sample prefix (auto-extracted from read1 filename if not provided)')
    validate_parser.add_argument('-o', '--output', default='result',
                               help='Output directory (default: result)')
    validate_parser.add_argument('-g', '--genome_path', required=True,
                               help='Directory containing genome FASTA files')
    validate_parser.add_argument('--hgt_table',
                               help='HGT events table (auto-found in output directory if not provided)')
    validate_parser.add_argument('-m', '--group_info', required=True,
                               help='Group info file')
    validate_parser.add_argument('-t', '--threads', type=int, default=1,
                               help='Number of threads (default: 1)')
    validate_parser.add_argument('--seed', type=int, default=42,
                               help='Random seed for reproducibility (default: 42)')
    validate_parser.add_argument('--sth', type=int, default=3,
                               help='Read span threshold (default: 3)')
    
    # HDMI merge
    merge_parser = subparsers.add_parser('merge', help='Merge and filter results from multiple samples')
    merge_parser.add_argument('-i', '--samples_dir',
                            help='Directory containing sample subdirectories (auto-found in output/02_validation if not provided)')
    merge_parser.add_argument('-hgt', '--hgt_events',
                            help='Original HGT events file (auto-found in output directory if not provided)')
    merge_parser.add_argument('-group', '--group_info', required=True,
                            help='Group info file')
    merge_parser.add_argument('-o', '--output', default='result',
                            help='Output directory (default: result)')
    merge_parser.add_argument('--threshold', type=float, default=1.0,
                            help='Abundance threshold (default: 1.0)')
    merge_parser.add_argument('--temp_dir',
                            help='Temporary directory for merged files')
    
    # HDMI index
    index_parser = subparsers.add_parser('index', help='Pre-build genome indices (run once for faster processing)')
    index_parser.add_argument('-g', '--genome_path', required=True,
                            help='Directory containing genome FASTA files')
    index_parser.add_argument('-m', '--group_info', required=True,
                            help='Group info file')
    index_parser.add_argument('-o', '--output',
                            help='Output directory (default: genome_folder parent/intermediate/index)')
    index_parser.add_argument('-t', '--threads', type=int, default=1,
                            help='Number of threads for bowtie2-build (default: 1)')
    
    # Removed test command as it's no longer needed
    
    # HDMI connect
    connect_parser = subparsers.add_parser('connect', help='Extract HGT sequences and generate simulated sequences')
    connect_parser.add_argument('-i', '--hgt_info',
                              help='HGT events information table (auto-found in output directory if not provided)')
    connect_parser.add_argument('-s', '--contig_seq',
                              help='Contig sequence file (auto-found in output directory if not provided)')
    connect_parser.add_argument('-o', '--output', default='./',
                              help='Output directory (default: ./)')
    
    # HDMI profile
    profile_parser = subparsers.add_parser('profile', help='Analyze read coverage for simulated sequences')
    profile_parser.add_argument('-r1', '--read1', required=True,
                              help='Path to read1 FASTQ file')
    profile_parser.add_argument('-r2', '--read2', required=True,
                              help='Path to read2 FASTQ file')
    profile_parser.add_argument('--prefix', '--sample_id',
                              help='Sample prefix (auto-extracted from read1 filename if not provided)')
    profile_parser.add_argument('-o', '--output', default='./',
                              help='Output directory (default: ./)')
    profile_parser.add_argument('-table_dir', '--HGT_table_path',
                              help='Path to HGT table (auto-found in output directory if not provided)')
    profile_parser.add_argument('-t', '--threads', type=int, default=1,
                              help='Number of threads (default: 1)')
    
    # HDMI summary
    summary_parser = subparsers.add_parser('summary', help='Calculate metagenomic evidence from profile results')
    summary_parser.add_argument('--profile_dir',
                                help='Directory containing profile results (auto-found in output/05_profile if not provided)')
    summary_parser.add_argument('-s', '--species_median',
                                help='Species median abundance file (auto-found in output directory if not provided)')
    summary_parser.add_argument('-o', '--output', default='output',
                                help='Output directory (default: output)')
    
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
        elif args.command == 'validate':
            cmd_validate(args)
        elif args.command == 'merge':
            cmd_merge(args)
        elif args.command == 'index':
            cmd_index(args)
        elif args.command == 'connect':
            cmd_connect(args)
        elif args.command == 'profile':
            cmd_profile(args)
        elif args.command == 'summary':
            cmd_summary(args)
    except Exception as e:
        print(f"\nERROR: {args.command} command failed - {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
