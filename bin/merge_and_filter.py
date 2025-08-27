#!/usr/bin/env python3
"""
Merge Sample Results and Filter HGT Events

This script combines Step 3 (merging) and Step 4 (validation & filtering):
1. Merges results from multiple samples
2. Validates and filters HGT events
3. Outputs final filtered HGT events

Author: Haoran Peng (penghr21@gmail.com)
"""

import os
import sys
import argparse
import pandas as pd
import glob
import subprocess
import shutil
from pathlib import Path


def merge_csv_files(file_pattern, output_file, merge_on='HGT_ID'):
    """
    Merge CSV files matching a pattern.
    
    Args:
        file_pattern: Glob pattern to find CSV files
        output_file: Output merged CSV file path
        merge_on: Column to merge on (default: 'HGT_ID')
        
    Returns:
        Merged DataFrame
    """
    csv_files = glob.glob(file_pattern)
    
    if not csv_files:
        print(f"Warning: No files found matching pattern: {file_pattern}")
        return None
    
    print(f"Found {len(csv_files)} files to merge:")
    for f in csv_files:
        print(f"  - {f}")
    
    # Read first file as base
    merged_df = pd.read_csv(csv_files[0])
    
    # Merge with remaining files
    for csv_file in csv_files[1:]:
        df = pd.read_csv(csv_file)
        merged_df = pd.merge(merged_df, df, on=merge_on, how='outer')
    
    # Save merged file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved: {output_file}")
    
    return merged_df


def merge_sample_results(samples_base_dir, output_dir, sample_pattern="sample_*"):
    """Merge results from all samples."""
    print("=== STEP 3: MERGING SAMPLE RESULTS ===")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge read_split files (loose)
    print("\nMerging read_split loose files...")
    loose_pattern = os.path.join(samples_base_dir, sample_pattern, "read_split", "*_loose.csv")
    loose_output = os.path.join(output_dir, "HGT_split_loose.csv")
    merge_csv_files(loose_pattern, loose_output)
    
    # Merge read_split files (strict)
    print("\nMerging read_split strict files...")
    strict_pattern = os.path.join(samples_base_dir, sample_pattern, "read_split", "*_strict.csv")
    strict_output = os.path.join(output_dir, "HGT_split_strict.csv")
    merge_csv_files(strict_pattern, strict_output)
    
    # Merge fraction files
    print("\nMerging fraction files...")
    fraction_pattern = os.path.join(samples_base_dir, sample_pattern, "fraction", "*_fraction.csv")
    fraction_output = os.path.join(output_dir, "HGT_cover_fraction.csv")
    merge_csv_files(fraction_pattern, fraction_output)
    
    # Merge abundance files
    print("\nMerging abundance files...")
    abundance_pattern = os.path.join(samples_base_dir, sample_pattern, "abundance", "*_median.csv")
    abundance_output = os.path.join(output_dir, "species_median.csv")
    merge_csv_files(abundance_pattern, abundance_output, merge_on='Group_id')
    
    # Check if all expected files were created
    expected_files = [
        'HGT_split_loose.csv',
        'HGT_cover_fraction.csv',
        'species_median.csv'
    ]
    
    missing_files = []
    for filename in expected_files:
        filepath = os.path.join(output_dir, filename)
        if not os.path.exists(filepath):
            missing_files.append(filename)
    
    if missing_files:
        raise FileNotFoundError(f"Missing merged files: {missing_files}")
    
    print(f"\nMerging complete. Files saved in: {output_dir}")
    return output_dir


def run_hgt_validation(hgt_events_file, group_info, merged_dir, output_dir, abundance_threshold=1.0):
    """Run HGT validation and filtering."""
    print("\n=== STEP 4: HGT VALIDATION AND FILTERING ===")
    
    # Get current script directory to find hgt_validation_filter.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    validation_script = os.path.join(script_dir, 'hgt_validation_filter.py')
    
    if not os.path.exists(validation_script):
        raise FileNotFoundError(f"Validation script not found: {validation_script}")
    
    # Define input file paths
    species_abundance = os.path.join(merged_dir, 'species_median.csv')
    read_split = os.path.join(merged_dir, 'HGT_split_loose.csv')
    coverage_fraction = os.path.join(merged_dir, 'HGT_cover_fraction.csv')
    
    # Check if required files exist
    required_files = [species_abundance, read_split, coverage_fraction]
    for filepath in required_files:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Required file not found: {filepath}")
    
    # Build command
    cmd = [
        'python', validation_script,
        '-hgt', hgt_events_file,
        '-group', group_info,
        '-abun', species_abundance,
        '-split', read_split,
        '-frac', coverage_fraction,
        '-o', output_dir,
        '-t', str(abundance_threshold)
    ]
    
    print(f"Running validation command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Validation failed with return code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        raise
    
    # Check for final output
    final_hgt_file = os.path.join(output_dir, 'HGT_events.csv')
    if not os.path.exists(final_hgt_file):
        raise FileNotFoundError(f"Final HGT events file not found: {final_hgt_file}")
    
    print(f"Validation complete. Final results: {final_hgt_file}")
    return final_hgt_file


def main():
    parser = argparse.ArgumentParser(description='Merge sample results and filter HGT events (Step 3 + 4)')
    
    parser.add_argument('-i', '--samples_dir', required=True,
                      help='Directory containing sample subdirectories from Step 2')
    parser.add_argument('-hgt', '--hgt_events', required=True,
                      help='HGT events CSV file from Step 1 (HGTdetect.py)')
    parser.add_argument('-group', '--group_info', required=True,
                      help='Group info file (Group_info_test.txt format)')
    parser.add_argument('-o', '--output_dir', required=True,
                      help='Output directory for final results')
    parser.add_argument('--sample_pattern', default='sample_*',
                      help='Pattern to match sample directories (default: sample_*)')
    parser.add_argument('-t', '--threshold', default=1.0, type=float,
                      help='Abundance threshold for species presence (default: 1.0)')
    parser.add_argument('--temp_dir', 
                      help='Directory for temporary merged files (default: output_dir/merged)')
    
    args = parser.parse_args()
    
    # Set up directories with proper structure
    if args.temp_dir:
        merged_dir = args.temp_dir
    else:
        merged_dir = os.path.join(args.output_dir, 'merged')
    
    # Create proper directory structure
    final_output_dir = args.output_dir
    os.makedirs(final_output_dir, exist_ok=True)
    os.makedirs(merged_dir, exist_ok=True)
    
    print("=== MERGE AND FILTER PIPELINE ===")
    print(f"Samples directory: {args.samples_dir}")
    print(f"HGT events file: {args.hgt_events}")
    print(f"Group info: {args.group_info}")
    print(f"Temp merged dir: {merged_dir}")
    print(f"Final output: {final_output_dir}")
    print(f"Abundance threshold: {args.threshold}")
    
    try:
        # Step 3: Merge sample results
        merge_sample_results(args.samples_dir, merged_dir, args.sample_pattern)
        
        # Step 4: Validate and filter HGT events
        final_hgt_file = run_hgt_validation(
            hgt_events_file=args.hgt_events,
            group_info=args.group_info,
            merged_dir=merged_dir,
            output_dir=final_output_dir,
            abundance_threshold=args.threshold
        )
        
        # Copy original HGT events to final results as raw file
        raw_hgt_file = os.path.join(final_output_dir, 'HGT_events_raw.csv')
        shutil.copy2(args.hgt_events, raw_hgt_file)
        print(f"Original HGT events copied to: {raw_hgt_file}")
        
        print(f"\n=== PIPELINE COMPLETED SUCCESSFULLY ===")
        print(f"Original HGT events: {raw_hgt_file}")
        print(f"Final filtered HGT events: {final_hgt_file}")
        print(f"Merged intermediate files: {merged_dir}")
        print(f"Final results directory: {final_output_dir}")
        
    except Exception as e:
        print(f"\nERROR: Pipeline failed - {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
