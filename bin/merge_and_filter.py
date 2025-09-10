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
import numpy as np
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


def calculate_custom_value(x, species_existence_status=None, coverage_fraction_status=None, is_query_element=None):
    """Calculate custom value based on read coverage data with validation filters"""
    if pd.isna(x) or x == '':
        return np.nan
    
    parts = [float(part) for part in str(x).split(",")[:3]] if len(str(x).split(",")) >= 3 else [np.nan, np.nan, np.nan]
    
    if len(parts) == 3:
        # HI = min(HI₁, HI₂) if min(HI₁, HI₂) ≥ 2 else 0
        HI = min(parts[:2]) if min(parts[:2]) >= 2 else 0
        
        # nHI = nHI if nHI ≥ 2 else 0
        nHI = parts[2] if parts[2] >= 2 else 0
        
        # If both HI and nHI are 0, return NA (reads ratio not sufficient)
        if HI == 0 and nHI == 0:
            result = np.nan
        # HGT presence = HI / (HI + nHI) > 0
        elif HI / (HI + nHI) > 0:
            result = 1.0
        # HGT absence = nHI / (HI + nHI) = 0 (when HI = 0 and nHI > 0)
        else:
            result = 0.0
    else:
        result = np.nan
    
    # Apply validation filters if provided
    if result == 1.0 and (species_existence_status is not None or coverage_fraction_status is not None):
        if is_query_element is not None:
            # Convert to string if needed, handle NaN
            species_status = str(species_existence_status) if species_existence_status is not None and not pd.isna(species_existence_status) else None
            coverage_status = str(coverage_fraction_status) if coverage_fraction_status is not None and not pd.isna(coverage_fraction_status) else None
            
            # For query elements: check first MAG and first fraction
            if is_query_element:
                # Check if first MAG exists (species_existence_status starts with '2_')
                if species_status is not None and not species_status.startswith('2_'):
                    result = np.nan
                # Check if first fraction is sufficient (coverage_fraction_status starts with '2_')
                elif coverage_status is not None and not coverage_status.startswith('2_'):
                    result = np.nan
            else:
                # For subject elements: check second MAG and second fraction
                # Check if second MAG exists (species_existence_status ends with '_2')
                if species_status is not None and not species_status.endswith('_2'):
                    result = np.nan
                # Check if second fraction is sufficient (coverage_fraction_status ends with '_2')
                elif coverage_status is not None and not coverage_status.endswith('_2'):
                    result = np.nan
        
    return result


def run_me_strict_analysis(profile_dir, species_median_file, output_file, filtered_hgt_file=None, validation_dir=None):
    """Run ME strict analysis from profile results"""
    
    # Read species median abundance data
    species_median = pd.read_csv(species_median_file)
    print(f"Loaded species median data: {species_median.shape}")
    
    # Load validation data if available
    species_existence_df = None
    coverage_fraction_df = None
    if validation_dir and os.path.exists(validation_dir):
        species_existence_file = os.path.join(validation_dir, 'species_existence.csv')
        coverage_fraction_file = os.path.join(validation_dir, 'merged', 'HGT_cover_fraction.csv')
        
        if os.path.exists(species_existence_file):
            species_existence_df = pd.read_csv(species_existence_file)
            print(f"Loaded species existence data: {species_existence_df.shape}")
        
        if os.path.exists(coverage_fraction_file):
            coverage_fraction_df = pd.read_csv(coverage_fraction_file)
            print(f"Loaded coverage fraction data: {coverage_fraction_df.shape}")
    
    # Read new_elements_info.csv to get element mapping
    filtered_hgt_ids = set()
    element_to_hgt_mapping = {}  # Maps element ID to (HGT_ID, Element_Type)
    elements_info_raw_file = os.path.join(os.path.dirname(os.path.dirname(profile_dir)), 'index', 'elements_info_raw.csv')
    
    if os.path.exists(elements_info_raw_file):
        elements_info_df = pd.read_csv(elements_info_raw_file)
        print(f"Loaded elements info: {elements_info_df.shape}")
        
        # Read filtered HGT events if provided
        if filtered_hgt_file and os.path.exists(filtered_hgt_file):
            filtered_hgt_df = pd.read_csv(filtered_hgt_file)
            # Check if file has data (more than just header)
            if len(filtered_hgt_df) > 0:
                filtered_hgt_event_ids = set(filtered_hgt_df['HGT_ID'].tolist())
                print(f"Filtered HGT event IDs: {len(filtered_hgt_event_ids)}")
                
                # Filter elements_info to only include filtered HGT events
                elements_info_df = elements_info_df[elements_info_df['HGT_ID'].isin(filtered_hgt_event_ids)]
                print(f"Filtered elements info: {elements_info_df.shape}")
            else:
                print("Filtered HGT events file is empty, using all elements for element-specific validation")
        else:
            print("No filtered HGT events file found, using all elements for element-specific validation")
        
        # Create mapping for all elements (for element-specific validation)
        for _, row in elements_info_df.iterrows():
            element_id = row['ID']
            hgt_event_id = row['HGT_ID']
            element_type = row['Element_Type']
            filtered_hgt_ids.add(element_id)
            element_to_hgt_mapping[element_id] = (hgt_event_id, element_type)
        
        print(f"Loaded {len(filtered_hgt_ids)} filtered HGT element IDs")
        print(f"Built element mapping for {len(element_to_hgt_mapping)} elements")
        print(f"Sample element mapping: {list(element_to_hgt_mapping.items())[:3]}")
    else:
        print("No elements_info_raw.csv found, processing all elements")
    
    # Find all profile result files
    profile_files = glob.glob(os.path.join(profile_dir, "*/temp/*.tsv"))
    print(f"Found {len(profile_files)} profile result files")
    
    if not profile_files:
        print("No profile result files found!")
        return
    
    merged_df = pd.DataFrame()
    
    for file in profile_files:
        base_name = os.path.basename(file)
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
        column_prefix = base_name.replace("_reads_counts_very_all.tsv", "")
        
        print(f"Processing {sample_name}: {base_name}")
        
        try:
            df = pd.read_csv(file, sep='\t', header=None, dtype=str)
            df.columns = ['HGT_ID'] + [f'{sample_name}']  # Use only sample_name, no duplication
            
            # Filter to only include elements that have mapping information
            if element_to_hgt_mapping:
                df = df[df['HGT_ID'].isin(element_to_hgt_mapping.keys())]
                print(f"Filtered to {len(df)} elements from {len(element_to_hgt_mapping)} mapped elements")
            
            # Apply validation filters if available
            if species_existence_df is not None and coverage_fraction_df is not None and element_to_hgt_mapping:
                # Find the correct column name in validation data (usually LLD_sample_name)
                validation_col = None
                for col in species_existence_df.columns:
                    if col != 'HGT_ID' and sample_name in col:
                        validation_col = col
                        break
                
                if validation_col:
                    # Add HGT_Event_ID and Element_Type to df
                    df['HGT_Event_ID'] = df['HGT_ID'].map(lambda x: element_to_hgt_mapping.get(x, (None, None))[0])
                    df['Element_Type'] = df['HGT_ID'].map(lambda x: element_to_hgt_mapping.get(x, (None, None))[1])
                    
                    # For species_existence, use HGT_Event_ID directly (HGT1, HGT2, ...)
                    # For coverage_fraction, we need to convert HGT_Event_ID to HGT_seq_X format
                    df['HGT_seq_ID'] = df['HGT_Event_ID'].map(lambda x: f"HGT_seq_{x[3:]}" if x and x.startswith('HGT') else None)
                    
                    # Merge with species existence data using HGT_Event_ID
                    df = df.merge(species_existence_df[['HGT_ID', validation_col]], left_on='HGT_Event_ID', right_on='HGT_ID', how='left', suffixes=('', '_existence'))
                    
                    # Merge with coverage fraction data using HGT_seq_ID
                    df = df.merge(coverage_fraction_df[['HGT_ID', validation_col]], left_on='HGT_seq_ID', right_on='HGT_ID', how='left', suffixes=('', '_coverage'))
                    
                    # Apply calculate_custom_value with validation filters
                    df[f'{sample_name}'] = df.apply(
                        lambda row: calculate_custom_value(
                            row[f'{sample_name}'], 
                            row.get(f'{validation_col}_existence'),
                            row.get(f'{validation_col}_coverage'),
                            row.get('Element_Type') == 'query'
                        ), axis=1
                    )
                    
                    # Drop validation columns
                    df = df.drop(columns=['HGT_Event_ID', 'Element_Type', 'HGT_seq_ID', f'{validation_col}_existence', f'{validation_col}_coverage'], errors='ignore')
                else:
                    print(f"Warning: Could not find validation column for sample {sample_name}")
                    # Apply calculate_custom_value without validation filters
                    df[f'{sample_name}'] = df[f'{sample_name}'].apply(calculate_custom_value)
            else:
                # Apply calculate_custom_value without validation filters
                df[f'{sample_name}'] = df[f'{sample_name}'].apply(calculate_custom_value)
            
            # Group by HGT_ID and take the maximum value (presence if any contig shows presence)
            # This handles the case where each HGT_ID has multiple rows (query and subject)
            df_grouped = df.groupby('HGT_ID')[f'{sample_name}'].max().reset_index()
            
            if merged_df.empty:
                merged_df = df_grouped
            else:
                merged_df = merged_df.merge(df_grouped, on='HGT_ID', how='outer')
                
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue
    
    # Save the result
    if not merged_df.empty:
        print(f"Final merged profile data: {merged_df.shape}")
        merged_df.to_csv(output_file, sep=',', index=False)
        print(f"Element table saved to: {output_file}")
    else:
        print("No data to save!")


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
        
        # === INTEGRATED SUMMARY ANALYSIS ===
        print("\n=== INTEGRATED SUMMARY ANALYSIS ===")
        try:
            # Find profile directory (should be in validation samples)
            profile_dir = args.samples_dir
            species_median_file = os.path.join(merged_dir, 'species_median.csv')
            
            if os.path.exists(species_median_file):
                print("Processing profile results for ME strict analysis...")
                element_table_file = os.path.join(final_output_dir, 'element_table.csv')
                filtered_hgt_file = os.path.join(final_output_dir, 'HGT_events.csv')
                run_me_strict_analysis(profile_dir, species_median_file, element_table_file, filtered_hgt_file, final_output_dir)
                print(f"Element table saved to: {element_table_file}")
                
                # Also copy element table to main output directory as final result
                main_output_dir = os.path.dirname(os.path.dirname(final_output_dir))  # Go up two levels from 03_final to output
                main_element_table = os.path.join(main_output_dir, 'element_table.csv')
                shutil.copy2(element_table_file, main_element_table)
                print(f"Element table also saved to main output: {main_element_table}")

                # Copy filtered HGT events to main output directory (overwrite any pre-profile file)
                main_hgt_events = os.path.join(main_output_dir, 'HGT_events.csv')
                if os.path.exists(filtered_hgt_file):
                    shutil.copy2(filtered_hgt_file, main_hgt_events)
                    print(f"Filtered HGT events also saved to main output: {main_hgt_events}")
                
                # Create final elements_info.csv with element details
                elements_info_raw_file = os.path.join(os.path.dirname(os.path.dirname(profile_dir)), 'index', 'elements_info_raw.csv')
                if os.path.exists(elements_info_raw_file):
                    elements_info_df = pd.read_csv(elements_info_raw_file)
                    # Filter to only include elements that have evidence in element_table
                    element_table_df = pd.read_csv(element_table_file, index_col=0)
                    filtered_elements = elements_info_df[elements_info_df['ID'].isin(element_table_df.index)]
                    main_elements_info = os.path.join(main_output_dir, 'element_info.csv')
                    filtered_elements.to_csv(main_elements_info, index=False)
                    print(f"Element info also saved to main output: {main_elements_info}")
                else:
                    print("Warning: elements_info_raw.csv not found, skipping elements_info.csv creation")
            else:
                print(f"Warning: Species median file not found: {species_median_file}")
                print("Skipping ME strict analysis...")
        except Exception as e:
            print(f"Warning: ME strict analysis failed: {e}")
            print("Continuing without element table generation...")
        
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
