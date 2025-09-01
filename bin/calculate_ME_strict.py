#!/usr/bin/env python3
"""
Calculate ME (Metagenomic Evidence) strict analysis
Process profile analysis results and species abundance to generate existence table

Author: Haoran Peng (penghr21@gmail.com)
"""

import pandas as pd
import glob
import os
import numpy as np
import argparse
from pathlib import Path

def calculate_custom_value(x):
    """Calculate custom value based on read coverage data"""
    if pd.isna(x) or x == '':
        return np.nan
    
    parts = [float(part) for part in str(x).split(",")[:3]] if len(str(x).split(",")) >= 3 else [np.nan, np.nan, np.nan]
    
    if len(parts) == 3:
        # HI = min(HI₁, HI₂) if min(HI₁, HI₂) ≥ 2 else 0
        HI = min(parts[:2]) if min(parts[:2]) >= 2 else 0
        
        # nHI = nHI if nHI ≥ 2 else 0
        nHI = parts[2] if parts[2] >= 2 else 0
        
        # If both HI and nHI are 0, return NA
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
        
    return result

def process_profile_results(profile_dir, species_median_file, output_dir):
    """Process profile results and generate ME strict analysis"""
    
    # If output_dir is a directory, create the output file path
    if os.path.isdir(output_dir):
        output_file = os.path.join(output_dir, 'ME_connect_Process_stricter.csv')
    else:
        output_file = output_dir
    
    # Read species median abundance data
    species_median = pd.read_csv(species_median_file)
    print(f"Loaded species median data: {species_median.shape}")
    
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
            df[f'{sample_name}'] = df[f'{sample_name}'].apply(calculate_custom_value)
            
            if merged_df.empty:
                merged_df = df
            else:
                merged_df = merged_df.merge(df, on='HGT_ID', how='outer')
                
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue
    
    # Merge with species median data
    if not merged_df.empty:
        # Get HGT_ID to species mapping from the original HGT events
        # This would need to be passed as a parameter or derived from the HGT events file
        print(f"Final merged profile data: {merged_df.shape}")
        
        # Save the result
        merged_df.to_csv(output_file, sep=',', index=False)
        print(f"ME strict analysis saved to: {output_file}")
    else:
        print("No data to save!")

def main():
    parser = argparse.ArgumentParser(description='Calculate ME strict analysis from profile results')
    parser.add_argument('-p', '--profile_dir', required=True, help='Directory containing profile results')
    parser.add_argument('-s', '--species_median', required=True, help='Species median abundance file')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    process_profile_results(args.profile_dir, args.species_median, args.output)

if __name__ == "__main__":
    main()
