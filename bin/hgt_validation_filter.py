#!/usr/bin/env python3
"""
HGT Validation and Filtering Module

This module provides the HGTValidator class for validating and filtering
HGT candidates based on abundance and coverage criteria.

Author: Haoran Peng (penghr21@gmail.com)
"""

import os
import pandas as pd
import numpy as np


class HGTValidator:
    """Class to handle HGT validation and filtering operations."""
    
    def __init__(self, hgt_events_file, group_mapping_file, species_abundance_file,
                 read_split_file, coverage_fraction_file, abundance_threshold=1.0):
        """
        Initialize HGT validator.
        
        Args:
            hgt_events_file: Path to HGT events CSV file
            group_mapping_file: Path to genome-to-group mapping file  
            species_abundance_file: Path to species abundance CSV file
            read_split_file: Path to merged read split CSV file
            coverage_fraction_file: Path to merged coverage fraction CSV file
            abundance_threshold: Minimum abundance threshold for species presence
        """
        self.hgt_events_file = hgt_events_file
        self.group_mapping_file = group_mapping_file
        self.species_abundance_file = species_abundance_file
        self.read_split_file = read_split_file
        self.coverage_fraction_file = coverage_fraction_file
        self.abundance_threshold = abundance_threshold
        
        # Initialize data containers
        self.hgt_df = None
        self.genome_to_cluster = {}
        self.species_abundance_df = None
        self.species_existence_df = None
        self.read_split_df = None
        self.coverage_fraction_df = None
        
    def load_data(self):
        """Load all required data files."""
        print("Loading data files...")
        
        # Load HGT events
        self.hgt_df = pd.read_csv(self.hgt_events_file)
        
        # Add HGT_ID column if it doesn't exist
        if 'HGT_ID' not in self.hgt_df.columns:
            self.hgt_df['HGT_ID'] = [f"HGT_seq_{i+1}" for i in range(len(self.hgt_df))]
            print(f"Added HGT_ID column to {len(self.hgt_df)} HGT events")
        else:
            print(f"Loaded {len(self.hgt_df)} HGT events with existing HGT_ID")
        
        # Load genome to cluster mapping
        if self.group_mapping_file.endswith('.txt'):
            # Handle Group_info_test.txt format: genome_file group_id (space-separated)
            mapping_df = pd.read_csv(self.group_mapping_file, header=None, names=['genome', 'primary_cluster'], delim_whitespace=True)
            mapping_df['genome'] = mapping_df['genome'].str.replace(r'\.(fa|fasta)$', '', regex=True)
        else:
            # Handle CSV format with headers
            mapping_df = pd.read_csv(self.group_mapping_file)
            mapping_df['genome'] = mapping_df['genome'].str.replace(r'\.(fa|fasta)$', '', regex=True)
        
        self.genome_to_cluster = dict(zip(mapping_df['genome'], mapping_df['primary_cluster']))
        print(f"Loaded mapping for {len(self.genome_to_cluster)} genomes")
        
        # Add group information to HGT events
        # Clean MAG names in HGT events to match the mapping dictionary
        self.hgt_df['MAG1_Group'] = self.hgt_df['MAG 1'].str.replace(r'\.(fa|fasta)$', '', regex=True).map(self.genome_to_cluster)
        self.hgt_df['MAG2_Group'] = self.hgt_df['MAG 2'].str.replace(r'\.(fa|fasta)$', '', regex=True).map(self.genome_to_cluster)
        
        # Load species abundance data
        self.species_abundance_df = pd.read_csv(self.species_abundance_file, index_col=0)
        self.species_existence_df = self.species_abundance_df >= self.abundance_threshold
        print(f"Loaded abundance data for {len(self.species_abundance_df)} species across {len(self.species_abundance_df.columns)} samples")
        
        # Load read split data
        self.read_split_df = pd.read_csv(self.read_split_file)
        print(f"Loaded read split data: {self.read_split_df.shape}")
        
        # Load coverage fraction data
        self.coverage_fraction_df = pd.read_csv(self.coverage_fraction_file)
        print(f"Loaded coverage fraction data: {self.coverage_fraction_df.shape}")
        
    def calculate_species_existence_status(self):
        """Calculate species existence status for each HGT event and sample."""
        print("Calculating species existence status...")
        
        sample_columns = self.species_abundance_df.columns
        result_data = {'HGT_ID': self.hgt_df['HGT_ID']}
        
        def calculate_existence_for_sample(row, sample):
            """Calculate existence status for a single sample."""
            mag1_group = row['MAG1_Group']
            mag2_group = row['MAG2_Group']
            
            # Handle missing group mappings
            if pd.isna(mag1_group) or pd.isna(mag2_group):
                return 'NA_NA'
            
            # Look for Group_X format in abundance data
            mag1_group_key = f"Group_{mag1_group}"
            mag2_group_key = f"Group_{mag2_group}"
            
            try:
                mag1_exists = self.species_existence_df.loc[mag1_group_key, sample]
                mag2_exists = self.species_existence_df.loc[mag2_group_key, sample]
                
                if mag1_exists and mag2_exists:
                    return '2_2'
                elif mag1_exists:
                    return '2_0'
                elif mag2_exists:
                    return '0_2'
                else:
                    return '0_0'
            except KeyError:
                # Try without Group_ prefix as fallback
                try:
                    mag1_exists = self.species_existence_df.loc[str(mag1_group), sample]
                    mag2_exists = self.species_existence_df.loc[str(mag2_group), sample]
                    
                    if mag1_exists and mag2_exists:
                        return '2_2'
                    elif mag1_exists:
                        return '2_0'
                    elif mag2_exists:
                        return '0_2'
                    else:
                        return '0_0'
                except KeyError:
                    return 'NA_NA'
        
        for sample in sample_columns:
            result_data[sample] = self.hgt_df.apply(
                lambda row: calculate_existence_for_sample(row, sample), axis=1
            )
        
        self.species_existence_status = pd.DataFrame(result_data)
        return self.species_existence_status
    
    def process_read_split(self):
        """Process read split data."""
        print("Processing read split data...")
        
        # Read split data is already in correct format (X_Y), no conversion needed
        return self.read_split_df.copy()
    
    def process_coverage_fraction(self):
        """Process coverage fraction data with threshold filtering."""
        print("Processing coverage fraction data...")
        
        def check_coverage_values(value):
            """Check if coverage values meet threshold (90%)."""
            try:
                num1, num2 = map(int, value.split('_'))
                if num1 >= 90 or num2 >= 90:
                    return '2_2'
                else:
                    return '0_0'
            except (ValueError, AttributeError):
                return '0_0'
        
        # Apply threshold check to all sample columns (skip HGT_ID column)
        processed_df = self.coverage_fraction_df.copy()
        sample_cols = [col for col in processed_df.columns if col != 'HGT_ID']
        processed_df[sample_cols] = processed_df[sample_cols].applymap(check_coverage_values)
        
        return processed_df
    
    def validate_read_split_with_coverage(self, processed_read_split_df, processed_coverage_df):
        """Validate read split results against coverage fraction."""
        print("Validating read split with coverage...")
        
        def validate_values(loose_val, coverage_val):
            """Validate loose read split against coverage fraction."""
            if loose_val == '2_2':
                if coverage_val == '2_2':
                    return loose_val
                elif coverage_val == '2_0':
                    return coverage_val
                elif coverage_val == '0_2':
                    return coverage_val
                else:
                    return '0_0'
            elif loose_val == '2_0':
                return loose_val if coverage_val.startswith('2_') else '0_0'
            elif loose_val == '0_2':
                return loose_val if coverage_val.endswith('_2') else '0_0'
            else:
                return loose_val
        
        validated_df = processed_read_split_df.copy()
        
        # Validate each sample column
        for col in processed_read_split_df.columns:
            if col != 'HGT_ID':
                validated_df[col] = [
                    validate_values(
                        processed_read_split_df.at[idx, col],
                        processed_coverage_df.at[idx, col]
                    )
                    for idx in range(len(processed_read_split_df))
                ]
        
        return validated_df
    
    def validate_with_species_existence(self, validated_df):
        """Final validation against species existence."""
        print("Validating with species existence...")
        
        def validate_final(loose_val, existence_val):
            """Final validation logic."""
            if loose_val == '2_2':
                if existence_val == '2_2':
                    return loose_val
                elif existence_val == '2_0':
                    return "2_NA"
                elif existence_val == '0_2':
                    return "NA_2"
                else:
                    return 'NA_NA'
            elif loose_val == '2_0':
                if existence_val == '2_2':
                    return "2_0"
                elif existence_val == '2_0':
                    return "2_NA"
                elif existence_val == '0_2':
                    return "NA_0"
                else:
                    return 'NA_NA'
            elif loose_val == '0_2':
                if existence_val == '2_2':
                    return "0_2"
                elif existence_val == '2_0':
                    return "0_NA"
                elif existence_val == '0_2':
                    return "NA_2"
                else:
                    return 'NA_NA'
            else:
                if existence_val == '2_2':
                    return "0_0"
                elif existence_val == '2_0':
                    return "0_NA"
                elif existence_val == '0_2':
                    return "NA_0"
                elif existence_val == '0_0':
                    return 'NA_NA'
                else:
                    return loose_val
        
        final_validated_df = validated_df.copy()
        
        # Final validation for each sample column
        for col in validated_df.columns:
            if col != 'HGT_ID':
                final_validated_df[col] = [
                    validate_final(
                        validated_df.at[idx, col],
                        self.species_existence_status.at[idx, col]
                    )
                    for idx in range(len(validated_df))
                ]
        
        return final_validated_df
    
    def filter_hgt_events(self, final_validated_df):
        """Filter HGT events based on validation results."""
        print("Filtering HGT events...")
        
        # Get data columns (exclude HGT_ID)
        data_cols = [col for col in final_validated_df.columns if col != 'HGT_ID']
        
        # Count occurrences of each validation status
        unique_values = final_validated_df[data_cols].stack().unique()
        counts_df = pd.DataFrame()
        
        for value in unique_values:
            counts_df[value] = final_validated_df[data_cols].apply(
                lambda row: sum(row.values == value), axis=1
            )
        
        # Add HGT_ID column
        counts_df.insert(0, 'HGT_ID', final_validated_df['HGT_ID'])
        
        # Apply filtering conditions (handle missing columns safely)
        def get_column_safe(df, col_name):
            if col_name in df.columns:
                return df[col_name]
            else:
                return pd.Series([0] * len(df), index=df.index)
        
        condition1 = (get_column_safe(counts_df, '2_NA') + get_column_safe(counts_df, '2_0')) > 0
        condition2 = (get_column_safe(counts_df, 'NA_2') + get_column_safe(counts_df, '0_2')) > 0
        condition3 = get_column_safe(counts_df, '2_2') > 0
        
        # Apply strict filtering conditions
        filtered_df = counts_df.loc[(condition1 & condition2) | condition3]
        
        # Get filtered HGT list
        hgt_list = filtered_df['HGT_ID'].tolist()
        final_validation_results = final_validated_df[final_validated_df["HGT_ID"].isin(hgt_list)]
        filtered_hgt_events = self.hgt_df[self.hgt_df["HGT_ID"].isin(hgt_list)]
        
        print(f"Filtered from {len(self.hgt_df)} to {len(filtered_hgt_events)} HGT events")
        
        return final_validation_results, filtered_hgt_events, counts_df
    
    def run_validation_pipeline(self, output_dir):
        """Run the complete validation pipeline."""
        print("=== Starting HGT Validation Pipeline ===")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Load all data
        self.load_data()
        
        # Calculate species existence status
        species_existence = self.calculate_species_existence_status()
        species_existence.to_csv(os.path.join(output_dir, "species_existence.csv"), index=False)
        
        # Process read split data
        processed_read_split = self.process_read_split()
        
        # Process coverage fraction
        processed_coverage = self.process_coverage_fraction()
        
        # Validate read split with coverage
        validated_read_split = self.validate_read_split_with_coverage(processed_read_split, processed_coverage)
        
        # Final validation with species existence
        final_validated = self.validate_with_species_existence(validated_read_split)
        
        # Filter HGT events
        validation_results, filtered_events, count_stats = self.filter_hgt_events(final_validated)
        
        # Save results
        validation_results.to_csv(os.path.join(output_dir, "HGT_Present_Table_Filter.csv"), index=False)
        filtered_events.to_csv(os.path.join(output_dir, "HGT_events.csv"), index=False)
        count_stats.to_csv(os.path.join(output_dir, "validation_count_stats.csv"), index=False)
        
        print(f"=== Validation Pipeline Complete ===")
        print(f"Results saved in: {output_dir}")
        print(f"Final filtered HGT events: {len(filtered_events)}")
        
        return validation_results, filtered_events


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='HGT validation and filtering pipeline')
    
    parser.add_argument('-hgt', '--hgt_events', required=True,
                      help='HGT events CSV file from HGTdetect.py')
    parser.add_argument('-group', '--group_mapping', required=True,
                      help='Genome to group mapping file (CSV or Group_info_test.txt format)')
    parser.add_argument('-abun', '--species_abundance', required=True,
                      help='Species abundance CSV file (species_median.csv)')
    parser.add_argument('-split', '--read_split', required=True,
                      help='Merged read split CSV file (HGT_split_loose.csv)')
    parser.add_argument('-frac', '--coverage_fraction', required=True,
                      help='Merged coverage fraction CSV file (HGT_cover_fraction.csv)')
    parser.add_argument('-o', '--output_dir', required=True,
                      help='Output directory for validation results')
    parser.add_argument('-t', '--threshold', default=1.0, type=float,
                      help='Abundance threshold for species presence (default: 1.0)')
    
    args = parser.parse_args()
    
    # Initialize validator
    validator = HGTValidator(
        hgt_events_file=args.hgt_events,
        group_mapping_file=args.group_mapping,
        species_abundance_file=args.species_abundance,
        read_split_file=args.read_split,
        coverage_fraction_file=args.coverage_fraction,
        abundance_threshold=args.threshold
    )
    
    # Run validation pipeline
    validation_results, filtered_events = validator.run_validation_pipeline(args.output_dir)


if __name__ == '__main__':
    main()
