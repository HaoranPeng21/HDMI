#!/usr/bin/env python3
"""
Connect Sequences Script

Extract HGT sequences and generate simulated sequences.

Author: Haoran Peng (penghr21@gmail.com)
GitHub: https://github.com/HaoranPeng21/HDMI
"""

import pandas as pd
import numpy as np
import importlib.util
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_positions_from_hgt_id(hgt_id, row):
    # For HGT_seq_X format, extract contig name and positions from the row
    if hgt_id.startswith('HGT_seq_'):
        # Use the contig_id and start/end positions from the row
        contig_name = row['contig_id']
        start = row['start']
        end = row['end']
        return contig_name, start, end
    else:
        # Original format: contig_name_start_end
        parts = hgt_id.rsplit('_', 2)  # Split from right, expecting 3 parts
        if len(parts) != 3:
            raise ValueError(f"Unexpected HGT_ID format: {hgt_id}")
        contig_name = parts[0]
        start, end = sorted(map(int, [parts[1], parts[2]]))
        return contig_name, start, end


def extract_and_simulate_sequences(elements_info_path, contigs_path, output_dir):
    #elements_info = pd.read_csv(elements_info_path)
    elements_info = elements_info_path
    contigs = SeqIO.to_dict(SeqIO.parse(contigs_path, "fasta"))

    new_info_rows = []
    sequences = {}  # Dictionary to hold sequences to be written
    added_contigs = set()

    for index, row in elements_info.iterrows():
        hgt_id = row['HGT_ID']
        contig_name, start, end = extract_positions_from_hgt_id(hgt_id, row)

        if contig_name not in added_contigs:
            original_seq = contigs[contig_name].seq
            original_seq_record = SeqRecord(Seq(str(original_seq)), id=contig_name, description="")
            sequences[contig_name] = original_seq_record
            added_contigs.add(contig_name)

        simulated_seq = original_seq[:start] + original_seq[end:]
        simulated_id = f"{hgt_id}_sim"
        simulated_seq_record = SeqRecord(Seq(str(simulated_seq)), id=simulated_id, description="")
        sequences[simulated_id] = simulated_seq_record  # This allows for multiple simulated sequences per contig
        
        # Update information for CSV
        new_info_rows.append({
            'HGT_ID': hgt_id,
            'contig_name': contig_name,
            'HGT_ID_sim': simulated_id,
            'contig_start': start,
            'contig_end': end,
            'HGT_ID_sim_position': start
        })

    # Write all sequences to a single FASTA file
    with open(f"{output_dir}/simi_sequences.fasta", "w") as combined_fasta:
        SeqIO.write(sequences.values(), combined_fasta, "fasta")

    new_info_df = pd.DataFrame(new_info_rows)
    new_info_df.to_csv(f"{output_dir}/new_elements_info.csv", index=False)

def check_query_presence(row):
    return row[1:].apply(lambda status: 1 if any(x in status for x in ['2_0', '2_NA', '2_2']) else 0 if any(x in status for x in ['0_2', '0_NA', '0_0']) else pd.NA)

def check_subject_presence(row):
    return row[1:].apply(lambda status: 1 if any(x in status for x in ['0_2', 'NA_2', '2_2']) else 0 if any(x in status for x in ['2_0', 'NA_0', '0_0']) else pd.NA)

def main():
    parser = argparse.ArgumentParser(description='simulate sequence without HGT')

    parser.add_argument('-i', '--hgt_info', required=True, help='HGT events information table.')
    parser.add_argument('-s', '--contig_seq', required=True, help='contig sequence from last step. eg. sequences_contig_all_FU.fa')
    parser.add_argument('-o', '--output', default="./", help='output path')
    args = parser.parse_args()

    hgt_clusters = pd.read_csv(args.hgt_info)

    query_info = hgt_clusters[['Sequence_query', 'MAG 1', 'query_congtig_id', 'q_start', 'q_end', 'query_length', 'Length']].copy()
    query_info['HGT_ID'] = query_info['Sequence_query']
    query_info.rename(columns={'MAG1_Group': 'MAG_info', 'MAG 1': 'MAG', 'query_congtig_id': 'contig_id', 'q_start': 'start', 'q_end': 'end', 'query_length': 'length'}, inplace=True)
    query_info['type'] = 'query'

    subj_info = hgt_clusters[['Sequence_subject', 'MAG 2', 'subject_contig_id', 's_start', 's_end', 'subject_length', 'Length']].copy()
    subj_info['HGT_ID'] = subj_info['Sequence_subject']
    subj_info.rename(columns={'MAG2_Group': 'MAG_info', 'MAG 2': 'MAG', 'subject_contig_id': 'contig_id', 's_start': 'start', 's_end': 'end', 'subject_length': 'length'}, inplace=True)
    subj_info['type'] = 'subject'

    elements_info = pd.concat([query_info,subj_info], ignore_index=True).drop_duplicates()
    elements_info_updated = elements_info.drop(['Sequence_subject', 'Sequence_query'], axis=1).drop_duplicates(subset=['HGT_ID'], keep='first')#.set_index('HGT_ID')
    elements_info_updated.to_csv(f"{args.output}/elements_info.csv", index=True)

    #hgt_info = "HGT_events/HGT_Clusters_info_filter_conseve_Final.csv"
    #contigs_path = "sequences_contig_all_FU1.fa"
    #output_dir = "output"
    extract_and_simulate_sequences(elements_info_updated, args.contig_seq, args.output)

if __name__ == "__main__":
    main()

