import importlib.util
import subprocess

def check_and_install_module(module_name):
    if importlib.util.find_spec(module_name) is None:
        print(f"{module_name} is not installed. Installing...")
        subprocess.check_call(["pip", "install", module_name])
    else:
        print(f"{module_name} is already installed.")

modules = ["os", "csv", "concurrent.futures", "argparse","biopython"]

for module in modules:
    check_and_install_module(module)

import os
import csv
import concurrent.futures
import argparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# read table get MAGs info
def read_mag_info(group_info):
    mag_info = {}
    with open(f"{group_info}", "r") as txtfile:
        reader = csv.reader(txtfile)
        for row in reader:
            mag_info[row[0]] = row[1]
    return mag_info

# read fasta file and store in dic
def read_fasta_files(mag_info, genome_path):
    fasta_dict = {}
    for mag in mag_info.keys():
        with open(os.path.join(genome_path, mag), "r") as fasta_file:
            # Store records in a dictionary using record.id as the key.
            fasta_dict[mag] = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    return fasta_dict

def check_match_direction(query_start,query_end,subject_start,subject_end):
    query_direction = query_end - query_start
    subject_direction = subject_end - subject_start

    same_match_direction = True
    if ((query_direction > 0) and (subject_direction < 0)) or ((query_direction < 0) and (subject_direction > 0)):
        same_match_direction = False

    return same_match_direction

def check_end_match(query_start,query_end,query_len,subject_start,subject_end,subject_len,best_hit_end_gap_len):
    
    best_hit_same_direction = check_match_direction(query_start,query_end,subject_start,subject_end)
    
    match_category = False

    if (best_hit_same_direction is True) and (query_len - query_end <= best_hit_end_gap_len) and (subject_start <= best_hit_end_gap_len):
        match_category = True

    # situation 2
    elif (best_hit_same_direction is True) and (query_start <= best_hit_end_gap_len) and (subject_len - subject_end <= best_hit_end_gap_len):
        match_category = True

    # situation 3
    elif (best_hit_same_direction is False) and (query_len - query_end <= best_hit_end_gap_len) and (subject_len - subject_start <= best_hit_end_gap_len):
        match_category = True

    # situation 4
    elif (best_hit_same_direction is False) and (query_start <= best_hit_end_gap_len) and (subject_end <= best_hit_end_gap_len):
        match_category = True
    
    return match_category

# Blast
def blast_comparison(mag_pair, genome_path, fasta_dict, pair_index, total_pairs,task_number,log_path):
    mag1, group1, mag2, group2 = mag_pair
    output_table = set()
    output_fasta_q = {}
    output_match_fasta_q = {}
    output_fasta_s = {}
    output_match_fasta_s = {}
    log_filename = f"{log_path}/log/output.log_{task_number}"
    with open(log_filename, "a") as log_file:
        log_file.write(f"Task {task_number} - Processing pair {pair_index + 1}/{total_pairs}: {mag1} vs {mag2}\n")

    blastn_cline = NcbiblastnCommandline(
        query=os.path.join(genome_path, mag1),
        subject=os.path.join(genome_path, mag2),
        outfmt=6)

    stdout, stderr = blastn_cline()
    blast_output = StringIO(stdout)

    reader = csv.reader(blast_output, delimiter="\t")
    for row in reader:
        query_seq_id, subject_seq_id, identity, length, q_start, q_end, s_start, s_end = row[0], row[1], float(row[2]), int(row[3]), int(row[6]), int(row[7]), int(row[8]), int(row[9])

        if identity >= 99.0 and length >= 500:
            # Access records directly using IDs
            query_record = fasta_dict[mag1].get(query_seq_id)
            subject_record = fasta_dict[mag2].get(subject_seq_id)

            # If either record is None, skip this result
            if query_record is None or subject_record is None:
                continue

            query_length = len(query_record)
            subject_length = len(subject_record)

            best_hit_end_gap_len = 100
            end_match = check_end_match(q_start, q_end, query_length, s_start, s_end, subject_length, best_hit_end_gap_len)

            full_match = (length >= query_length*0.90) or (length >= subject_length*0.90)


            output_table.add((f"{query_seq_id}_{q_start}_{q_end}", f"{subject_seq_id}_{s_start}_{s_end}" ,query_seq_id,subject_seq_id, mag1, mag2, q_start, q_end, s_start, s_end, query_length, subject_length, length, identity, end_match, full_match))
            
            output_fasta_q[query_seq_id] = query_record

            output_fasta_s[subject_seq_id] = subject_record

            query_direction = q_end - q_start
            subject_direction = s_end - s_start

            if query_direction > 0:
                matched_seq_q = query_record.seq[q_start - 1:q_end]
            else:
                matched_seq_q = query_record.seq[q_end - 1:q_start]

            matched_record_q = SeqRecord(Seq(str(matched_seq_q)), id=f"{query_seq_id}_{q_start}_{q_end}", description="Matched sequence")
            
            if subject_direction >0:
                matched_seq_s = subject_record.seq[s_start - 1:s_end]
            else:
                matched_seq_s = subject_record.seq[s_end - 1:s_start]

            matched_record_s = SeqRecord(Seq(str(matched_seq_s)), id=f"{subject_seq_id}_{s_start}_{s_end}", description="Matched sequence")

            output_match_fasta_q[query_seq_id + "_"+ str(q_start) + "_" + str(q_end)] = matched_record_q
            output_match_fasta_s[subject_seq_id + "_" + str(s_start) + "_" + str(s_end)] = matched_record_s

    return output_table, output_fasta_q,output_fasta_s, output_match_fasta_q,output_match_fasta_s


def main():
    parser = argparse.ArgumentParser(description='HGT detection pipeline step1')

    parser.add_argument('-i', '--path', required=True, help='Path to genomes FASTA.')
    parser.add_argument('-o', '--output', default="output", help='the output folder')
    parser.add_argument('-m', '--table', required=True, help='genomes groups information')
    parser.add_argument("-number", "--task_number", type=int, default=1, help="Task number (1-indexed)")
    parser.add_argument("-total", "--total_tasks", type=int, default=1, help="Total number of tasks")
    args = parser.parse_args()
    task_number = args.task_number
    total_tasks = args.total_tasks
    mag_info = read_mag_info(args.table)
    fasta_dict = read_fasta_files(mag_info, args.path)

    # mkdir dir
    os.makedirs(f"{args.output}", exist_ok=True)
    os.makedirs(f"{args.output}/log", exist_ok=True)
    mag_list = list(mag_info.items())
    mag_pairs = [(mag1, group1, mag2, group2) for (mag1, group1) in mag_list for (mag2, group2) in mag_list if mag1 < mag2 and group1 != group2]
    
    # split mag_pairs list
    chunk_size = len(mag_pairs) // total_tasks
    start_index = chunk_size * (task_number - 1)
    end_index = start_index + chunk_size if task_number < total_tasks else None
    selected_mag_pairs = mag_pairs[start_index:end_index]

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(lambda pair_index_and_pair: blast_comparison(pair_index_and_pair[1], args.path, fasta_dict, pair_index_and_pair[0], len(selected_mag_pairs), task_number, args.output), enumerate(selected_mag_pairs)))

    final_output_table = set()
    final_output_fasta_q = {}
    final_matched_sequences_q = {}
    final_output_fasta_s = {}
    final_matched_sequences_s = {}

    for output_table, output_fasta_q,output_fasta_s,output_match_fasta_q,output_match_fasta_s in results:
        final_output_table |= output_table
        final_output_fasta_q.update(output_fasta_q)
        final_matched_sequences_q.update(output_match_fasta_q)
        final_output_fasta_s.update(output_fasta_s)
        final_matched_sequences_s.update(output_match_fasta_s)
    # write result
    output_filename = f"{args.output}/New_HGT_table{task_number}.csv"
    output_fasta_q_filename = f"{args.output}/sequences_contig_q{task_number}.fa"
    output_matched_q_filename = f"{args.output}/sequences_matched_seq_q{task_number}.fa"
    output_fasta_s_filename = f"{args.output}/sequences_contig_s{task_number}.fa"
    output_matched_s_filename = f"{args.output}/sequences_matched_seq_s{task_number}.fa"
    

    with open(output_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sequence_query","Sequence_subject", "query_congtig_id", "subject_contig_id","MAG 1", "MAG 2", "q_start", "q_end", "s_start", "s_end", "query_length", "subject_length", "Length", "Identity", "End Match", "Full Match"])
        writer.writerows(final_output_table)

    SeqIO.write(final_output_fasta_q.values(), output_fasta_q_filename, "fasta")
    SeqIO.write(final_matched_sequences_q.values(), output_matched_q_filename, "fasta")
    SeqIO.write(final_output_fasta_s.values(), output_fasta_s_filename, "fasta")
    SeqIO.write(final_matched_sequences_s.values(), output_matched_s_filename, "fasta")
    print(f"Task {task_number} Done!")


if __name__ == "__main__":
    main()
