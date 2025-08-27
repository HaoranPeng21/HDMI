# HDMI (HGT Detection from MAGs in Individual) Pipeline

A comprehensive pipeline for detecting and analyzing horizontal gene transfer (HGT) events in metagenomic data.

Author: Haoran Peng (penghr21@gmail.com)
GitHub: https://github.com/HaoranPeng21/HDMI

## Overview

The HDMI pipeline consists of 7 main steps:

1. **Index** - Build genome indices for faster processing
2. **Detect** - Find HGT candidates using BLAST analysis
3. **Validate** - Validate HGT events for individual samples
4. **Merge** - Merge and filter results from multiple samples
5. **Connect** - Extract HGT sequences and generate simulated sequences
6. **Profile** - Analyze read coverage for simulated sequences
7. **Summary** - Generate final element table with metagenomic evidence

## Installation

### Prerequisites

- Python 3.7+
- BLAST+
- Bowtie2
- Samtools
- Conda (recommended)

### Installation Methods

#### Method 1: Quick Install (Recommended)
```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Run the installation script
./install.sh
```

#### Method 2: Manual Installation
```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate HDMI

# Install HDMI
pip install -e .
```

#### Method 3: Direct Use
```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Make HDMI script executable
chmod +x bin/HDMI

# Use directly
./bin/HDMI --help
```

## Quick Start

### Complete Pipeline Example

```bash
# Step 1: Build genome indices
HDMI index -g genome_folder -m Group_info_test.txt -o output

# Step 2: Detect HGT candidates (batch processing)
HDMI detect -i genome_folder -o output -m Group_info_test.txt -number 1 -total 2
HDMI detect -i genome_folder -o output -m Group_info_test.txt -number 2 -total 2

# Step 3: Validate samples
HDMI validate -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz --prefix sample1 -o output -g genome_folder -m Group_info_test.txt --threads 5
HDMI validate -r1 data/sample2_R1.fq.gz -r2 data/sample2_R2.fq.gz --prefix sample2 -o output -g genome_folder -m Group_info_test.txt --threads 5

# Step 4: Merge and filter results
HDMI merge -o output -group Group_info_test.txt

# Step 5: Connect sequences
HDMI connect -o output

# Step 6: Profile analysis
HDMI profile -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz --prefix sample1 -o output --threads 5
HDMI profile -r1 data/sample2_R1.fq.gz -r2 data/sample2_R2.fq.gz --prefix sample2 -o output --threads 5

# Step 7: Generate final summary
HDMI summary -o output/element_table.csv
```

## Command Reference

### HDMI index

Build genome indices for faster processing.

```bash
HDMI index -g <genome_folder> -m <group_info> -o <output>
```

**Parameters:**
- `-g, --genome_path`: Path to genome folder
- `-m, --group_info`: Group information file
- `-o, --output`: Output directory (auto-creates index subfolder)

### HDMI detect

Find HGT candidates using BLAST analysis.

```bash
HDMI detect -i <genome_folder> -o <output> -m <group_info> [-number <batch_num> -total <total_batches>]
```

**Parameters:**
- `-i, --genome_path`: Path to genome folder
- `-o, --output`: Output directory (auto-creates intermediate/01_detection)
- `-m, --group_info`: Group information file
- `-number, --task_number`: Batch number (for batch processing)
- `-total, --total_tasks`: Total number of batches
- `--count-only`: Show genome pair count and performance estimates only

### HDMI validate

Validate HGT events for a single sample.

```bash
HDMI validate -r1 <read1> -r2 <read2> [--prefix <sample_prefix>] -o <output> -g <genome_folder> -m <group_info> [--threads <threads>]
```

**Parameters:**
- `-r1, --read1`: Read 1 file (supports .fq.gz, .fastq.gz, .fq, .fastq)
- `-r2, --read2`: Read 2 file
- `--prefix`: Sample prefix (auto-extracted from filename if not provided)
- `-o, --output`: Output directory (auto-creates intermediate/02_validation/sample_name)
- `-g, --genome_path`: Path to genome folder
- `-m, --group_info`: Group information file
- `--threads`: Number of threads (default: auto-detected)

### HDMI merge

Merge and filter results from multiple samples.

```bash
HDMI merge -o <output> -group <group_info> [--threshold <threshold>]
```

**Parameters:**
- `-o, --output`: Output directory (auto-creates intermediate/03_final)
- `-group, --group_info`: Group information file
- `--threshold`: Abundance threshold (default: 1.0)

### HDMI connect

Extract HGT sequences and generate simulated sequences.

```bash
HDMI connect -o <output>
```

**Parameters:**
- `-o, --output`: Output directory (auto-creates intermediate/04_connect)

### HDMI profile

Analyze read coverage for simulated sequences.

```bash
HDMI profile -r1 <read1> -r2 <read2> [--prefix <sample_prefix>] -o <output> [--threads <threads>]
```

**Parameters:**
- `-r1, --read1`: Read 1 file
- `-r2, --read2`: Read 2 file
- `--prefix`: Sample prefix (auto-extracted from filename if not provided)
- `-o, --output`: Output directory (auto-creates intermediate/05_profile/sample_name)
- `--threads`: Number of threads (default: auto-detected)

### HDMI summary

Generate final element table with metagenomic evidence.

```bash
HDMI summary -o <output_file>
```

**Parameters:**
- `-o, --output`: Output file path (e.g., output/element_table.csv)

## Output Structure

```
output/
├── HGT_events.csv              # Final filtered HGT events
├── HGT_events_raw.csv          # Raw HGT events from detection
├── element_table.csv           # Final element table with metagenomic evidence
├── sequences_contig_combined.fa # Combined contig sequences
├── sequences_contig_q.fa       # Query contig sequences
├── sequences_contig_s.fa       # Subject contig sequences
├── sequences_matched_seq_q.fa  # Matched query sequences
├── sequences_matched_seq_s.fa  # Matched subject sequences
├── index/                      # Genome indices
├── intermediate/               # Intermediate files
│   ├── 01_detection/           # Batch processing files
│   ├── 02_validation/          # Sample validation results
│   ├── 03_final/               # Merge and filter results
│   ├── 04_connect/             # Sequence connection results
│   └── 05_profile/             # Profile analysis results
└── log/                       # Log files
```

## File Formats

### Group Information File

Tab-separated file with columns:
- Genome name
- Group number
- Representative flag (1 for representative, 0 for non-representative)

Example:
```
bin.1	1	1
bin.2	2	1
bin3	3	1
bin4	4	1
```

### HGT Events File

CSV file with columns:
- HGT_ID: Unique identifier
- Query_contig: Query contig name
- Subject_contig: Subject contig name
- Query_genome: Query genome name
- Subject_genome: Subject genome name
- Query_group: Query genome group
- Subject_group: Subject genome group
- Identity: BLAST identity
- Coverage: BLAST coverage
- E_value: BLAST E-value

## Troubleshooting

### Common Issues

1. **BLAST not found**: Ensure BLAST+ is installed and in PATH
2. **Bowtie2 not found**: Ensure Bowtie2 is installed and in PATH
3. **Memory issues**: Use batch processing for large datasets
4. **File not found errors**: Check file paths and permissions

### Performance Tips

1. Use batch processing for large genome sets
2. Pre-build indices for faster validation
3. Use appropriate thread counts for your system
4. Monitor disk space for intermediate files

## Version

**v1.0**: Initial release with complete HGT detection pipeline

## Citation

If you use HDMI in your research, please cite:
```
HDMI: HGT Detection from MAGs in Individual
Haoran Peng, 2024
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

- GitHub: [https://github.com/HaoranPeng21/HDMI](https://github.com/HaoranPeng21/HDMI)
- Issues: [https://github.com/HaoranPeng21/HDMI/issues](https://github.com/HaoranPeng21/HDMI/issues)

