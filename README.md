# HDMI (HGT Detection from MAGs in Individual)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-available-green.svg)](https://docs.conda.io/)

A comprehensive and user-friendly pipeline for detecting and analyzing horizontal gene transfer (HGT) events in metagenomic data from MAGs (Metagenome-Assembled Genomes).

**Author**: Haoran Peng (penghr21@gmail.com)  
**GitHub**: https://github.com/HaoranPeng21/HDMI  
**Documentation**: This README

## 🚀 Quick Start

```bash
# 1. Clone and install
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI
mamba env create -f environment.yml
conda activate HDMI
pip install -e .

# 2. Run complete pipeline
HDMI index -g genome_folder -m Group_info_test.txt -o output
HDMI detect -i genome_folder -o output -m Group_info_test.txt
HDMI validate -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output -g genome_folder -m Group_info_test.txt
HDMI merge -o output -group Group_info_test.txt
HDMI connect -o output
HDMI profile -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output
HDMI summary -o output
```

## 📋 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start Guide](#quick-start-guide)
- [Detailed Usage](#detailed-usage)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Performance Tips](#performance-tips)
- [FAQ](#faq)
- [Citation](#citation)

## 🔍 Overview

HDMI is a comprehensive pipeline designed to detect horizontal gene transfer (HGT) events in metagenomic data. The pipeline consists of 7 main steps that work together to identify, validate, and analyze HGT events:

### Pipeline Steps

1. **🔍 Index** - Pre-build genome indices for faster processing
2. **🎯 Detect** - Find HGT candidates using BLAST analysis
3. **✅ Validate** - Validate HGT events for individual samples
4. **🔄 Merge** - Merge and filter results from multiple samples
5. **🔗 Connect** - Extract HGT sequences and generate simulated sequences
6. **📊 Profile** - Analyze read coverage for simulated sequences
7. **📋 Summary** - Generate final element table with metagenomic evidence

### Key Features

- ✅ **Easy Installation**: One-command setup with conda/mamba
- ✅ **User-Friendly**: Simple commands with automatic file detection
- ✅ **Batch Processing**: Support for large datasets with parallel processing
- ✅ **Comprehensive Output**: Detailed results with multiple validation levels
- ✅ **Clean File Names**: Final output files without confusing numbers
- ✅ **Auto-Detection**: Automatically finds required files and directories

## 📦 Installation

### Prerequisites

- **Python**: 3.7 or higher
- **Conda/Mamba**: For environment management
- **System Tools**: BLAST+, Bowtie2, Samtools (automatically installed via conda)

### Installation Methods

#### 🎯 Method 1: Mamba Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Create environment with mamba (faster than conda)
mamba env create -f environment.yml

# Activate environment
conda activate HDMI

# Install HDMI package for command-line access
pip install -e .
```

#### 🔧 Method 2: Conda Installation

```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate HDMI

# Install HDMI package
pip install -e .
```

#### ⚡ Method 3: Direct Use (No Installation)

```bash
# Clone the repository
git clone https://github.com/HaoranPeng21/HDMI.git
cd HDMI

# Create environment
conda env create -f environment.yml
conda activate HDMI

# Use directly without installation
python HDMI.py --help
```

### Verification

After installation, verify that HDMI is working:

```bash
# Check if HDMI command is available
HDMI --help

# Should show available commands:
# detect, validate, merge, index, connect, profile, summary
```

## 🚀 Quick Start Guide

### Step-by-Step Example

Here's a complete example using the provided test data:

```bash
# 1. Build genome indices (run once)
HDMI index -g genome_folder -m Group_info_test.txt -o output

# 2. Detect HGT candidates
HDMI detect -i genome_folder -o output -m Group_info_test.txt

# 3. Validate a sample
HDMI validate -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output -g genome_folder -m Group_info_test.txt

# 4. Merge and filter results
HDMI merge -o output -group Group_info_test.txt

# 5. Connect sequences
HDMI connect -o output

# 6. Profile analysis
HDMI profile -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output

# 7. Generate final summary
HDMI summary -o output
```

### Expected Results

After running the complete pipeline, you should see:

- **42 original HGT events** detected
- **18 filtered HGT events** after validation
- **36 final HGT elements** in the element table
- **Clean file names** without confusing numbers

## 📖 Detailed Usage

### HDMI index

Build genome indices for faster processing. This step only needs to be run once per genome set.

```bash
HDMI index -g <genome_folder> -m <group_info> -o <output>
```

**Parameters:**
- `-g, --genome_path`: Path to folder containing genome FASTA files
- `-m, --group_info`: Group information file (see format below)
- `-o, --output`: Output directory (will create index subfolder)

**Example:**
```bash
HDMI index -g genomes/ -m groups.txt -o results/
```

### HDMI detect

Find HGT candidates using BLAST analysis. Supports batch processing for large datasets.

```bash
HDMI detect -i <genome_folder> -o <output> -m <group_info> [-number <batch_num> -total <total_batches>]
```

**Parameters:**
- `-i, --genome_path`: Path to genome folder
- `-o, --output`: Output directory
- `-m, --group_info`: Group information file
- `-number, --task_number`: Batch number (for parallel processing)
- `-total, --total_tasks`: Total number of batches
- `--count-only`: Show genome pair count without running detection

**Examples:**
```bash
# Single batch (default)
HDMI detect -i genomes/ -o results/ -m groups.txt

# Batch processing (2 batches)
HDMI detect -i genomes/ -o results/ -m groups.txt -number 1 -total 2
HDMI detect -i genomes/ -o results/ -m groups.txt -number 2 -total 2
```

### HDMI validate

Validate HGT events for individual samples using read mapping.

```bash
HDMI validate -r1 <read1> -r2 <read2> [--prefix <sample_prefix>] -o <output> -g <genome_folder> -m <group_info> [--threads <threads>]
```

**Parameters:**
- `-r1, --read1`: Read 1 file (supports .fq.gz, .fastq.gz, .fq, .fastq)
- `-r2, --read2`: Read 2 file
- `--prefix`: Sample prefix (auto-extracted from filename if not provided)
- `-o, --output`: Output directory
- `-g, --genome_path`: Path to genome folder
- `-m, --group_info`: Group information file
- `--threads`: Number of threads (default: 8)

**Examples:**
```bash
# With auto-extracted prefix
HDMI validate -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/ -g genomes/ -m groups.txt

# With custom prefix
HDMI validate -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz --prefix my_sample -o results/ -g genomes/ -m groups.txt --threads 16
```

### HDMI merge

Merge and filter results from multiple samples.

```bash
HDMI merge -o <output> -group <group_info> [--threshold <threshold>]
```

**Parameters:**
- `-o, --output`: Output directory
- `-group, --group_info`: Group information file
- `--threshold`: Abundance threshold (default: 1.0)

**Example:**
```bash
HDMI merge -o results/ -group groups.txt --threshold 0.5
```

### HDMI connect

Extract HGT sequences and generate simulated sequences for coverage analysis.

```bash
HDMI connect -o <output>
```

**Parameters:**
- `-o, --output`: Output directory

**Example:**
```bash
HDMI connect -o results/
```

### HDMI profile

Analyze read coverage for simulated sequences.

```bash
HDMI profile -r1 <read1> -r2 <read2> [--prefix <sample_prefix>] -o <output> [--threads <threads>]
```

**Parameters:**
- `-r1, --read1`: Read 1 file
- `-r2, --read2`: Read 2 file
- `--prefix`: Sample prefix (auto-extracted if not provided)
- `-o, --output`: Output directory
- `--threads`: Number of threads (default: 1)

**Example:**
```bash
HDMI profile -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/ --threads 8
```

### HDMI summary

Generate final element table with metagenomic evidence. This is the final step that produces the main results.

```bash
HDMI summary -o <output>
```

**Parameters:**
- `-o, --output`: Output directory (default: output)

**Example:**
```bash
HDMI summary -o results/
```

**Note**: This command automatically finds all required files and generates the final `element_table.csv`.

## 📁 Input File Formats

### Group Information File

A tab-separated file that defines genome groups and representatives.

**Format:**
```
<genome_name>	<group_number>
```

**Example (Group_info_test.txt):**
```
bin.1.fa	1
bin.2.fa	2
bin3.fa	3
bin4.fa	4
```

**Description:**
- **genome_name**: Name of the genome FASTA file (without path)
- **group_number**: Integer group ID (1, 2, 3, etc.)

### Genome Files

- **Format**: FASTA (.fa, .fasta)
- **Location**: Specified genome folder
- **Naming**: Should match names in group info file

**Example genome folder structure:**
```
genome_folder/
├── bin.1.fa
├── bin.2.fa
├── bin3.fa
└── bin4.fa
```

### Read Files

- **Format**: FASTQ (.fq.gz, .fastq.gz, .fq, .fastq)
- **Type**: Paired-end reads
- **Naming**: HDMI auto-extracts sample prefix from filename

**Example read files:**
```
data/
├── sample1_R1.fq.gz
├── sample1_R2.fq.gz
├── sample2_R1.fq.gz
└── sample2_R2.fq.gz
```

## 📊 Output Files

### Main Output Files

```
output/
├── element_table.csv           # 🎯 MAIN RESULT: Final HGT elements with evidence
├── HGT_events.csv              # Filtered HGT events after validation
├── HGT_events_raw.csv          # Raw HGT events from detection
├── sequences_contig_combined.fa # Combined contig sequences
├── sequences_contig_q.fa       # Query contig sequences
├── sequences_contig_s.fa       # Subject contig sequences
├── sequences_matched_seq_q.fa  # Matched query sequences
├── sequences_matched_seq_s.fa  # Matched subject sequences
└── index/                      # Genome indices (for faster processing)
```

### Intermediate Files

```
output/intermediate/
├── 01_detection/               # HGT detection results
├── 02_validation/              # Sample validation results
├── 03_final/                   # Merge and filter results
├── 04_connect/                 # Sequence connection results
└── 05_profile/                 # Profile analysis results
```

### Output File Descriptions

#### element_table.csv (Main Result)
The final output file containing HGT elements with metagenomic evidence.

**Format:**
```csv
HGT_ID,sample1,sample2,...
NODE_71_length_18390_cov_19.696046_13958_14594,1.0,0.0,...
```

**Columns:**
- **HGT_ID**: Unique HGT element identifier
- **sample1, sample2, ...**: Evidence scores for each sample (1.0 = evidence present, 0.0 = no evidence)

#### HGT_events.csv
Filtered HGT events after validation and quality control.

**Format:**
```csv
Sequence_query,Sequence_subject,query_congtig_id,subject_contig_id,MAG 1,MAG 2,q_start,q_end,s_start,s_end,query_length,subject_length,Length,Identity,End Match,Full Match,HGT_ID,MAG1_Group,MAG2_Group
```

## 📈 Examples

### Example 1: Single Sample Analysis

```bash
# Complete pipeline for one sample
HDMI index -g genomes/ -m groups.txt -o results/
HDMI detect -i genomes/ -o results/ -m groups.txt
HDMI validate -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/ -g genomes/ -m groups.txt
HDMI merge -o results/ -group groups.txt
HDMI connect -o results/
HDMI profile -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/
HDMI summary -o results/
```

### Example 2: Multiple Samples

```bash
# Build indices once
HDMI index -g genomes/ -m groups.txt -o results/

# Detect HGT candidates once
HDMI detect -i genomes/ -o results/ -m groups.txt

# Validate each sample
HDMI validate -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/ -g genomes/ -m groups.txt
HDMI validate -r1 sample2_R1.fq.gz -r2 sample2_R2.fq.gz -o results/ -g genomes/ -m groups.txt
HDMI validate -r1 sample3_R1.fq.gz -r2 sample3_R2.fq.gz -o results/ -g genomes/ -m groups.txt

# Merge all samples
HDMI merge -o results/ -group groups.txt

# Connect sequences
HDMI connect -o results/

# Profile each sample
HDMI profile -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o results/
HDMI profile -r1 sample2_R1.fq.gz -r2 sample2_R2.fq.gz -o results/
HDMI profile -r1 sample3_R1.fq.gz -r2 sample3_R2.fq.gz -o results/

# Generate final summary
HDMI summary -o results/
```

### Example 3: Batch Processing for Large Datasets

```bash
# For large genome sets, use batch processing
HDMI detect -i genomes/ -o results/ -m groups.txt -number 1 -total 4
HDMI detect -i genomes/ -o results/ -m groups.txt -number 2 -total 4
HDMI detect -i genomes/ -o results/ -m groups.txt -number 3 -total 4
HDMI detect -i genomes/ -o results/ -m groups.txt -number 4 -total 4
```

## 🔧 Troubleshooting

### Common Issues and Solutions

#### 1. "HDMI: command not found"
**Problem**: HDMI command is not available after installation.

**Solution**:
```bash
# Make sure you're in the HDMI environment
conda activate HDMI

# Install HDMI package
pip install -e .

# Verify installation
HDMI --help
```

#### 2. "BLAST not found" or "Bowtie2 not found"
**Problem**: Required tools are not installed.

**Solution**:
```bash
# Reinstall environment
conda deactivate
conda env remove -n HDMI
mamba env create -f environment.yml
conda activate HDMI
pip install -e .
```

#### 3. "Profile directory not found"
**Problem**: Summary command can't find profile results.

**Solution**:
```bash
# Make sure you've run all previous steps
HDMI profile -r1 reads_R1.fq.gz -r2 reads_R2.fq.gz -o output
HDMI summary -o output
```

#### 4. Memory Issues
**Problem**: Out of memory errors with large datasets.

**Solution**:
```bash
# Use batch processing
HDMI detect -i genomes/ -o output/ -m groups.txt -number 1 -total 4
HDMI detect -i genomes/ -o output/ -m groups.txt -number 2 -total 4
# ... continue for all batches
```

#### 5. File Permission Errors
**Problem**: Cannot write to output directory.

**Solution**:
```bash
# Check permissions
ls -la output/

# Fix permissions if needed
chmod 755 output/
```

### Error Messages and Solutions

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `HDMI: command not found` | Package not installed | Run `pip install -e .` |
| `BLAST not found` | Missing dependencies | Reinstall environment |
| `Profile directory not found` | Missing profile step | Run `HDMI profile` first |
| `Memory error` | Dataset too large | Use batch processing |
| `File not found` | Wrong file path | Check file locations |

## ⚡ Performance Tips

### 1. Use Mamba Instead of Conda
```bash
# Install mamba first
conda install mamba -n base -c conda-forge

# Then use mamba for faster environment creation
mamba env create -f environment.yml
```

### 2. Batch Processing for Large Datasets
```bash
# For large genome sets, split into batches
HDMI detect -i genomes/ -o output/ -m groups.txt -number 1 -total 8
HDMI detect -i genomes/ -o output/ -m groups.txt -number 2 -total 8
# ... continue for all 8 batches
```

### 3. Optimize Thread Usage
```bash
# Use appropriate thread counts
HDMI validate -r1 reads_R1.fq.gz -r2 reads_R2.fq.gz -o output/ -g genomes/ -m groups.txt --threads 16
HDMI profile -r1 reads_R1.fq.gz -r2 reads_R2.fq.gz -o output/ --threads 8
```

### 4. Monitor Disk Space
```bash
# Check available space
df -h

# Clean up intermediate files if needed
rm -rf output/intermediate/01_detection/*.bam
```

### 5. Use SSD Storage
For better performance, run HDMI on SSD storage when possible.

## ❓ FAQ

### Q: How long does the pipeline take to run?
**A**: Runtime depends on dataset size:
- Small dataset (4 genomes, 1 sample): ~30 minutes
- Medium dataset (20 genomes, 5 samples): ~2-3 hours
- Large dataset (100 genomes, 20 samples): ~8-12 hours

### Q: Can I run steps in parallel?
**A**: Yes! You can run validation and profile steps for different samples in parallel:
```bash
# Run multiple samples simultaneously
HDMI validate -r1 sample1_R1.fq.gz -r2 sample1_R2.fq.gz -o output/ -g genomes/ -m groups.txt &
HDMI validate -r1 sample2_R1.fq.gz -r2 sample2_R2.fq.gz -o output/ -g genomes/ -m groups.txt &
wait
```

### Q: What if I have single-end reads?
**A**: HDMI currently requires paired-end reads. Consider using tools to convert single-end to paired-end or contact the developers for single-end support.

### Q: How do I interpret the element_table.csv results?
**A**: The element_table.csv contains HGT elements with evidence scores:
- **1.0**: Strong evidence for HGT
- **0.0**: No evidence for HGT
- Each column represents a sample

### Q: Can I use my own genome annotations?
**A**: HDMI works with genome sequences (FASTA files). Annotations are not required but can be used for downstream analysis.

### Q: What's the difference between HGT_events.csv and element_table.csv?
**A**: 
- **HGT_events.csv**: Raw HGT events with detailed BLAST information
- **element_table.csv**: Final HGT elements with metagenomic evidence scores

## 📚 Citation

If you use HDMI in your research, please cite:

```bibtex
@software{hdmi2024,
  title={HDMI: HGT Detection from MAGs in Individual},
  author={Peng, Haoran},
  year={2024},
  url={https://github.com/HaoranPeng21/HDMI}
}
```

## 🤝 Contributing

We welcome contributions! Please feel free to:
- Report bugs via GitHub issues
- Suggest new features
- Submit pull requests
- Improve documentation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

- **Author**: Haoran Peng (penghr21@gmail.com)
- **GitHub**: [https://github.com/HaoranPeng21/HDMI](https://github.com/HaoranPeng21/HDMI)
- **Issues**: [https://github.com/HaoranPeng21/HDMI/issues](https://github.com/HaoranPeng21/HDMI/issues)
- **Documentation**: This README

---

**Happy HGT Detection! 🧬**

