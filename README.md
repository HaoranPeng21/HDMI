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

# 2. Run complete pipeline (4 optimized steps)
HDMI detect -i genome_folder -o output -m Group_info_test.txt -t 10
HDMI index -g genome_folder -m Group_info_test.txt -o output -t 10
HDMI profile -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output -g genome_folder -m Group_info_test.txt -t 10
HDMI summary -o output -group Group_info_test.txt
```

## 📋 Table of Contents

- [Overview](#-overview)
- [Installation](#-installation)
- [Detailed Usage](#-detailed-usage)
- [Input File Formats](#-input-file-formats)
- [Output Files](#-output-files)
- [Citation](#-citation)
- [Contributing](#-contributing)
- [License](#-license)
- [Contact](#-contact)

## 🔍 Overview

HDMI is a comprehensive pipeline designed to detect horizontal gene transfer (HGT) events in metagenomic data. The pipeline consists of 4 optimized steps that work together to identify, validate, and analyze HGT events:

1. **HDMI detect**: Identifies HGT candidates between genomes using BLAST-based similarity search
2. **HDMI index**: Builds Bowtie2 indices for contigs contain HGT
3. **HDMI profile**: Profiles and validates events through read spanning analysis
4. **HDMI summary**: Merges results across samples, filters validated events, and generates final element table

The pipeline is optimized for high-throughput analysis and provides comprehensive validation of HGT events through multiple criteria including read spanning, coverage fraction, and species abundance.

## 🛠 Installation

### Prerequisites

- Python 3.7 or higher
- Conda or Mamba package manager
- At least 16GB RAM (32GB recommended for large datasets)
- Sufficient disk space for intermediate files (typically 2-3x the size of input data)

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/HaoranPeng21/HDMI.git
   cd HDMI
   ```

2. **Create and activate the conda environment**:
   ```bash
   # Option 1: Using conda (slower but more compatible)
   conda env create -f environment.yml
   
   # Option 2: Using mamba (faster, recommended)
   # mamba env create -f environment.yml
   
   conda activate HDMI
   ```

3. **Install HDMI (development install)**:
   ```bash
   pip install -e .
   ```

### Install via Conda (pending...)

```bash
conda install -c conda-forge -c bioconda hdmi
# or
mamba install -c conda-forge -c bioconda hdmi
```

4. **Verify installation**:
   ```bash
   HDMI --help
   ```
5. **Test Dataset installation**:
   ```
   https://zenodo.org/records/17515026
   ```

## 📖 Detailed Usage

### Step 1: HDMI detect - HGT Candidate Detection

**Purpose**: Identifies potential HGT events between MAGs using BLAST-based similarity search. 

**Note**: Highly recommend using high-quality MAGs (e.g., from co-assembly [fairy](https://github.com/bluenote-1577/fairy), passing [GUNC](https://grp-bork.embl-community.io/gunc/) quality control)

**Command**:
```bash
HDMI detect -i <genome_folder> -o <output_dir> -m <group_info> [options]
```

**Required Parameters**:
- `-i, --input`: Directory containing MAG files (FASTA format)
- `-o, --output`: Output directory for results
- `-m, --mapping`: Group information file (MAG to species mapping)

**Optional Parameters**:
- `-t, --threads`: Number of threads (default: 10)
- `-n, --number`: Batch number for parallel processing (default: 1)
- `--total`: Total number of batches (default: 1)
- `--count-only`: Only count genome pairs without running BLAST (useful for estimating runtime)

**Example**:
```bash
# Full detection
HDMI detect -i genome_folder -o output -m Group_info_test.txt -t 10

# Count-only mode (estimate pairs and runtime - highly recommended before full run)
HDMI detect -i genome_folder -o output -m Group_info_test.txt --count-only
```

**Performance Data** (test_real dataset):
- **Runtime**: 2.3 hours per batch (112,000 genome pairs - use `--count-only` to estimate)
- **Memory**: ~16GB peak usage
- **Output**: HGT_events_raw.csv with ~32,000 HGT candidates
- **Batch processing**: 10 batches total for 1.12 million genome pairs

**Output Files**:
- `HGT_events_raw.csv`: Raw HGT candidates with similarity scores
- `sequences_contig_q.fa`: Query contig sequences
- `sequences_contig_s.fa`: Subject contig sequences
- `sequences_matched_seq_q.fa`: Query HGT region sequences
- `sequences_matched_seq_s.fa`: Subject HGT region sequences
### Step 2: HDMI index - Index Building and Sequence Extraction

**Purpose**: Builds Bowtie2 indices for HGT regions and extracts HGT sequences.

**Command**:
```bash
HDMI index -g <genome_folder> -m <group_info> -o <output_dir> [options]
```

**Required Parameters**:
- `-g, --genome_folder`: Directory containing MAG files
- `-m, --mapping`: Group information file
- `-o, --output`: Output directory

**Optional Parameters**:
- `-t, --threads`: Number of threads (default: 10)

**Example**:
```bash
HDMI index -g genome_folder -m Group_info_test.txt -o output -t 10
```

**Performance Data** (test_real dataset):
- **Runtime**: ~50 minutes
- **Memory**: ~8GB peak usage
- **Processing**: ~32,000 HGT events

**Output Files**:
- `elements_info_raw.csv`: Raw element information with HGT_ID and Element_Type
- Bowtie2 indices for contigs

### Step 3: HDMI profile - Read Mapping and Validation

**Purpose**: Maps reads to HGT regions and validates events through read spanning analysis.

**Command**:
```bash
HDMI profile -r1 <R1_fastq> -r2 <R2_fastq> -o <output_dir> -g <genome_folder> -m <group_info> [options]
```

**Required Parameters**:
- `-r1, --read1`: Forward reads (FASTQ format)
- `-r2, --read2`: Reverse reads (FASTQ format)
- `-o, --output`: Output directory
- `-g, --genome_folder`: Directory containing MAG files
- `-m, --mapping`: Group information file

**Optional Parameters**:
- `-t, --threads`: Number of threads (default: 10)
- `--sth`: Span threshold for read validation (default: 2)

**Example**:
```bash
HDMI profile -r1 data/sample1_R1.fq.gz -r2 data/sample1_R2.fq.gz -o output -g genome_folder -m Group_info_test.txt -t 10
```

**Performance Data** (test_real dataset):
- **Runtime**: ~32 minutes per sample (average)
- **Memory**: ~12GB peak usage
- **Input file size**: 1.5GB per sample (FASTQ.gz)
- **Output file size**: 2.6MB per sample
- **Processing**: 3 samples with ~32,000 HGT events each


### Step 4: HDMI summary - Result Merging and Final Analysis

**Purpose**: Merges results across samples, filters validated events, and generates final element table (integrates merge functionality).

**Command**:
```bash
HDMI summary -o <output_dir> -group <group_info> [options]
```

**Required Parameters**:
- `-o, --output`: Output directory
- `-group, --group_info`: Group information file

**Optional Parameters**:
- `--threshold`: Abundance threshold for species presence (default: 1.0)

**Example**:
```bash
HDMI summary -o output -group Group_info_test.txt
```


**Output Files**:
- `HGT_events.csv`: Filtered and validated HGT events
- `elements_info.csv`: Detailed information for validated HGT elements
- `element_table.csv`: Final element presence/absence table

## 📁 Input File Formats

### Genome Folder Structure
```
genome_folder/
├── MAG_001.fasta
├── MAG_002.fasta
├── MAG_003.fasta
└── ...
```

### Group Information File
Tab-separated file with MAG to species mapping:
```
MAG_ID	Species_ID	Group_ID
MAG_001	Species_A	Group_1
MAG_002	Species_B	Group_1
MAG_003	Species_C	Group_2
```

### Read Files
Paired-end FASTQ files:
- `sample_R1.fq.gz`: Forward reads
- `sample_R2.fq.gz`: Reverse reads

## 📊 Output Files

### Main Output Files (in `output/` directory)

1. **`HGT_events.csv`**: Filtered and validated HGT events
   - Includes validation status for each sample

2. **`elements_info.csv`**: Detailed information for HGT elements
   - Includes HGT_ID, Element_Type, and genomic coordinates

3. **`element_table.csv`**: Final element presence/absence table
   - Shows presence (1.0), absence (0.0), or insufficient data (NA) for each sample

### Intermediate Files (in `output/intermediate/` directory)

- **`01_detection/`**: BLAST results and HGT candidates
- **`02_validation/`**: Per-sample validation results
- **`03_final/`**: Merged and filtered results

## 📚 Citation

If you use HDMI in your research, please cite:

```
@article{fu2024longitudinal,
  title={Longitudinal Gut Microbiota Tracking Reveals the Persistent Spread of Mobile Genes and HGT-Driven Community Stabilization},
  author={Fu, Jingyuan and Peng, Haoran and Andreu-S{\'a}nchez, Sergio and Ruiz-Moreno, Angel and others},
  journal={Research Square},
  year={2024},
  doi={10.21203/rs.3.rs-6509357/v1}
}
```

**Preprint**: https://doi.org/10.21203/rs.3.rs-6509357/v1

## 🤝 Contributing

We welcome contributions to HDMI! Please feel free to:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

For major changes, please open an issue first to discuss what you would like to change.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

- **Author**: Haoran Peng
- **Email**: penghr21@gmail.com
- **GitHub**: https://github.com/HaoranPeng21/HDMI

For questions, bug reports, or feature requests, please open an issue on GitHub or contact the author directly.
