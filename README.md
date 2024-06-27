# HDMI ("HGT Detection from MAGs in Individuals") workflow

## Overview

HDMI workflow is a workflow to detect recent HGT events based on metagenome-assembled-genomes (MAGs). HGTfinder captured individual-specific, recent HGT (0-10000 years) with their bacterial host genomes. 

If you have any questions about [HGTfinder](https://github.com/HaoranPeng21/HGTfinder), feel free to contact me (h.peng@umcg.nl)

## Software requirement

* [Samtool](https://www.htslib.org/)
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [BLAST](https://doi.org/10.1186/1471-2105-10-421)

## Quick Start

### Prepare software

```
conda create -n HGTfinder python=3.7 blast bowtie2 samtools
```



## Step1: Detect HGT from MAGs

```
usage: HGTdetect.py [-h] -i PATH [-o OUTPUT] -m TABLE [-number TASK_NUMBER]
                    [-total TOTAL_TASKS]

HGT detection pipeline step1

optional arguments:
  -h, --help            show this help message and exit
  -i PATH, --path PATH  Path to genomes FASTA.
  -o OUTPUT, --output OUTPUT
                        the output folder
  -m TABLE, --table TABLE
                        genomes groups information
  -number TASK_NUMBER, --task_number TASK_NUMBER
                        Task number (1-indexed)
  -total TOTAL_TASKS, --total_tasks TOTAL_TASKS
                        Total number of tasks
```

#### Input



#### Output



#### Example:

```
python HGTdetect.py -i genome_folder -m Group_info_test.txt
```



## Step 2: HGT profiling in individuals

```
usage: HGTfinder.py [-h] -r1 READ1 -r2 READ2 [-i SAMPLE_ID] [-o OUTPUT]
                    [-mag_dir GENOME_PATH] [-table_dir HGT_TABLE_PATH]
                    [-threads THREADS] [-sth STH]

HGT detection pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r1 READ1, --read1 READ1
                        Path to read1 FASTA.
  -r2 READ2, --read2 READ2
                        Path to read2 FASTA.
  -i SAMPLE_ID, --sample_id SAMPLE_ID
                        Sample ID.
  -o OUTPUT, --output OUTPUT
                        Working directory.
  -mag_dir GENOME_PATH, --genome_path GENOME_PATH
                        Directory containing genomes.
  -table_dir HGT_TABLE_PATH, --HGT_table_path HGT_TABLE_PATH
                        Path to the HGT table.
  -threads THREADS, --threads THREADS
                        Number of threads.
  -sth STH, --sth STH   reads span sites number
```



#### Input



#### Output





#### Example

```
python HGTfinder.py -r1 $read1 -r2 $read2 -i ${i} -o ./result -mag_dir genome_folder -table_dir output/New_HGT_table1.csv -threads 5
```









