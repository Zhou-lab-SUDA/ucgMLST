# ucgMLST

**Universal Core Gene Multilocus Sequence Typing**

A comprehensive bioinformatics toolkit for bacterial genome analysis, taxonomic profiling, and phylogenetic inference using universal core genes (UCG) and multilocus sequence typing (MLST) approaches.


## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Dependencies](#dependencies)
  - [Setup](#setup)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Database Construction](#database-construction)
  - [Query Analysis](#query-analysis)
  - [Phylogenetic Analysis](#phylogenetic-analysis)
  - [Profiling and Compilation](#profiling-and-compilation)
- [Modules](#modules)
- [Input/Output Formats](#inputoutput-formats)
- [Examples](#examples)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

ucgMLST is a powerful toolkit designed for high-throughput bacterial genome analysis. It leverages universal core genes to perform taxonomic classification, strain-level identification, and phylogenetic reconstruction from both assembled genomes and metagenomic sequencing data.

The package integrates multiple approaches including:
- MinHash-based genome sketching for rapid similarity assessment
- Core gene identification and annotation
- Read mapping and taxonomic profiling
- Phylogenetic tree construction
- Multi-sample comparative analysis

## Features

- **Flexible Database Construction**: Build custom reference databases from NCBI, GTDB, or user-provided genomes
- **Rapid Taxonomic Profiling**: Identify bacterial species and strains from sequencing reads with ANI-based classification
- **Core Gene Analysis**: Identify and analyze universal single-copy core genes (USCG)
- **Phylogenetic Inference**: Generate phylogenetic trees from core gene alignments
- **Multi-sample Comparison**: Compile and compare taxonomic profiles across multiple samples
- **High Performance**: Multi-threaded processing with efficient memory usage
- **Comprehensive Outputs**: SAM/BAM alignments, JSON profiles, consensus sequences, and phylogenetic trees

## Installation

### Prerequisites

- Python 3.7 or higher
- Linux/Unix operating system
- At least 16GB RAM (32GB+ recommended for large databases)
- Sufficient disk space for reference databases

### Dependencies

**Python packages:**
```bash
numpy
pandas
click
ete3
pyarrow  # for feather format support
```

**External tools** (automatically detected if in PATH or `externals/` directory):
- [bindash](https://github.com/zhaoxiaofei/bindash) - MinHash-based genome distance estimation
- [minimap2](https://github.com/lh3/minimap2) - Sequence alignment
- [samtools](https://github.com/samtools/samtools) - SAM/BAM file processing
- [diamond](https://github.com/bbuchfink/diamond) - Protein sequence alignment
- [hmmsearch](http://hmmer.org/) - HMM-based sequence search (HMMER3)
- [getorf](http://emboss.sourceforge.net/) - ORF prediction (EMBOSS)
- [iqtree](http://www.iqtree.org/) - Phylogenetic tree reconstruction
- pigz or gzip - File compression

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/ucgMLST.git
cd ucgMLST
```

2. Install Python dependencies:
```bash
pip install numpy pandas click ete3 pyarrow
```

3. Download and place external tools in the `externals/` directory or ensure they are available in your system PATH.

4. Verify installation:
```bash
python configure.py
```

## Quick Start

```bash
# 1. Build a reference database
python build.py \
  -d mydb \
  -m bacteria \
  -g genomes.csv \
  -T bac120_taxonomy.tsv.gz \
  -G taxdump.tar.gz \
  assembly_summary.txt

# 2. Query sequencing reads against the database
python genoQuery.py \
  -d mydb \
  -r reads_R1.fastq.gz reads_R2.fastq.gz \
  -o results

# 3. View the taxonomic profile
cat results/profile.json
```

## Usage

### Database Construction

#### Building the Main Database (`build.py`)

Create a reference database from genome assemblies:

```bash
python build.py \
  --dbname <database_path> \
  --module <module_name> \
  --genome <genome_list.csv> \
  --cutoff 0.99 \
  --threads 80 \
  --gtdb bac120_taxonomy.tsv.gz \
  --genbank taxdump.tar.gz \
  assembly_summary_1.txt assembly_summary_2.txt
```

**Parameters:**
- `-d, --dbname`: Path to the database directory (required)
- `-m, --module`: Name of the module/dataset (required)
- `-g, --genome`: CSV file mapping accessions to genome files (required)
- `-c, --cutoff`: ANI cutoff for redundancy filtering (default: 0.99)
- `-t, --threads`: Number of threads (default: 80)
- `-T, --gtdb`: GTDB taxonomy file for bacterial classification
- `-G, --genbank`: NCBI taxonomy dump file
- `queries`: One or more NCBI assembly summary files

**Genome list format (`genome_list.csv`):**
```
GCF_000001405.40,/path/to/GCF_000001405.40.fna
GCF_000002985.6,/path/to/GCF_000002985.6.fna
```

#### Building Core Gene Database (`build_ucgDB.py`)

Extract and index universal core genes:

```bash
python build_ucgDB.py \
  --dbname <database_path> \
  --domain bacteria \
  --n_proc 40
```

**Parameters:**
- `-d, --dbname`: Path to the database directory
- `-D, --domain`: Taxonomic domain (bacteria, archaea, or eukaryota)
- `-n, --n_proc`: Number of parallel processes

#### Building User-Specific Database (`build_userDB.py`)

Add custom genomes to the database:

```bash
python build_userDB.py \
  --dbname <database_path> \
  --module <module_name> \
  --genomes genome1.fasta genome2.fasta \
  --n_proc 20
```

### Query Analysis

#### Taxonomic Profiling (`genoQuery.py`)

Identify bacterial species and strains from sequencing reads:

```bash
python genoQuery.py \
  --dbname <database_path> \
  --reads sample_R1.fastq.gz sample_R2.fastq.gz \
  --output results_dir \
  --n_proc 40 \
  --min_identity 0.93 \
  --min_prevalence 0.01
```

**Parameters:**
- `-d, --dbname`: Path to the database
- `-r, --reads`: Input FASTQ files (single or paired-end)
- `-o, --output`: Output directory
- `-n, --n_proc`: Number of processes (default: 40)
- `-i, --min_identity`: Minimum identity threshold (default: 0.93)
- `-p, --min_prevalence`: Minimum prevalence for reporting (default: 0.01)
- `--block_size`: Fragment size for analysis (default: 500)
- `--min_depth`: Minimum depth for consensus calling (default: 3)

**Outputs:**
- `profile.json`: Taxonomic profile with abundances and ANI values
- `primary.bam`: Read alignments to core genes
- Consensus sequences for detected species

#### Resolving Specific Taxa (`genoResolve.py`)

Perform detailed analysis of specific taxonomic groups:

```bash
python genoResolve.py \
  --dbname <database_path> \
  --reference "Escherichia coli" \
  --reads sample_R1.fastq.gz sample_R2.fastq.gz \
  --output results_dir \
  --genus
```

**Parameters:**
- `-d, --dbname`: Database path
- `-R, --reference`: Reference genome or taxon name
- `-r, --reads`: Input reads
- `-o, --output`: Output directory
- `-G, --genus`: Include all genomes from the same genus
- `--min_identity`: Minimum ANI (default: 0.93)

### Phylogenetic Analysis

#### Building Phylogenetic Trees (`genoPhylo.py`)

Construct phylogenetic trees from core genes:

```bash
python genoPhylo.py \
  --dbname <database_path> \
  --reference "GCF_000005845.2" \
  --reads sample1_R1.fq.gz sample1_R2.fq.gz \
  --output phylo_results \
  --n_proc 20
```

**Parameters:**
- `-d, --dbname`: Database path
- `-r, --reference`: Reference genome accession or taxon
- `-R, --reads`: Input sequencing reads
- `-o, --output`: Output directory
- `-n, --n_proc`: Number of processes
- `--min_identity`: Minimum alignment identity
- `-G, --genus`: Include genus-level genomes

**Outputs:**
- Core gene alignments
- Phylogenetic tree (Newick format)
- SNP matrices

#### Prepare Resolve Database (`build_resolveDB.py`)

Create a specialized database for high-resolution analysis:

```bash
python build_resolveDB.py \
  --dbname <main_database> \
  --reference "Salmonella enterica" \
  --outdir resolve_db \
  --genus
```

### Profiling and Compilation

#### Compile Multi-Sample Results (`genoCompile.py`)

Compare taxonomic profiles across multiple samples:

```bash
python genoCompile.py \
  --output comparison_report \
  --min_rpkm 0.001 \
  --min_ani 0.96 \
  sample1/profile.json sample2/profile.json sample3/profile.json
```

**Parameters:**
- `-o, --output`: Output prefix
- `-m, --min_rpkm`: Minimum RPKM threshold (default: 0.001)
- `-M, --min_ani`: Minimum ANI threshold (default: 0.96)
- `-p, --no_profile`: Skip profile comparison
- `-u, --otu`: Include OTU-level comparison
- `infiles`: Input profile JSON files

**Outputs:**
- `<output>.profile`: Taxa abundance matrix
- `<output>.otus`: OTU abundance matrix (if requested)

## Modules

### Core Modules

| Module | Description |
|--------|-------------|
| `build.py` | Main database construction with genome sketching and taxonomy integration |
| `build_ucgDB.py` | Universal core gene identification and indexing |
| `build_userDB.py` | User genome integration into existing databases |
| `build_resolveDB.py` | High-resolution database preparation for specific taxa |

### Analysis Modules

| Module | Description |
|--------|-------------|
| `genoQuery.py` | Taxonomic profiling from sequencing reads |
| `genoResolve.py` | Detailed strain-level analysis |
| `genoPhylo.py` | Phylogenetic tree construction |
| `genoCompile.py` | Multi-sample comparison and reporting |
| `genoEffectors.py` | Effector protein and virulence factor analysis |

### Support Modules

| Module | Description |
|--------|-------------|
| `configure.py` | Configuration and external tool management |
| `SRA_Funcs.py` | Read processing and output generation functions |

## Input/Output Formats

### Input Formats

**Genome List (CSV):**
```
accession,genome_path
GCF_000001405.40,/data/genomes/GCF_000001405.40.fna.gz
GCF_000002985.6,/data/genomes/GCF_000002985.6.fna
```

**Assembly Summary (TSV):**
Standard NCBI assembly summary format with columns:
- accession
- organism_name
- taxonomy (optional)
- assembly_level
- refseq_category
- relation_to_type_material

**Sequencing Reads:**
- FASTQ format (.fastq, .fq)
- Gzip compressed (.fastq.gz, .fq.gz)
- Single-end or paired-end

### Output Formats

**Profile JSON:**
```json
{
  "profile": [
    [RPKM, read_count, ANI, "species_name", "accession", "taxonomy"],
    ...
  ],
  "OTU": [
    [RPKM, read_count, ANI, "species_name", ["references"], "taxonomy"],
    ...
  ]
}
```

**BAM/SAM:**
Standard alignment format with core gene references

**Feather Database:**
Binary DataFrame format for efficient database storage

## Examples

### Example 1: Build a Database from NCBI Genomes

```bash
# Download GTDB taxonomy
wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz

# Download NCBI taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# Download assembly summaries
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# Create genome list from downloaded genomes
ls /data/genomes/*.fna.gz | awk -F'/' '{print substr($NF,1,15)","$0}' > genomes.csv

# Build database
python build.py \
  -d bacteria_db \
  -m refseq \
  -g genomes.csv \
  -T bac120_taxonomy.tsv.gz \
  -G taxdump.tar.gz \
  -t 80 \
  assembly_summary.txt

# Build core gene index
python build_ucgDB.py -d bacteria_db -D bacteria -n 40
```

### Example 2: Taxonomic Profiling of Metagenomic Samples

```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
  python genoQuery.py \
    -d bacteria_db \
    -r ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    -o results/${sample} \
    -n 40
done

# Compile results
python genoCompile.py \
  -o comparison \
  -m 0.01 \
  -M 0.95 \
  results/*/profile.json
```

### Example 3: Phylogenetic Analysis of E. coli Isolates

```bash
# Create specialized database
python build_resolveDB.py \
  -d bacteria_db \
  -r "Escherichia coli" \
  -o ecoli_db \
  --genus

# Generate phylogenetic tree
python genoPhylo.py \
  -d ecoli_db \
  -r GCF_000005845.2 \
  -R isolate1_R1.fq.gz isolate1_R2.fq.gz \
       isolate2_R1.fq.gz isolate2_R2.fq.gz \
  -o phylo_tree \
  -n 20
```

## Performance Considerations

- **Memory**: Database construction requires 16-32GB RAM; queries need 8-16GB
- **Storage**: Plan for ~100GB per 10,000 genomes (including sketches and core genes)
- **Threading**: Scale threads based on available CPU cores; optimal performance at 40-80 threads
- **Batch Processing**: For large sample sets, process in batches to manage resources

## Troubleshooting

**Common Issues:**

1. **Missing external tools**: Ensure all required tools are in PATH or `externals/` directory
2. **Memory errors**: Reduce batch size or increase available RAM
3. **No results**: Check minimum identity and prevalence thresholds
4. **Taxonomy mismatches**: Verify GTDB and NCBI taxonomy files are current

## Citation


## License

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions, issues, or contributions:

- **Issues**: [GitHub Issues](https://github.com/yourusername/ucgMLST/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/ucgMLST/discussions)

## Acknowledgments

ucgMLST integrates several excellent bioinformatics tools:
- bindash for genome sketching
- minimap2 for sequence alignment
- HMMER for profile searches
- DIAMOND for protein alignment
- IQ-TREE for phylogenetic inference

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

**Version**: 1.0.0
**Last Updated**: January 2026
