# Bioinformatics Pipelines for Nostoc co-culture medium optimization for phycocyanin production
16S rRNA processing and comparative genomics workflows for Nostoc co-culture.

This repository contains the custom scripts, data processing pipelines, and bioinformatic methodologies used in the manuscript:
**"Scalable optimization of phycocyanin production in a photosynthetic co-culture using high-throughput definitive screening design."**

## Repository Contents

This repository includes the following documentation and code:

*   **`16S_pipeline_methods.md`**: Details the 16S rRNA amplicon sequence processing pipeline. It includes R scripts and Bash commands for quality control via FastQC, ASV inference using DADA2, and taxonomy assignment against the SILVA reference database. It also covers contaminant removal using the `decontam` package and data integration via `phyloseq`.
*   **`Comparative_Genomics_Methods.md`**: Outlines the comparative genomics workflow. It includes custom Python scripts designed to evaluate the metabolic cross-feeding potential among three metagenome-assembled genomes (MAGs) affiliated with *Nostoc*, *Erythrobacter*, and *Allorhizobium*. The provided pipeline demonstrates how to filter KEGG annotations and construct a presence/absence matrix focused on cobalamin (vitamin B12) biosynthesis.

## Software Dependencies

To execute the pipelines described in the provided documents, the following environments and packages are required:

**R Environment (16S Amplicon Pipeline)**
*   R statistical environment
*   `dada2` (v1.28.0)
*   `phyloseq` (v1.44.0)
*   `decontam` (v1.20.0)
*   `ggplot2`
*   `readxl`

**Python Environment (Comparative Genomics)**
*   Python 3
*   `pandas`
*   `numpy`

## Data Availability

The sequencing data associated with this study have been deposited in the NCBI Sequence Read Archive (SRA):
*   **Raw 16S rRNA sequencing files:** BioProject PRJNA1513703
*   **Raw PacBio sequencing reads and assembled MAGs:** BioProject PRJNA1405787
