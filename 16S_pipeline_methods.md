# 16S rRNA Amplicon Sequence Processing Pipeline

## 1. Methodology
Raw FASTQ files underwent quality control via FastQC (v0.11.9) (Andrews, 2010) and were processed using the DADA2 pipeline (v1.28.0) (Callahan et al., 2016) in the R statistical environment (R Core Team, 2023). Forward and reverse reads were truncated at 220 bp to remove low-quality tails. Following dereplication, reads were denoised to infer exact Amplicon Sequence Variants (ASVs), merged, and purged of chimeras using the consensus method. Taxonomy was assigned against the SILVA reference database (release 138.1) (Quast et al., 2013). To ensure data integrity, putative contaminants were identified via negative controls and removed using the decontam R package (v1.20.0) (Davis, Proctor, Holmes, Relman, & Callahan, 2018), while a defined mock community was processed in parallel to validate pipeline accuracy. For community analysis, the filtered ASV table, taxonomic assignments, and metadata were integrated into a phyloseq object (v1.44.0) (McMurdie & Holmes, 2013). 

**Required Tabular Data:**
The metadata required to execute this process (including sample IDs, negative control indicators, and mock community labels) is provided in the accompanying file: `try_paired_end_template_nephele.xlsx`.

---

## 2. Custom Code

### 2.1. Quality Control (Bash)
```bash
# Create directory for FastQC results
mkdir -p fastqc_results

# Run FastQC on all raw FASTQ files
fastqc *.fastq.gz -o fastqc_results/ -t 4
```

### 2.2. DADA2, Decontam, and Phyloseq Pipeline (R)
```R
# ==========================================
# 0. Load Libraries & Setup
# ==========================================
library(dada2)     # v1.28.0
library(decontam)  # v1.20.0
library(phyloseq)  # v1.44.0
library(ggplot2)
library(readxl)

# Define path to FASTQ files
path <- "path/to/your/fastq/files"

# Get forward and reverse read file names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names from filenames
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# ==========================================
# 1. Filter and Trim
# ==========================================
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Truncate at 220 bp and filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(220, 220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# ==========================================
# 2. Learn Error Rates
# ==========================================
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# ==========================================
# 3. Denoise and Infer ASVs
# ==========================================
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# ==========================================
# 4. Merge Paired Reads
# ==========================================
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

# ==========================================
# 5. Remove Chimeras
# ==========================================
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# ==========================================
# 6. Assign Taxonomy
# ==========================================
taxa <- assignTaxonomy(seqtab.nochim, "path/to/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# ==========================================
# 7. Integrate into Phyloseq
# ==========================================
# Load Metadata from the provided tabular data file
metadata_file <- "path/to/try_paired_end_template_nephele.xlsx"
samdf <- as.data.frame(read_excel(metadata_file))

# Set rownames to match sample IDs (assuming Sample IDs are in the first column)
rownames(samdf) <- samdf[, 1]

# Construct Phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# ==========================================
# 8. Identify & Remove Contaminants (Decontam)
# ==========================================
# Define negative controls based on metadata (adjust "Sample_Type" as necessary)
sample_data(ps)$is.neg <- sample_data(ps)$Sample_Type == "Negative_Control"

# Identify contaminants using prevalence method
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

# Prune contaminants from the dataset
ps.clean <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.clean <- subset_samples(ps.clean, Sample_Type != "Negative_Control")

# ==========================================
# 9. Mock Community Validation
# ==========================================
ps.mock <- subset_samples(ps.clean, Sample_Type == "Mock")
ps.mock <- prune_taxa(taxa_sums(ps.mock) > 0, ps.mock)

mock.taxa <- tax_table(ps.mock)
print(mock.taxa)

# ==========================================
# 10. Export Data (ASV and Taxonomy Tables)
# ==========================================
# Save ASV table as a .csv file
write.csv(seqtab.nochim, "ASV_table.csv", quote=FALSE)

# Save Taxonomy table as a .csv file
write.csv(taxa, "Taxonomy_table.csv", quote=FALSE)
```
