# Comparative Genomic Analysis of Vitamin B12 Metabolism

## 1. Methodology
To evaluate the potential for metabolic cross-feeding within the co-culture, the genomes of the three assembled MAGs (*Nostoc*, *Erythrobacter*, and *Allorhizobium*) were systematically assessed for their metabolic capacities, with a specific focus on cobalamin (vitamin B12) biosynthesis. Coding sequences (CDS) were predicted and subsequently annotated to assign KEGG Orthology (KO) identifiers. Genes responsible for relevant metabolic pathways were identified by searching the translated CDS against the KEGG database. Hit filtering was rigorously applied using an E-value cutoff of $\le 1 \times 10^{-5}$, a minimum query coverage of 70%, and a sequence identity of $\ge 30\%$. Finally, a comparative genomic matrix was constructed to map the presence and absence of specific KO identifiers across the three MAGs, cataloging the specific protein locus tags when present.

---

## 2. Custom Code (Python)
```python
import pandas as pd
import numpy as np

# ==========================================
# 1. Load or Simulate Raw Annotation Data
# ==========================================
# Assuming you have a raw BLAST/Annotation output file (e.g., "raw_kegg_annotations.tsv")
# Here, we simulate a mock dataframe to demonstrate the pipeline
data = {
    'MAG': ['Nostoc', 'Nostoc', 'Erythrobacter', 'Allorhizobium', 'Allorhizobium', 'Allorhizobium'],
    'Locus_Tag': ['lcl|JBWBQA01_310', 'lcl|JBWBQA01_307', 'lcl|JBWBPZ01_2390', 'lcl|JBWBPY01_1156', 'lcl|JBWBPY01_1983', 'lcl|JBWBPY01_1456'],
    'KO_ID': ['K02006', 'K02007', 'K00432', 'K00432', 'K00595', 'K00768'],
    'Evalue': [1e-10, 1e-20, 1e-45, 1e-50, 1e-30, 1e-25],
    'QCov': [95.0, 88.0, 99.0, 100.0, 75.0, 80.0],
    'PIdent': [45.5, 50.2, 60.1, 85.0, 40.0, 35.5]
}
df_raw = pd.DataFrame(data)

# Load the KO descriptions (this could also be loaded from a separate metadata CSV)
ko_desc_data = {
    'Metabolism_Group': ['01_Vitamin B12 Metabolism'] * 5,
    'KEGG_KO_ID': ['K00432', 'K00595', 'K00768', 'K02006', 'K02007'],
    'KO_Description': [
        'gpx, btuE, bsaA; glutathione peroxidase [EC:1.11.1.9]',
        'cobL-cbiET; precorrin-6B C5,15-methyltransferase',
        'E2.4.2.21, cobU, cobT; nicotinate-nucleotide...',
        'cbiO; cobalt/nickel transport system ATP-binding protein',
        'cbiM; cobalt/nickel transport system permease protein'
    ]
}
df_ko = pd.DataFrame(ko_desc_data)

# ==========================================
# 2. Hit Filtering
# ==========================================
# Apply thresholds: E-value <= 1e-5, Query Coverage >= 70%, Percent Identity >= 30%
filtered_df = df_raw[
    (df_raw['Evalue'] <= 1e-5) & 
    (df_raw['QCov'] >= 70.0) & 
    (df_raw['PIdent'] >= 30.0)
]

# ==========================================
# 3. Aggregate Multiple Hits
# ==========================================
# If a single MAG has multiple locus tags for the same KO, join them with a comma
aggregated_df = filtered_df.groupby(['MAG', 'KO_ID'])['Locus_Tag'].apply(lambda x: ', '.join(x)).reset_index()

# ==========================================
# 4. Construct Presence/Absence Matrix (Pivot)
# ==========================================
pivot_df = aggregated_df.pivot(index='KO_ID', columns='MAG', values='Locus_Tag')

# Define the target MAGs exactly as they should appear in the final table
target_mags = ['Nostoc', 'Erythrobacter', 'Allorhizobium']

# Ensure all target MAG columns exist; add them as NaN if a MAG has no hits at all
for mag in target_mags:
    if mag not in pivot_df.columns:
        pivot_df[mag] = np.nan

# Fill empty cells with 'Absent'
pivot_df = pivot_df[target_mags].fillna('Absent').reset_index()

# ==========================================
# 5. Merge with KEGG Metadata
# ==========================================
final_table = pd.merge(df_ko, pivot_df, left_on='KEGG_KO_ID', right_on='KO_ID', how='left')

# Clean up any remaining NaNs for KOs that were in the reference but not found in any MAG
for mag in target_mags:
    final_table[mag] = final_table[mag].fillna('Absent')

# Drop the redundant KO_ID column
final_table = final_table.drop(columns=['KO_ID'])

# ==========================================
# 6. Export Results
# ==========================================
output_file = "GenomeComparative.csv"
final_table.to_csv(output_file, index=False)
print(f"Table successfully saved to {output_file}")
```
