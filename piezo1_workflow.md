# Pipeline overview

This pipeline was applied independently for all taxonomic orders included in this study, besides sections indicated by "[ ]", however this document focuses on the primate order.

The main workflow consists of:

1. Data collection (NCBI, UniProt, DNAZoo).
2. Quality assessment of protein sequences and isoform selection.
3. De novo prediction of PIEZO1 sequences
4. Dataset validation
5. Extraction of QE regions
6. Correlational analysis



# 1. Data collection

## Description

Script:`PIEZO1_sequence_collection.ipynb`

All available PIEZO1 ortholog sequences for each taxonomic order were collected from three databases: NCBI RefSeq, UniProt and DNAZoo.
NCBI is queried by gene name, taxon and filtered by RefSeq entries.
Uniprot also uses gene name and taxon, with the additional fixed parameters of an annotation score of at least 2 and level of `existence:3` (inferred from homolog<>).
For DNAZoo instead, all species names present in NCBI Taxonomy are used to pull inferred whole-proteome files (`HiC.fasta_v2.functional.proteins.fasta.gz`). These are downloaded, and filters records whose description contains the gene name PIEZO1.
Results are then merged into a single FASTA file per order, with source database tagged in each header.

## Parameters

Inputs: taxonomic order name (optional: protein name. Also usable for other proteins)

Output: Fasta file with all collected sequences, optional CSV summary with species, accession and source.

## Code used in workflow

```bash
python PIEZO1_sequence_collection.ipynb --taxon primates --csv
```

## Results

Primates total sequences: 158 (NCBI:103, UniProt:18, DNAZoo:36), 52 species
Rodents total sequences: 170 (NCBI:132, UniProt:6, DNAZoo:32), 56 species
Carnivora total sequences: 266 (NCBI:158, UniProt:3, DNAZoo:108), 68 species
Cetartiodactyla total sequences: 278 (NCBI:190, UniProt:12, DNAZoo:76) 70 species

A total of 246 species were represented in the initial dataset. Meaning that on average each ortholog presented around 3 alternative sequences.

# 2. Quality assessment and isoform selection of PIEZO1 orthologs

## 2.1. Description: Quality control algorithm

Script:`QA_protein_sequences.sh`

PIEZO1 ortholog sequences were then processed to filter out potentially erroneous sequences and to select a single version of PIEZO1 per species.

All sequences are blasted against the reference sequence and filtered based on their percentage identity (pident), query coverage (qcov), ratio length and whether a canonical methionine is present at the start of the sequence.
When multiple sequences were available for the same species, the best matching one was kept to avoid a redundant dataset. Selection was based on a scoring metric considering pident, qcov and a length ratio penalty.

Reference sequence selection: For each taxonomic order, a known protein was used (NP) if available. If not, a predicted model protein sequence (XP) was selected based on the quality of the support for the predicted transcript (XM) (preferably supported at 99-100% by RNA data, including support for introns).

### References used

|Taxon          | Accession       | Species                  |
|---------------|-----------------|--------------------------|
|Primates       | NP_001136336.2 | *Homo sapiens*           |
|Rodentia       | NP_001344278.1 | *Mus musculus*           |
|Carnivora      | XP_005620631.1 | *Canis lupus familiaris* |
|Cetartiodactyla| XP_015331270.2 | *Bos taurus*             |

### Parameters

**Inputs**: Multi-FASTA file of protein sequences and accession number of reference

**Outputs**:

```
(input_file_name)_QA_results
|   |-- INTERMEDIATE_results
|   |   |-- failed_species.txt
|   |   |-- final_selected_sequences.fasta
|   |   `-- qa_results.tsv
|   |-- LIGHT_results
|   |   |-- failed_species.txt
|   |   |-- final_selected_sequences.fasta
|   |   `-- qa_results.tsv
|   |-- STRICT_results
|   |   |-- failed_species.txt
|   |   |-- final_selected_sequences.fasta
|   |   `-- qa_results.tsv
|   |-- quality_thresholds.tsv
|   |-- blastp_raw_results.tsv
|   `-- sensitivity_summary.tsv

```

### Code used in workflow

```bash
bash ./QA_protein_sequences.sh ./piezo1_primates_available_sequences.fasta NP_001136336.2
```

## 2.2. Validation of threshold-based filtering

### 2.2.1. Sensitivity analysis

To determine the robustness of quality assessment, three thresholds were set for percentage identity, query coverage and length ratio. Stable numbers of species retention across thresholds suggest that filtering criteria are not too sensitive to parameter value.

Primates

| Threshold    | Species_with_Sequences | Species_Missing | Species_Total |
|--------------|------------------------|-----------------|---------------|
| LIGHT        | 33                     | 19              | 52            |
| INTERMEDIATE | 31                     | 21              | 52            |
| STRICT       | 31                     | 21              | 52            |

Rodents

| Threshold    | Species_with_Sequences | Species_Missing | Species_Total |
|--------------|------------------------|-----------------|---------------|
| LIGHT        | 46                     |11               | 56            |
| INTERMEDIATE | 43                     |14               | 56            |
| STRICT       | 41                     |16               | 56            |

Carnivora

| Threshold    | Species_with_Sequences | Species_Missing | Species_Total |
|--------------|------------------------|-----------------|---------------|
| LIGHT        | 39                     | 29              | 68            |
| INTERMEDIATE | 38                     | 30              | 68            |
| STRICT       | 35                     | 33              | 68            |

Cetartiodactyla

| Threshold    | Species_with_Sequences | Species_Missing | Species_Total |
|--------------|------------------------|-----------------|---------------|
| LIGHT        | 45                     | 25              | 70            |
| INTERMEDIATE | 45                     | 2              | 70            |
| STRICT       | 44                     | 26              | 70            |


The small variations between different thresholds indicates a consistent and robust filtering. The STRICT threshold set of sequences were used for downstream anaylsis. This conservative choice was made considering the possibility of retrieving failed sequences via *de novo* prediction in Section 3.

Preliminary ortholog dataset was composed of 157 sequences (species) across all 4 orders.

### 2.2.2. Description: Pairwise comparison of quality metric thresholds

Script: `QA_metric_plot_analysis.Rmd`

To determine whether the passed and failed sequences correspond to two real and distinct populations a pairwise comparisons between all metric was carried out using `ggpairs` function from GGally R package (Schloerke B *et al*, 2025). 

Pairwise scatterplots, per-metric density distributions and correlation coefficients were generated, with points coloured by QA status (PASS/FAIL) to assess the separation between retained and discarded sequences along each metric.

#### Parameters

**Inputs:**

Merged file of TSV output file from QA script for all three sensitivity levels (`QA_results_merged_<order>.tsv`)

**Outputs:**

A matrix of scatterplots, density plots and correlation tests between metrics (`QA_metrics_assessment_<order>.png`)


# 3. De novo prediction of Piezo1 sequences

## 3.1. Description: genome assemblies

Genome of species whose initial collected sequences failed the QA step were downloaded from NCBI or DNAZoo, depending on whether an assembly was available in either one with a preference for NCBI assemblies. Additionally, some genome assemblies available on NCBI for species that were not initially present in the dataset were also included to further expand the dataset.
All assemblies from NCBI were fetched using the ncbi "datasets" tool.

### Code used in workflow

```bash
# Download genomes of interest from NCBI
datasets download genome accession --inputfile NCBI_primate_list.txt
```

## 3.2. Description: references

Each order was divided into taxa that best separated all the selected species into big enough groups in order to use best suited references for inference, but keeping the amount of groups to a pragmatic amount. For primates the taxa were hominoidea, hylobatidae, cercopithecinae, colobinae, lemuriformes, lorisiformes, and platyrrhini and tarsiidaes as one single group.

Annotation and genome assemblies of reference sequences for each taxa were fetched using the ncbi "datasets" tool. Then the specific annotation of Piezo1 was extracted using command line tools.

Where no annotated reference was available for a given clade (e.g. Tarsiidae for primates), species were paired with the reference of the closest available sister clade (Platyrrhini in this case). This grouping is specific to the prediction step; the clade annotations used in downstream visualisations (Section 4) follow a taxonomic partition based on NCBI Taxonomy.

### Code used in workflow

```bash
# Download genomes and annotations for primate references from NCBI

datasets download genome accession --inputfile NCBI_ann_<order>_GeMoMa.txt --include genome,gff3

# Extract Piezo1 annotation 

cat <accession>_<assembly_name>_genomic.gff | grep "XM_012812057.2" > <clade>_piezo1_annotation.gff

cat <accession>_<assembly_name>_genomic.gff | grep "ID=gene-PIEZO1" >> <clade>_piezo1_annotation.gff
```

## 3.3. Description: GeMoMa protein prediction

Script: `run_GeMoMa_multiple_genomes.sh`

New primary sequences were predicted using the homology-based algorithm GeMoMa (version 1.9) \(J. Keilwagen *et al* 2016, 2018\).

For each order, target genomes were organised into the clade-level groups and each clade was paired with its corresponding reference genome and Piezo1-specific GFF annotation file. The GeMoMaPipeline module was run parallel across all targets using GNU `parallel` (O. Tange 2018). To automate the clade-specific reference selection, target genomes and their corresponding reference genome-annotation pairs were organised into matching clade subdirectories (e.g. `hominoidea/`, `cercopithecinae/`), allowing each target to be paired programmatically with the reference of its assigned clade. The maximum intron length parameter (`GeMoMa.m`) was increased from default to 50'000 bp to account for long introns of the PIEZO1 genomic locus and avoid fragmentation.

### Parameters

**Inputs**: Taxonomic order name. However correct folder set-up is required.

References folder set-up

reference/
|-- primates/
|   |-- cercopithecinae
|   |   |-- cercopithecinae_GCF_049350105.2_Macaca_mulatta.fna
|   |   |-- cercopithecinae_piezo1_annotation.gff
|   |-- colobinae
|   |-- hominoidea
|   |-- hylobatidae
|   |-- lemuriformes
|   |-- lorisiformes
|   |-- platyrrhini_tarsiidae
|-- cetartiodactyla
|-- carnivora
|-- rodentia

Target genomes folder set-up

genomes
|-- primates/
|   |-- cercopithecinae
|   |   |-- GCA_000000001.1_species_name.fna
|   |   |-- Species_name_HiC.fna
|   |   |-- ...
|   |-- colobinae
|   |-- hominoidea
|   |-- hylobatidae
|   |-- lemuriformes
|   |-- lorisiformes
|   |-- platyrrhini_tarsiidae
|-- cetartiodactyla
|-- carnivora
|-- rodentia

**Outputs**:

results_<TAXON>/
|-- <species>_GeMoMa_results/
|   |-- final_annotation.gff
|   |-- predicted_proteins.fasta
|   `-- protocol_*.txt
|-- <species>_error.log
~/piezo1_pipeline/GeMoMa/gemoma_manifest_<TAXON>.tsv

### Code used in workflow

```bash
./run_GeMoMa_multiple_genomes.sh primates
```
## 3.4. Description: QA of new sequences and final Piezo1 ortholog dataset

All primary sequences predicted using GeMoMa were collected into as single multi-fasta file for each order. The quality of the resulting sequences was assessed using the same script and rationale as for the initial collection. The GeMoMa file of each order with the passed sequences was then merged into a new file with the corresponding one with sequences from NCBI, UniProt and DNAZoo.

### Code used in worflow

```bash

# Gather results for GeMoMa

find ./results_<order> -maxdepth 1 -type d -name "*_GeMoMa_results" | sort | parallel '
  sp=$(basename {} _GeMoMa_results | sed -E "s/^(GCA|GCF)_[0-9]+\.[0-9]+_//; s/_HiC$//")
  genus=$(echo "$sp" | cut -d_ -f1)
  species=$(echo "$sp" | cut -d_ -f2)
  awk -v g="$genus" -v spec="$species" '\''
    /^>/ {match($0, /aa=[0-9]+/); print ">" g "_" spec " " substr($0, RSTART, RLENGTH) " [" g " " spec "]"; next} 1
  '\'' {}/predicted_proteins.fasta
' | awk '
  /^>/ {
    split($0, a, " ")
    sp = substr(a[1], 2)
    count[sp]++
    sub(/^>[^ ]+/, ">" sp "_pred" count[sp])
  }
  {print}
' > all_GeMoMa_<order>_sequences.fasta

# Quality Assessment

./QA_protein_sequences.sh ./GeMoMa/all_GeMoMa_<order>_sequences.fasta NP_001136336.2

# Add IDs to predictions
 awk '/^>/ {n++; printf ">GeP_%03d %s\n", n, substr($0,2); next} 1' ./QA_protein_sequences/all_GeMoMa_<order>_sequences.fasta_results/STRICT_results/final_selected_sequences.fasta > ./GeMoMa/<order>_GeMoMa_final_piezo1_predictions.fasta

# Merge GeMoMa and UniProt, NCBI and DNAZoo files

# Get GeMoMa predictions
cat ~/piezo1_pipeline/GeMoMa/<order>_GeMoMa_final_piezo1_predictions.fasta > ~/piezo1_pipeline/final_piezo1_sequences/<order>_piezo1_set.fasta

# Add NCBI, UniProt and DNAZoo sequences
cat ~/piezo1_pipeline/QA_protein_sequences/piezo1_<order>_available_sequences.fasta_results/STRICT_results/final_selected_sequences.fasta >> ~/piezo1_pipeline/final_piezo1_sequences/<order>_piezo1_set.fasta
```
# 4. Dataset annotation and validation

This section covers two steps to qualitatively validate the final ortholog set produced in Section 3.4. First, is necessary to assign each species to a taxonomic clade for downstream annotation, and then validate the *de novo* GeMoMa predictions tnrough pairwise sequence identity analysis.

## 4.1. Description: Clade assignment of species

To group species into taxonomically meaningful clades for downstream visualisation and analysis, full taxonomic lineafes were retrieved from NCBI Taxonomy via Entrez Direct (`esearch`/`efetch`/`xtract`; Kans, 2025). The lineage strings were manually inspected at successive taxonomic levels to identify a partition that separated the sequences into a small number of well-populated clades and with simple interpretability.

For Carnivora, Rodentia and Cetartiodactyla, a single lineage level produced a satisfactory partition (e.g. Canidae, Musteloidea, Pinnipedia, Ursidae, Felidae, Eupleridae, Herpestidae, Hyaenidae and Viverridae for Carnivora). For primates, no single level partitioned the order satisfactorily — the most populated level placed two-thirds of species into a single clade (Simiiformes) — so clades were defined across two adjacent lineage levels, with the strepsirrhine groups (Lemuriformes, Lorisiformes, Chiromyiformes) collapsed into a single Strepsirrhini clade and Simiiformes broken down into Hominoidea, Cercopithecoidea and Platyrrhini, alongside Tarsiidae. Each species was then assigned to its corresponding clade, producing a two-column mapping file used to annotate species in downstream visualisations.

### Parameters:

**Inputs:** Species lists per order (`all_<order>_species_initial.tsv`); selected clade list per order (`<order>_selected_clades.tsv`).

**Outputs:** Two-column species-to-clade mapping file per order (`<order>_species_clade_mapping.tsv`).

### Code used in workflow

```bash
# Retrieve full lineage strings from NCBI Taxonomy 

cat ./data_analysis/all_<order>_species_initial.tsv | parallel -j 3 --delay 0.4 'esearch -db taxonomy -query {} | efetch -format xml | xtract -pattern Taxon -element ScientificName Lineage' > initial_<order>_taxonomy.tsv

# Inspect lineage strings to identify a level partitioning the order into clades

cat initial_<order>_taxonomy.tsv | cut -d';' -f25 | sort | uniq -c

# Assign each species to its selected clade

cat <order>_selected_clades.tsv | parallel "grep {} initial_<order>_taxonomy.tsv | awk -F'\t' -v g={} '{print \$1 \"\t\" g}' " > <order>_species_clade_mapping.tsv

```

## 4.2. Description: Generation of pairwise-identity matrix

Script: `pairwise_matrix.sh`

The merged sequence set (NCBI + UniProt + DNAZoo + GeMoMa) was used to build an all-against-all pairwise BLASTp databse. For each pair of species the percentage identity of reciprocal sequence-level hits were averaged into a single species-pair value. These were then assembled into a symmetric matrix.

### Parameters

**Inputs:** Final merged FASTA from Section 3.4 (`<order>_piezo1_dataset.fasta`)

**Outputs:** All-vs-all BLASTp pairwise output (`species_pairwise_matrix.tsv`)

## 4.3. Description: Validation of ortholog datasets

Script: `clustered_heatmap.rmd`

The matrix (Section 4.2) was converted to a dissimilarity matrix (100-identity) and visualised as a clustered heatmap using the `ComplexHeatmap` R package (Gu, 2022), with hierarchical clustering by average linkage. Species were annotated by clade-level taxonomic group (based on Section 4.1) via coloured bars and silhouettes from PhyloPic (`rphylopic` package) were added alongisde the heatmap for each taxonomic group. 

### Parameters

**Inputs:** All-vs-all BLASTp pairwise output (`species_pairwise_matrix.tsv`); clade-level taxonomy (`<order>_species_clade_mapping.tsv`)

**Outputs:**  Clade-annotated clustered heatmap with silhouette panel (`<order>_piezo1_heatmap_with_silhouettes.png`).



# 5. Extraction of QE regions

## 5.1. Description: Multiple Sequence Alignment

Validated ortholog datasets were used to carry out multiple sequence alignment (MSA) using MAFFT (v7.526) L-INS-i algorithm.

### Parameters

**Inputs:**  Final Piezo1 ortholog dataset (`<order>_piezo1_dataset.fasta`)

**Outputs:**  Multiple sequence alignment file (`<order>_final_MSA.fasta`)

### Code used in workflow

```bash
linsi `<order>_piezo1_dataset.fasta` `<order>_final_MSA.fasta`
```

## 5.2. Description: Identifying regions of interest

An initial visual inspection of an MSA carried out with some primate Piezo1 sequences revealed possible amino acid repeat, or low-complexity, regions (LCR). 
Three LCRs were initially identified in human PIEZO1 by submitting NP_001136336.2 to PlaToLoCo (Jarnot et al., 2020) and using the SEG-intermediate algorithm output.

SEG intermediate results:

>NCBI|NP_001136336.2 NP_001136336.2 piezo-type mechanosensitive ion channel component 1 [Homo sapiens]|738|755
EEQQEHQQQQQEEEEEE
>NCBI|NP_001136336.2 NP_001136336.2 piezo-type mechanosensitive ion channel component 1 [Homo sapiens]|1454|1475
AQAVLRRRQQEQEQARQEQAG
>NCBI|NP_001136336.2 NP_001136336.2 piezo-type mechanosensitive ion channel component 1 [Homo sapiens]|1881|1910
ARKGAAAIEAEDREEEEGEEEKEAPTGRE

Two interval definitions were used in parallel. The LCR-anchored intervals were defined by conserved amino-acid motifs flanking each SEG-identified Q/E-rich stretch in the cross-order MSA. The domain-extended intervals were defined by the boundaries of the corresponding cytoplasmic domain in the human PIEZO1 UniProt annotation (Q92508), which encloses each LCR. The LCR-anchored intervals are used as the primary analysis, with the domain-extended intervals reported as a pre-specified sensitivity check.

To define analogous regions in each of the four orders, the human LCRs were located within the primate MSA and conserved flanking motifs were identified by visual inspection in JalView (2.11.5.1) combined with consensus-frequency analysis across all four order MSAs. Anchor motifs were retained only if they were conserved at ≥50% consensus frequency in every order MSA, ensuring consistent interval definition across orders.

The shared anchor motifs identified were: PLL and KWG flanking Region 1; AYQAWVTN and RVLS flanking Region 2; and SLRF (or SIRF in rodentia) and AATDV flanking Region 3. The S↔I substitution at the second position of the Region 3 upstream anchor reflects a conservative substitution fixed across rodentia. UniProt Q92508 was consulted to confirm that all three regions fall within cytoplasmic loops or the C-terminal cytoplasmic domain (CTD) of PIEZO1. For Region 3 specifically, no close-in downstream motif was conserved at ≥50% across all four orders, so the AATDV anchor was used; this places Region 3 as the entire C-terminal cytoplasmic extension containing LCR3, consistent with the UniProt CTD annotation.

Motif anchors:

Region 1. PLL-KWG
Region 2. QTVL-RVLS
Region 3. S[L/I]R[F/L]-DILH

Intervals:

# primate_intervals.txt
935:987:QE_region_1
1792:1858:QE_region_2
2250:2335:QE_region_3

# rodentia_intervals.txt
878:933:QE_region_1
1612:1686:QE_region_2
2087:2205:QE_region_3

# carnivora_intervals.txt
825:858:QE_region_1
1565:1617:QE_region_2
2016:2238:QE_region_3

# cetartiodactyla_intervals.txt
1134:1164:QE_region_1
1853:1901:QE_region_2
2302:2392:QE_region_3

## 5.3. Description: Extracting LCR metrics for analysis

Script: `Extract_QE_regions.sh`

### Code used in workflow

```bash
# Get metrics of interest for LCRs from untrimmed MSA

./Extract_QE_regions.sh ./<order>_final_MSA.fasta ./<order>_intervals_untrimmed.txt

# Get metrics of interest for LCRs from trimmed MSA

./Extract_QE_regions.sh ./<order>_trimmed_MSA.fasta ./<order>_intervals_trimmed.txt
```

## 5.4. Description: Generating final datasets with body weight and QE metrics

Script: `dataset_generation.Rmd`

The per-region per-species QE metrics produced in section 5.4 (E_count, QE_sum, E_ratio, Neg_charge, Net_charge, and Net_charge_no_H, for each of the three regions) were merged with species-level life-history traits to produce the analysis ready dataset used in the correlational analysis. Two parallel datasets were generated using the same procedure, one based on the untrimmed MSAs and one based on the 10% occupancy-trimmed MSAs (Section 5.2), enabling comparison across the two interva definitions.

Body mass and family assignments were obtained from the COMBINE trait dataset (Soria et al. 2021). Since this dataset contains multiple entries per species for some traits, these were deduplicated by averaging log10-transformed body mass values. The resulting trait table was joined to the per-order QE table by species name.

Two categorical covariates were derived from the family assignment for use in the DPG-interaction models. A binary indicator (DPG_binary) flagged species belonging to families with reduced erythrocyte 2,3-DPG sensitivity of haemoglobin: Feliformia (Felidae, Viverridae, Herpestidae, Hyaenidae, Eupleridae, Nandiniidae) and Ruminantia (Bovidae, Cervidae, Giraffidae, Antilocapridae, Moschidae, Tragulidae). A three-level factor (DPG_non_binary) further distinguished Feliformia and Ruminantia as separate categories, with all other species coded as Sensitive, allowing slope and intercept differences between the two natural-control groups to be tested independently.

he TimeTree-derived phylogeny used in Section 6 carries tip labels that
reflect older mammalian taxonomy. To allow caper::comparative.data to match
trait records to tips, the data are renamed to the tree's nomenclature
(rather than the reverse, since the tree cannot be edited downstream).
Every rename below is a documented species-level synonymy in the ASM
Mammal Diversity Database (MDD; Burgin et al. 2018; Mammaldiversity.org,
accessed YYYY-MM-DD).
"the tree's tip labels reflect older synonymies; we harmonised the trait dataset to the tree's nomenclature so that comparative.data could match them, with the synonymies cross-referenced to the ASM Mammal Diversity Database (Burgin et al. 2018; mammaldiversity.org)
 Drops were unevenly distributed across orders, with Primates losing 4/75 species (5.3%) and Carnivora retaining the full sample. Within Cetartiodactyla, the three dropped species fell at the smaller end of the body-mass range (mean log10 body mass 4.71 vs 5.39 for retained), although the retained range remained fully bracketed by other species at both extremes.

Additionally, two haematological response variables relevant to the mechanistic hypothesis, mass-specific basal metabolic rate (MS-BMR, W/g) and haemoglobin oxygen affinity (P50, pH 7.2, 20°C), were obtained from Udroiu (2023) and joined to the dataset by species. These traits are sparsely sampled at the species level (approximately 92 of 276 species have BMR, with similar coverage for P50): However these were not used for analysis in this thesis.

Additionally, body mass distribution across the four orders was visually inspected on both raw and log10-transformed scales (Supplementary Fig. SX) to confirm that log10-transformation symmetrises the right-skewed raw distributions, justifying the use of log10_body_mass_g as the predictor in all PGLS models. Per-order × per-region summary statistics for all six QE metrics (n, mean, SD, min, max) are reported in Supplementary Table SX.¨


# 6. Correlational analysis

Datasets generated according to the description in Section 5.5. were used to test for linear relationships between metrics describing QE regions composition and mass of each correspondent species. Phylogenetic Generalised Least Squares (PGLS) was implemented to account to non-independence of species (data points) because of shared evolutionary history. This works by incorporating an expected covariance structure derived from phylogeny to an ordinary least squares regression framework.
The strength of phylogenetic signal was estimated from the data using Pagel's Lambda by maximum likelihood optimisation within bounded interval [0.01-0.99] for better interpretation. (look into this aspect)

Initially, relationship was investigated only within orders, however, a pooled version was carried out to investigate whether within order analysis provide to few independent events to actually produce significant results. Additionally, since mammalian orders have different ways of modulating hemoglobin affinity, interactions with 2,3-DPG-insensitivity were also tested to investigate its impact on haemoglobin affinity modulation.

## 6.1. Description: Per-order PGLS models

Scripts: `PGLS_analysis_primates.Rmd, PGLS_analysis_rodentia.Rmd, PGLS_analysis_carnivora.Rmd, PGLS_analysis_cetartiodactyla.Rmd`

For each order, a time-calibrated species-level phylogeny was generated from the list of complete species in all orders using TimeTree (`ll_species_timetree.nwk`). For individual order analysis this was pruned using `comparative.data` function from caper. The same set of linear models was fitted in every order, with `log10_body_mass_g`as the sole predictor and either one of the selected QE metrics as response: QE_sum, E_count, E_ration, Neg_charge, Net_charge and Net_charge_no_H for each region. Additionally metrics describing the sum of all regions were also tested as QE_sum_total and E_sum_total.

For each fitted model, the slope, intercept, ML-estimated lambda, standard error, 95% confidence interval, p-value and R-squared were saved into a per-order summary table.

The quality of fitted models was assessed using standard linear regression diagnostics (Mundry et al. 2014). Hence, for each model Q-Q plots of phylogenetic residuals, residuals-vs-fitted plots, lambda profile and leave-one-out refitting procedure were carried out and results were reported in table. each species was iteratively dropped, the model refitted with λ re-estimated by maximum likelihood, and standardised slope-change values (dfbetas) computed; species with |dfbeta| > 1 were flagged as influential

### Parameters

**Inputs:** Dataset with per-species trait values and QE values (`main_dataset.tsv`), per-order QE long-format (`<order>_QE_regions_290426.tsv`), clade mapping file (`<order>_species_clades_mapping.tsv`), time-calibrated tree (`all_species_timetree.nwk`).

**Output:** Per-order PGLS summary table (`<order>_PGLS_all_results.xlsx`); diagnostic plots (`diag_<metric>_<region>.png`), per-metric scatter plots with phylogentic regression lines and 95% confidence bands for eacc region and coloured by clade (`PGLS_<order>_<metric>.png`).

## 6.2. Description: Pooled four-order PGLS models

Script: `PGLS_analysis_Pooled_orders.Rmd`

To test whether negative results observed in within order analysis, PGLS for all four orders together was carried out. Some species were harmonised to match phylogeny tip labels by 

Same models were fitted as in Section 6.1 on the pooled data. Diagnostic was also carried out similarly to Section 6.1, with the addition of colouring residulas by order and reporting also the variances of residual by order to test for order-specific heteroscedascity. The LOO procedure was parallelised using the package `furrr`. 

### Parameters

**Inputs:** See section 6.1

**Outputs:** Pooled baseline PGLS summary table (`pooled_PGLS_all_results.xlsx`). Diagnostic panels (`diag_<metric>_<region>.png`); pooled PGLS plots as in Section 6.1 (`PGLS_pooled_<metric>.png`).

# 6.3. Description: Pooled PGLS with DPG-sensitivity interaction

Script: `PGLS_analysis_pooled_DPG.Rmd`

To test whether the pooled mass-QE relationship was modulated by 2,3-DPG sensitivity of hemoglobin, the same models were refitted on the pooled dataset with the three-level interaction term `DPG_factor` (Sensitive, Feliformia and Ruminantia), with sensitive group set as the reference level, so the slopes of DPG insensitive groups were tested to see if they changed compared to the sensitive group.
For each paramater the correspondent relevant values were reported as in Section 6.1 in a summary table.

Diagnositc panels match those of Section 6.2. but with residuals coloured and grouped by `DPG_factor` rather than order, and per-group residual variances reported separately for Sensitive, Feliformia and Ruminantia. For visualisation, separate confidence bands were computed for each DPG group from the model variance-covariance matrix and overlaid as three coloured regression lines per panel (solid for Sensitive, dashed for Feliformia, dotted for Ruminantia), faceted by region.

### Parameters

**Inputs:** Same pooled inputs as Section 6.2, with the additional requirement that DPG_factor is present and ordered with Sensitive as the reference level.

**Outputs:** Pooled DPG-interaction PGLS summary table (pooled_PGLS_DPG_all_results.xlsx); diagnostic panels (diag_<metric>_<region>_DPG.png); pooled DPG-interaction plots with three group-specific regression lines and confidence bands per region (PGLS_pooled_DPG_int_<metric>.png).