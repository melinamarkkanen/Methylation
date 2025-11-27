# Methylation project

## *Preprint:*

## *Graphical abstract*

## Contents
- [Data](#Data)
- Analysis steps:
  - [Synthetic community](./HAMBI_analyses/README.md)
  - [Wastewater community](./WW_analyses/README.md)
- HiFi read generation and metagenomic assembly:
  -  [Synthetic community](./workflow/Snakefile_HAMBI_preanalysis)
  -  [Wastewater community](./workflow/Snakefile_WW_preanalysis)
- Preliminary analysis for methylation detection:
  -  [Synthetic community](./workflow/Snakefile_HAMBI_methylation_analysis)
  -  [Wastewater community](./workflow/Snakefile_WW_methylation_analysis)
- Generation of Position Weight Matrices (PWMs)
  -  [Synthetic community](./workflow/Snakefile_HAMBI_preanalysis)
  -  [Wastewater community](./workflow/Snakefile_WW_preanalysis)
- Analyses and visualization of PWMs
  - Synthetic community:
    - [Random Forest classifier](./notebooks/Random_Forest_analysis_HAMBI.ipynb)
    - [Sequence logos](?)
    - [Uniform Manifold Approximation and Projection (UMAP)](./HAMBI_analyses/UMAP.md)
      - [Generation of attached data](./HAMBI_analyses/UMAP.md#generating-attached-data)
      - [Binning clustered contigs](./HAMBI_analyses/UMAP.md#binning-clustered-contigs)
  - Wastewater community:
    - [Uniform Manifold Approximation and Projection (UMAP)](./WW_analyses/UMAP.md)
      - [Generation of attached data](./WW_analyses/UMAP.md#generating-attached-data)
      - [Binning clustered contigs](./WW_analyses/UMAP.md#binning-clustered-contigs)
- Genetic context analyses:
  - [Latent class D 2 beta-lactamase](./WW_analyses/Genetic_context_analyses.md#class-d-2-beta-lactamase)
  - [Latent class C beta-lactamase](./WW_analyses/Genetic_context_analyses.md#class-c-beta-lactamase)
  - [*bla*OXA-129](./WW_analyses/Genetic_context_analyses.md#blaoxa-129)
  - [*sul1*](./WW_analyses/Genetic_context_analyses.md#sul1)
  - [*erm*(F)](./WW_analyses/Genetic_context_analyses.md#ermf)


## Data
### Synthethic community ("HAMBI") of 36 species by Partanen *et al.* 2025 [🔗](https://doi.org/10.1093/ismeco/ycaf113)
#### PacBio HiFi metagenomic long-reads
| Sample code | Sample | Description |
| -------------------- | ------------ | ------------ |
| bcAd1023T | tc_0.2_15_5  | 0.2 µg/ml tetracycline, 15 °C   |
| bcAd1037T | su_20_25_3   | 20 µg/ml sulphamethazine, 25 °C |
| bcAd1039T | su_20_25_4   | 20 µg/ml sulphamethazine, 25 °C |
| bcAd1046T | tc_0.02_37_3 | 0.02 µg/ml tetracycline, 37 °C  |
| bcAd1063T | no_0_37_3    | no antibiotic, 37 °C            |

### Synthethic community (HAMBI) of 36 species by Partanen *et al.* 2025 [🔗](https://doi.org/10.1093/ismeco/ycaf113)
### Assemblies for available WGS data (n=18 strains) by Hogle *et al.* 2024 [🔗](https://doi.org/10.1128/mra.00111-24)

| Strain ID | Species | Element |
| ------------- | ------------- | ------------- |
| HAMBI_0006 | *Pseudomonas_E putida* | chromosome |
| HAMBI_0105 | *Agrobacterium tumefaciens* | chromosome 1 |
| HAMBI_0105 | *Agrobacterium tumefaciens* | chromosome 2 |
| HAMBI_0105 | *Agrobacterium tumefaciens* | plasmid|
| HAMBI_0262 | *Brevundimonas bullata* | chromosome* |
| HAMBI_0403 | *Comamonas testosteroni_C* | chromosome |
| HAMBI_0403 | *Comamonas testosteroni_C* | plasmid|
| HAMBI_1287 | *Citrobacter_B koseri* | chromosome |
| HAMBI_1287 | *Citrobacter_B koseri* | plasmid1 |
| HAMBI_1287 | *Citrobacter_B koseri* | plasmid2 |
| HAMBI_1292 | *Morganella morganii* | chromosome |
| HAMBI_1299 | *Kluyvera intermedia* | chromosome |
| HAMBI_1299 | *Kluyvera intermedia* | plasmid|
| HAMBI_1842 | *Sphingobium yanoikuyae* | chromosome |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid1* |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid2 |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid3 |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid4 |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid5* |
| HAMBI_1842 | *Sphingobium yanoikuyae* | plasmid6 |
| HAMBI_1896 | *Sphingobacterium spiritivorum* | chromosome |
| HAMBI_1972 | *Aeromonas caviae* | chromosome |
| HAMBI_1977 | *Pseudomonas_E chlororaphis* | chromosome |
| HAMBI_2159 | *Trinickia caryophylli* | chromosome1* |
| HAMBI_2159 | *Trinickia caryophylli* | chromosome2* |
| HAMBI_2160 | *Bordetella avium* | chromosome |
| HAMBI_2164 | *Cupriavidus oxalaticus* | chromosome1 |
| HAMBI_2164 | *Cupriavidus oxalaticus* | chromosome2 |
| HAMBI_2443 | *Paracoccus denitrificans* | chromosome |
| HAMBI_2443 | *Paracoccus denitrificans* | plasmid1 |
| HAMBI_2443 | *Paracoccus denitrificans* | plasmid2 |
| HAMBI_2494 | *Paraburkholderia kururiensis_A* | chromosome* |
| HAMBI_2494 | *Paraburkholderia kururiensis_A* | plasmid* |
| HAMBI_2659 | *Stenotrophomonas maltophilia* | chromosome |
| HAMBI_3237 | *Microvirga lotononidis* | chromosome |
| HAMBI_3237 | *Microvirga lotononidis* | plasmid1 |
| HAMBI_3237 | *Microvirga lotononidis* | plasmid2 |
| HAMBI_3237 | *Microvirga lotononidis* | plasmid3* |
#### * not detected by BLASTn
#### ** no methylation signals

### Wastewater
#### PacBio HiFi metagenomic long-reads
| Sample code | Sample |
| ------------- | ------------- |
| EFF1 | effluent |
| EFF2 | effluent |
| EFF3 | effluent |
| INF1 | influent |
| INF2 | influent |
| INF3 | influent |
| SLU1 | sludge |
| SLU2 | sludge |
| SLU3 | sludge |