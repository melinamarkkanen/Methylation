# Random Forest analysis
&nbsp;
## HAMBI
### Number of contigs and analyzed for methylation
| sample | contigs | GFF files | label 'NA' |
|------------------------|-------|----------|-------------------|
| bcAd1023T   | 305   | 300      | 90                |
| bcAd1037T   | 336   | 327      | 123               |
| bcAd1039T   | 713   | 698      | 171               |
| bcAd1046T   | 263   | 262      | 159               |
| bcAd1063T   | 373   | 364      | 126               |
&nbsp;
### Unique taxa
|---------| All data | After filtering >50 |
|---------|----------|-----------------|
| Species |    15    | 11 |
| Genus   |    15    | 11 |
| Family  |    12    | 9  |
| Order   |    9     | 7  |
| Class   |    3     | 2  |
| Phylum  |    2     | 1  |
&nbsp;
## Random Forest
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.870786 | 0.893405 | 0.870786 | 0.858709 |
| **Domain to Genus**   | 1.0 | 0.870786 | 0.893405 | 0.870786 | 0.858709 |
| **Domain to Family**  | 1.0 | 0.854748 | 0.834308 | 0.854748 | 0.834842 |
| **Domain to Order**   | 1.0 | 0.866666 | 0.898901 | 0.866666 | 0.836683 |
| **Domain to Class**   | 1.0 | 0.934065 | 0.936121 | 0.934065 | 0.930331 |
| **Domain to Phylum**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
&nbsp;
## Random Forest **Exploring sufficient modification data amount**:
### Above 50 lines in .gff
#### Number of Taxa: 11
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.884353 | 0.863946 | 0.884353 | 0.869720 |
| **Domain to Genus**   | 1.0 | 0.884353 | 0.863946 | 0.884353 | 0.869720 |
| **Domain to Family**  | 1.0 | 0.918367 | 0.896655 | 0.918367 | 0.903183 |
| **Domain to Order**   | 1.0 | 0.905405 | 0.911239 | 0.905405 | 0.895936 |
| **Domain to Class**   | 1.0 | 0.953020 | 0.953301 | 0.953020 | 0.951368 |
| **Domain to Phylum**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
&nbsp;
### Above 100 lines in .gff
#### Number of Taxa: 10
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.946902 | 0.950707 | 0.946902 | 0.944577 |
| **Domain to Genus**   | 1.0 | 0.946902 | 0.950707 | 0.946902 | 0.944577 |
| **Domain to Family**  | 1.0 | 0.955752 | 0.954380 | 0.955752 | 0.951657 |
| **Domain to Order**   | 1.0 | 0.931034 | 0.884082 | 0.931034 | 0.905217 |
| **Domain to Class**   | 1.0 | 0.948275 | 0.951259 | 0.948275 | 0.943871 |
| **Domain to Phylum**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
&nbsp;
### Above 200 lines in .gff
#### Number of Taxa: 6
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.987654 | 0.989711 | 0.987654 | 0.986980 |
| **Domain to Genus**   | 1.0 | 0.987654 | 0.989711 | 0.987654 | 0.986980 |
| **Domain to Family**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
| **Domain to Order**   | 1.0 | 0.965517 | 0.969172 | 0.965517 | 0.961190 |
| **Domain to Class**   | 1.0 | 0.977528 | 0.978076 | 0.977528 | 0.976262 |
| **Domain to Phylum**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
&nbsp;
### Above 500 lines in .gff
#### Number of Taxa: 4
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
| **Domain to Genus**   | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
| **Domain to Family**  | 1.0 | 0.984615 | 0.985164 | 0.984615 | 0.984643 |
| **Domain to Order**   | 1.0 | 0.984615 | 0.985164 | 0.984615 | 0.984643 |
| **Domain to Class**   | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
| **Domain to Phylum**  | 1.0 | 1.0      | 1.0      | 1.0      | 1.0      |
&nbsp;
&nbsp;
&nbsp;
&nbsp;
## Mapping samples bcAd1039, bcAd1046 and bcAd1063 to Paracoccus GCF_034627565.1_ASM3462756v1
|  | bcAd1039 | bcAd1046 | bcAd1063 |
| ------------------------------------- | -------- | -------- | -------- |
| **Mapped Reads**                          | 15390    | 39375    | 82361    |
| **Alignments**                            | 19877    | 42325    | 87874    |
| **Mapped Bases**                          | 49095177 | 308123371| 733049617|
| **Mean Gap-Compressed Sequence Identity** | 86.9052% | 96.1005% | 96.5226% |
| **Max Mapped Read Length**                | 30408    | 38483    | 40843    |
| **Mean Mapped Read Length**               | 2469.95  | 7279.94  | 8342.05  |
&nbsp;
&nbsp;
&nbsp;
&nbsp;
## Wastewater
### Number of contigs and analyzed for methylation
| sample | contigs | .gff files | after filtering 50 lines in .gff |
|--------|---------|------------|----------------------------------|
| EFF1   | 60537   | 60374      | 60373                            |
| EFF2   | 59860   | 59697      |                 |
| EFF3   | 54991   | 54838      |                 |
| INF1   | 65422   | 64745      |                 |
| INF2   | 62366   | 62229      |                 |
| INF3   | 62921   | 62795      |                 |
| SLU1   | 74572   | 74388      |                 |
| SLU2   | 76311   | 76128      |                 |
| SLU3   | 73811   | 73646      |                 |

