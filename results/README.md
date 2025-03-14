# Jotain
&nbsp;
## HAMBI
### Number of contigs and analyzed for methylation
| sample | contigs | .gff files | label 'NA' |
|------------------------|-------|----------|-------------------|
| bcAd1023T--bcAd1023T   | 305   | 300      | 90                |
| bcAd1037T--bcAd1037T   | 336   | 327      | 123               |
| bcAd1039T--bcAd1039T   | 713   | 698      | 171               |
| bcAd1046T--bcAd1046T   | 263   | 262      | 159               |
| bcAd1063T--bcAd1063T   | 373   | 364      | 126               |

### Number of lines in .gff
#### Before filtering: 1951
#### After filtering (No methylation detected or less than 20 detected since those are converted to 0 filled matrices): 1282

### Missing labels after filtering: 372

### Unique taxa at
#### Domain to species: 15
#### Domain to genus: 15
#### Domain to family: 12
#### Domain to order: 9
#### Domain to class: 3
#### Domain to phylum: 2




## Random Forest **All data**:
### All features
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.870787 | 0.871968 | 0.870787 | 0.859151 |
| **Domain to Genus**   | 1.0 | 0.848315 |          |          |          |
| **Domain to Family**  | 1.0 | 0.893855 |          |          |          |
| **Domain to Order**   | 1.0 | 0.905556 |          |          |          |
| **Domain to Class**   | 1.0 | 0.934066 |          |          |          |
| **Domain to Phylum**  | 1.0 | 1.0      |          |          |          |

#### Testing different number of features
| Number of features | Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|------|-----------------------|-----|----------|----------|----------|----------|
| 1    | **Domain to Species** | 0.522472 | 0.449438 | 0.474459 | 0.449438 | 0.40248 |
| 3    | **Domain to Species** | 1.0 | 0.662921 | 0.652217 | 0.662921 | 0.656572 |
| 5    | **Domain to Species** | 1.0 | 0.707865 | 0.709828 | 0.707865 | 0.708207 |
| 10   | **Domain to Species** | 1.0 | 0.848315 | 0.851571 | 0.848315 | 0.847413 |
| 20   | **Domain to Species** | 1.0 | 0.870787 | 0.872965 | 0.870787 | 0.867311 |
| 50   | **Domain to Species** | 1.0 | 0.882022 | 0.88702  | 0.882022 | 0.871105 |
| 60   | **Domain to Species** | 1.0 | 0.882022 | 0.900687 | 0.882022 | 0.86989  |
| 70   | **Domain to Species** | 1.0 | 0.893258 | 0.893939 | 0.893258 | 0.883861 |
| 80   | **Domain to Species** | 1.0 | 0.86989  | 0.877542 | 0.882022 | 0.870674 |
| 90   | **Domain to Species** | 1.0 | 0.88764  | 0.889939 | 0.889939 | 0.879    |
| 100  | **Domain to Species** | 1.0 | 0.870787 | 0.890658 | 0.870787 | 0.858113 |
| 120  | **Domain to Species** | 1.0 | 0.882022 | 0.898927 | 0.882022 | 0.867911 |
| 150  | **Domain to Species** | 1.0 | 0.893258 | 0.909739 | 0.893258 | 0.878497 |
| 180  | **Domain to Species** | 1.0 | 0.870787 | 0.898134 | 0.898134 | 0.855244 |
| 200  | **Domain to Species** | 1.0 | 0.848315 | 0.821896 | 0.848315 | 0.829591 |
| 250  | **Domain to Species** | 1.0 | 0.870787 | 0.871968 | 0.870787 | 0.859151 |

## Random Forest **Chromosome data**:
### All features
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|----------|----------|----------|----------|----------|
| **Domain to Species** | 0.752464 | 0.585152 | 0.700569 | 0.585152 | 0.618785 |
| **Domain to Genus**   | 0.748083 | 0.606986 | 0.712237 | 0.606986 | 0.637655 |
| **Domain to Family**  | 0.735229 | 0.698689 | 0.774919 | 0.698689 | 0.709561 |
| **Domain to Order**   | 0.842794 | 0.782608 | 0.823232 | 0.782608 | 0.754668 |
| **Domain to Class**   | 0.955240 | 0.913043 | 0.915420 | 0.913043 | 0.899911 |
| **Domain to Phylum**  | 1.0      | 1.0      | 1.0      | 1.0      | 1.0      |

#### Testing different number of features
| Number of features | Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|------|-----------------------|-----|----------|----------|----------|----------|
| 5    | **Domain to Species** | 0.764374 | 0.563063 | 0.693150 | 0.563063 | 0.594407 |
| 10   | **Domain to Species** | 0.764374 | 0.612613 | 0.743913 | 0.612613 | 0.646788 |
| 20   | **Domain to Species** | 0.764374 | 0.621622 | 0.752343 | 0.621622 | 0.660354 |
| 50   | **Domain to Species** | 0.764374 | 0.635135 | 0.767735 | 0.635135 | 0.676777 |
| 70   | **Domain to Species** | 0.764374 | 0.635135 | 0.767255 | 0.635135 | 0.674047 |
| 120  | **Domain to Species** | 0.764374 | 0.626126 | 0.763572 | 0.626126 | 0.669832 |
| 200  | **Domain to Species** | 0.764374 | 0.603604 | 0.730455 | 0.603604 | 0.642021 |
| 250  | **Domain to Species** | 0.764374 | 0.639640 | 0.777630 | 0.639640 | 0.682537 |



## Mapping samples bcAd1039, bcAd1046 and bcAd1063 to Paracoccus GCF_034627565.1_ASM3462756v1

| ------------------------------------- | bcAd1039 | bcAd1046 | bcAd1063 |
| ------------------------------------- | -------- | -------- | -------- |
| Mapped Reads                          | 15390    | 39375    | 82361    |
| Alignments                            | 19877    | 42325    | 87874    |
| Mapped Bases                          | 49095177 | 308123371| 733049617|
| Mean Gap-Compressed Sequence Identity | 86.9052% | 96.1005% | 96.5226% |
| Max Mapped Read Length                | 30408    | 38483    | 40843    |
| Mean Mapped Read Length               | 2469.95  | 7279.94  | 8342.05  |








## Wastewater
### Number of contigs and analyzed for methylation
| sample | contigs | .gff files | after filtering |
|--------|---------|------------|-----------------|
| EFF1   | 60537   | 60374      | ?               |
| EFF2   | 59860   | 59697      |                 |
| EFF3   | 54991   | 54838      |                 |
| INF1   | 65422   | 64745      |                 |
| INF2   | 62366   |            |                 |
| INF3   | 62921   |            |                 |
| SLU1   | 74572   | 74388      |                 |
| SLU2   | 76311   |            |                 |
| SLU3   | 73811   |            |                 |

