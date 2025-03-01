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




# Random Forest **All data**:
## All features
| Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|-----------------------|-----|----------|----------|----------|----------|
| **Domain to Species** | 1.0 | 0.870787 | 0.871968 | 0.870787 | 0.859151 |
| **Domain to Genus**   | 1.0 | 0.848315 |          |          |          |
| **Domain to Family**  | 1.0 | 0.893855 |          |          |          |
| **Domain to Order**   | 1.0 | 0.905556 |          |          |          |
| **Domain to Class**   | 1.0 | 0.934066 |          |          |          |
| **Domain to Phylum**  | 1.0 | 1.0      |          |          |          |

## Testing different number of features
| Number of features | Taxa level | Accuracy Train | Accuracy Test | Precision | Recall | F1 Score |
|------|-----------------------|-----|----------|----------|----------|----------|
| 5    | **Domain to Species** | 1.0 | 0.707865 | 0.709828 | 0.707865 | 0.708207 |
| 10   | **Domain to Species** | 1.0 | 0.848315 | 0.851571 | 0.848315 | 0.847413 |
| 20   | **Domain to Species** | 1.0 | 0.870787 | 0.872965 | 0.870787 | 0.867311 |
| 50   | **Domain to Species** | 1.0 | 0.882022 | 0.88702  | 0.882022 | 0.871105 |
| 70   | **Domain to Species** | 1.0 | 0.893258 | 0.893939 | 0.893258 | 0.883861 |
| 120  | **Domain to Species** | 1.0 | 0.882022 | 0.898927 | 0.882022 | 0.867911 |
| 200  | **Domain to Species** | 1.0 | 0.848315 | 0.821896 | 0.848315 | 0.829591 |
| 250  | **Domain to Species** | 1.0 | 0.870787 | 0.871968 | 0.870787 | 0.859151 |

## Random Forest **Filtered data**:

## Random Forest **Chromosome data**:

## Random Forest **Chromosome Filtered data**:














## Wastewater
### Number of contigs and analyzed for methylation
| sample | contigs | .gff files | after filtering |
|--------|---------|------------|-----------------|
| EFF1   | 60537   | 60374      | ?               |
| EFF2   | 59860   |            |                 |
| EFF3   | 54991   |            |                 |
| INF1   | 65422   | 64745      |                 |
| INF2   | 62366   |            |                 |
| INF3   | 62921   |            |                 |
| SLU1   | 74572   | 74388      |                 |
| SLU2   | 76311   |            |                 |
| SLU3   | 73811   |            |                 |