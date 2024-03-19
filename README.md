# miRTSS workflow
See [notebook](https://nbviewer.org/github/tzhu-bio/microRice/blob/main/miRTSS/miRTSS_training.ipynb) for CNN model training.

After the model has been trained, follow the steps below to get the specific TSS.

## get_nearest_candidate_TSS.py

### Background

We trained four different models using H3K4me3 data from different tissues, Each model scores candidate genome regions, providing a likelihood of being a transcription start site (TSS). This script aggregates all scoring results, allowing for the option to either select all of them or a subset based on biological situation. Then, a custome cutoff (defaulting to 0.6) is applied, where regions scoring above this threshold are identified as TSS candidate sites. Next, we select the closest candidate TSS to each miRNA,  temporarily considering it as the miRNA's potential TSS. 

---

### Example command

```python
python get_nearest_candidate_TSS.py -a miRNA.bed -b *H3K4me3_predicted_score_predict_miRNA* -cutoff 0.6
```

---

### Input parameters

- **`-a`** or **`--miRNAfile`**:
    
     .bed file containing genome location of pre-miRNAs.
    
- **`-b`** or **`--resultfiles`**:
    
    TSS predicted score result outputs from our models. 
    
- **`-cutoff`** or **`--score_cutoff`**:
    
    Threshold score for TSS identification. Regions with scores exceeding this cutoff will be considered as TSS, with a default value of 0.6.
    

### Output

**`tss_candidate.csv`** contains details of candidate TSS for miRNAs. Each row corresponds to a miRNA and includes the following columns:

1. **`Chromosome`**: Chromosomal location of the miRNA.
2. **`Start`**: Starting position of the miRNA.
3. **`End`**: Ending position of the miRNA.
4. **`Name`**: Name of the miRNA.
5. **`Start_b`**: Starting position of the candidate TSS.
6. **`End_b`**: Ending position of the candidate TSS.
7. **`Strand_b`**: Strand direction of the miRNA.
8. **`Distance`**: Distance from the candidate TSS to the miRNA.
9. **`Score`**: TSS score assigned to the tss region

## miRNA_classify.py

### Background

Intragenic miRNAs and their host genes either share the promoter or have independent transcription starting site [(Ref)](https://pubmed.ncbi.nlm.nih.gov/30785618/). It is widely acknowledged that intragenic miRNAs are generally transcribed along with their host genes. Therefore, if there are no candidate TSSs within 5 kb upstream intragenic miRNA, the region within 1 kb upstream of the TSS of the host gene is considered as the miRNA TSS. This script aims to identify intragenic miRNAs (excluding antisense miRNAs, as they have opposite transcriptional orientations to their host genes and are unlikely to share promoters) along with their respective host gene genomic locations.

![Untitled](https://github.com/tzhu-bio/microRice/blob/main/png/1.png)

---

### Example command

```python
python miRNA_classify.py -a miRNA.bed -b gene.bed
```

---

### Input parameters

- **`-a`** or **`--miRNAfile`**:
    
     .bed file containing genome locations of pre-miRNAs. Make sure that this file includes strand information and adheres to the BED format specifications.
    
- **`-b`** or  **`--genefile`**:
    
    .bed file containing genome locations of protein-coding genes, long non-coding RNAs, or both, depending on biological situation. Make sure that this file includes strand information and adheres to the BED format specifications. 
    

### Output

**`miRNA_classification.txt`**  contains miRNAs classification. Each row corresponds to a miRNA and includes the following columns:

1. Name of the miRNA.
2. classification of the miRNA.

**`intragenic_miRNA_samestrand_hostgene.csv`** contains details of intragenic miRNAs located on the same strand as their host genes. Each row represents a miRNA and includes the following columns:

1. **`Chromosome`**: Chromosomal location of the miRNA.
2. **`Start`**: Starting position of the miRNA.
3. **`End`**: Ending position of the miRNA.
4. **`Name`**: Name of the miRNA.
5. **`Strand`** : Strand direction of the miRNA.
6. **`Start_b`**: Starting position of the host gene.
7. **`End_b`**: Ending position of the host gene.
8. **`Strand_b`**: Strand direction of the host gene.
9. **`Name_b`**: Name of the host gene.

## modify_candidate_TSS.py

### Background

This script aims to modify the output results obtained from **`get_nearest_candidate_TSS.py`**. It will assign host gene TSSs to intragenic miRNAs located on the same strand as their host genes and lacking candidate TSSs within 5kb upstream.  This distance can be customized.  Before running this script, make sure to run **`get_nearest_candidate_TSS.py`** and **`miRNA_classify.py`** first.

---

### Example command

```python
python modify_candidate_TSS.py --distance_threshold 5000
```

---

### Input parameters

- **`--distance_threshold`**:
    
    The maximum distance from intragenic miRNAs (excluding antisense miRNAs) to candidate TSSs. Intragenic miRNAs with candidate TSSs exceeding this distance will be assigned the TSS of their host gene.
    

### Output

**`miRNA_TSS.bed`** contains details of candidate TSS for miRNAs after modifying intragenic miRNAs. Each row corresponds to a miRNA and includes the following columns:

1. Chromosomal location of the predicted miRNA TSS. 
2. Starting position of the predicted miRNA TSS. 
3. Ending position of the predicted miRNA TSS. 
4. Name of the miRNA.
5. place holder column according to .bed file format
6. Strand direction of the miRNA. 
7. Starting position of the miRNA.
8. Ending position of the miRNA.
9. Distance from the candidate TSS to the miRNA.
10. TSS score of this predicted TSS region. 
11. TSS category."H3K4me3_enriched_zone" indicates that the predicted TSS is the closest candidate TSS to this miRNA, which is also a H3K4me3-enriched region. "hostgene_TSS" indicates that the predicted TSS is host gene TSS of this intragenic miRNA.
