# 16-ass-gastrit
Repository for code created at the ITMO bioinformatics hackathon. 16S-analysis project. "My bacteria save me"
<p align="center">
  <img width="500" alt="image" src="https://github.com/ailiskab-hub/16-ass-gastrit/blob/main/pics/photo_2024-06-25_16-12-01.jpg">
</p> 


### Data
Almost all the data are avilable in the folder ``` /data ```. Data of the health person obtainen from *American Gut Project* and it can be dowloaded [there](https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00000596/pipelines/5.0/file/ERP012803_taxonomy_abundances_SSU_v5.0.tsv)

### Raw data processing
The raw ```fastq``` reads were processed using the DADA2 pipeline.

Script can be launched using following command
```r
Rscript Dada2_scr.R -d Central   # for the central samples
```

### Subsequent analysis
Our hypotheses
- H0: "There is no difference between samples obtained from different regions"
- H1: 'There is some difference between samples obtained in different regions'

### Results
Based on our analysis, we accepted H0. Thus, there is no significant difference between microbiota in gastritis in different regions

More details on the analysis methods, results and conclusions of our study will be available at the [link](https://docs.google.com/document/d/1bZTrj9Nc1Jj4KxGpigHii50RLc0wTp1ZZgyn0IoR1_4/delite_this_part) that will be updated later
