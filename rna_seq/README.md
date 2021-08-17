# Cop1 and Cebpd KO RNA-seq analysis

The RIMA pipeline is being used to analyze the RNA-seq data. 

```bash
$ cd RIMA_Kraken
$ sbatch run.sh
```

## Samples and meta-information

The RNA-seq data was done in two batches, one for the Cop1 KO experiments and a second for the Cebpd KO experiments. There are generally two replicates per condition. Although because the Rosa26 conditon (control guide RNA) was repeated in both batches, there are four replicates for Rosa26. In addition, conditions may be with and without IFNG treatment. The sample information is recorded in the metasheet.txt file in the RIMA_Kraken directory. All fastq files can be found in the raw_data directory, and the sample conditions for each files is recorded in the config.yaml file in the RIMA_Kraken directory.


## RIMA pipeline

The RIMA pipeline was cloned from github on 5/23/2021. However, as long as you clone the RIMA_kraken branch you should be fine.

```bash
$ git clone https://github.com/liulab-dfci/RIMA_pipeline
```

I filled out the metasheet.txt and config.yaml according to conversations with Stan and Xiaoqing. This included figuring out the RNA library type used for this analysis is paired-end unstranded ("fr-unstranded"). You will need to copy the files in this directory into the RIMA folder cloned from github.

I first ran the level 1 jobs by marking "true" in the execution.yaml file, which mostly consists of read mapping and QC. I then ran level2 by submitting a second time, but with level 1 marked as "false" and level 2 as "true". Level 2 consists mostly of differential expression, batch effect removal and gene set enricment.

## Additional notes

I performed the DEG analysis for the Cop1 data using only data from the original batch. For the cebpd analysis, I included all samples, and added the batch as a covariate. However, only the results for the Cebpd KO vs Rosa are useful.

I ultimately had to run most of the level 2 analysis in separate steps. I used it before DEG analysis using Dseq2, but had to create my own GSEA script ("scripts/mygsea.R") so as to avoid using the raw log2FC for GSEA analysis. Lastly, I ran the batch removal step on all the samples in a separate run of the RIMA pipeline. For using run.sh, don't forget to edit the metasheet to contain the samples that you want to run.

```bash
$ # Edit execution.yaml to only run deseq2
$ sbatch run.sh
$ # Use degrader env on noah
$ ./run_gsea.sh
$ # edit execution.yaml to run only batch removal for all samps
$ sbatch run.sh
```
