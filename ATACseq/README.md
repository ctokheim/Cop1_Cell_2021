# Downstream analysis of Cop1 ATACseq data

## Dependencies 

- deeptools 3.5.1
- htseq 0.13.5
- homer 4.11
- DESeq2 1.30.1

## Runing code
Put `bam` files, `bigwig` files, and `peak` files to the corresponding folders:

- `data/bam` for bam files. bam files have been uploaded to GEO.
- `data/bw` for bigwig files. bigwig files have been uploaded to <a href="https://data.mendeley.com/datasets/y4zpyjgjbj/draft?a=36b1f8d9-8fdc-4464-9e87-6fcfa7d4d98f">Mendeley</a>.
- `data/peak` for peak files. peak files can be found in the `data` folder.

Change working directory to `script` directory

``` bash
bash run_pipeline.sh
``` 


