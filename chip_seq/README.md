# ChIP-seq Analysis

This code is meant to analyze the Cebpd ChIP-seq dataset.

## Installation

You will need to install two environments. One is to run the "chips" pipeline and the second is to subsequently analyze the results.

To install the environment for basic ChIP-seq processing, please use the following conda command:

```bash
$ conda env create -f cidc_chips/environment.yml 
```

To install the environment to subsequently analyze the ChIP-seq peak calls please use the following conda command:

```bash
$ conda env create -f environment.yml
```

## Running code

To run the ChIP-seq pipeline all the way to peak calling, first download the Cebpd ChIP-seq FASTQ files and correctly note their downloaded path in the configuration file ("config.yaml"). Then issue the following commands:

```bash
$ conda activate cidc_chips
$ snakemake -s cidc_chips/chips.snakefile -p -j 16 --latency-wait 100 --notemp
```

Where the -j parameter specifies the number of cores that you want to use.

After the pipeline run finishes, issue the following commands:

```bash
$ conda activate post_analysis
$ ./peak_analysis.sh analysis
```

You should find results in the "analysis" folder, with further post-processing results in the "analysis/union" folder.
