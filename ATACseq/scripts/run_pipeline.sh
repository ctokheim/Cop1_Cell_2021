### cop1_ifn_vs_rosa26_ifn

# step1: merge peaks 
python merge_peak.py -n 4 -o ../results/cop1_ifn_vs_rosa26_ifn -p ../data/peak/Cop1sg1_IFN_r1.rep1_sorted_peaks.narrowPeak,../data/peak/Cop1sg1_IFN_r2.rep1_sorted_peaks.narrowPeak,../data/peak/Cop1sg2_IFN_r1.rep1_sorted_peaks.narrowPeak,../data/peak/Cop1sg2_IFN_r2.rep1_sorted_peaks.narrowPeak,../data/peak/Rosa26_IFN_r1.rep1_sorted_peaks.narrowPeak,../data/peak/Rosa26_IFN_r2.rep1_sorted_peaks.narrowPeak

# step2: tranlate bed file to gtf file and count reads using htseq-count
python count_reads.py -b ../results/cop1_ifn_vs_rosa26_ifn/merge_count/merge.sort.unique.peak.bed \
-ba ../data/bam/Cop1sg1_IFN_r1_unique.sorted.bam,../data/bam/Cop1sg1_IFN_r2_unique.sorted.bam,../data/bam/Cop1sg2_IFN_r1_unique.sorted.bam,../data/bam/Cop1sg2_IFN_r2_unique.sorted.bam,../data/bam/Rosa26_IFN_r1_unique.sorted.bam,../data/bam/Rosa26_IFN_r2_unique.sorted.bam \
-n 4 \
-o ../results/cop1_ifn_vs_rosa26_ifn

# step 3: identify differential peaks, prefix ignore .count.txt
/usr/local/bin/Rscript diff_peak_DEseq.R -c Rosa26_IFN_r1_unique.sorted,Rosa26_IFN_r2_unique.sorted \
-t Cop1sg1_IFN_r1_unique.sorted,Cop1sg1_IFN_r2_unique.sorted,Cop1sg2_IFN_r1_unique.sorted,Cop1sg2_IFN_r2_unique.sorted -o ../results/cop1_ifn_vs_rosa26_ifn -p 0.05 -f 1 


# setp 4: deeptools: correlation, heatmap
python deeptools.py -n 4 -o ../results/cop1_ifn_vs_rosa26_ifn -bw ../data/bw/Cop1sg1_IFN_r1.rep1_treat_pileup.bw,../data/bw/Cop1sg1_IFN_r2.rep1_treat_pileup.bw,../data/bw/Cop1sg2_IFN_r1.rep1_treat_pileup.bw,../data/bw/Cop1sg2_IFN_r2.rep1_treat_pileup.bw,../data/bw/Rosa26_IFN_r1.rep1_treat_pileup.bw,../data/bw/Rosa26_IFN_r2.rep1_treat_pileup.bw

# step 5.1: homer motif
# install mm9 genome for homer
perl /Users/ziyili/Documents/binbin/tools/homer/.//configureHomer.pl -install mm9

python homer_motif.py -o ../results/cop1_ifn_vs_rosa26_ifn -r mm9 

# step 5.2: scatter plot of enriched motif
Rscript mortif_scatter.R -f ../results/cop1_ifn_vs_rosa26_ifn/homer/up_peak/knownResults.txt -p up_peak -o ../results/cop1_ifn_vs_rosa26_ifn
Rscript mortif_scatter.R -f ../results/cop1_ifn_vs_rosa26_ifn/homer/down_peak/knownResults.txt -p down_peak -o ../results/cop1_ifn_vs_rosa26_ifn

### cop1_vs_rosa26
# step1: merge peaks 
python merge_peak.py -n 4 -o ../results/cop1_vs_rosa26 -p \
../data/peak/Cop1sg1_r1.rep1_sorted_peaks.narrowPeak,\
../data/peak/Cop1sg1_r2.rep1_sorted_peaks.narrowPeak,\
../data/peak/Cop1sg2_r1.rep1_sorted_peaks.narrowPeak,\
../data/peak/Cop1sg2_r2.rep1_sorted_peaks.narrowPeak,\
../data/peak/Rosa26_r1.rep1_sorted_peaks.narrowPeak,\
../data/peak/Rosa26_r2.rep1_sorted_peaks.narrowPeak

# step2: translate bed file to gtf file and count reads using htseq-count
python count_reads.py -b \
../results/cop1_vs_rosa26/merge_count/merge.sort.unique.peak.bed \
-ba ../data/bam/Cop1sg1_r1_unique.sorted.bam,\
../data/Cop1sg1_r2_unique.sorted.bam,\
../data/Cop1sg2_r1_unique.sorted.bam,\
../data/Cop1sg2_r2_unique.sorted.bam,\
../data/Rosa26_r1_unique.sorted.bam,\
../data/Rosa26_r2_unique.sorted.bam \
-n 4 \
-o ../results/cop1_vs_rosa26

# step 3: identify differential peaks, prefix ignore .count.txt
Rscript diff_peak_DEseq.R -c Rosa26_r1_unique.sorted,Rosa26_r2_unique.sorted \
-t Cop1sg1_r1_unique.sorted,Cop1sg1_r2_unique.sorted,Cop1sg2_r1_unique.sorted,Cop1sg2_r2_unique.sorted -o ../results/cop1_vs_rosa26 -p 0.05 -f 1 

# setp 4: deeptools: correlation, heatmap
python deeptools.py -n 4 -o ../results/cop1_vs_rosa26 -bw ../data/bw/Cop1sg1_r1.rep1_treat_pileup.bw,\
../data/bw/Cop1sg1_r2.rep1_treat_pileup.bw,\
../data/bw/Cop1sg2_r1.rep1_treat_pileup.bw,\
../data/bw/Cop1sg2_r2.rep1_treat_pileup.bw,\
../data/bw/Rosa26_r1.rep1_treat_pileup.bw,\
../data/bw/Rosa26_r2.rep1_treat_pileup.bw

# step 5.1: homer motif
python scripts/homer_motif.py -o ../results/cop1_vs_rosa26 -r mm9 

# step 5.2: plot enriched motif
Rscript mortif_scatter.R -f ../results/cop1_vs_rosa26/homer/up_peak/knownResults.txt -p up_peak -o ../results/cop1_vs_rosa26
Rscript mortif_scatter.R -f ../results/cop1_vs_rosa26/homer/down_peak/knownResults.txt -p down_peak -o ../results/cop1_vs_rosa26