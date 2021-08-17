output_dir=$1
mkdir -p $output_dir/union/5fold
mkdir -p $output_dir/union/5fold/individual
# compile links for all bed files
for d in `pwd`/$output_dir/peaks/* ; do
    mybase=`basename $d`
    if [ -d $d ]
    then
    ln -s $d/`echo $mybase`_5fold_peaks.bed $output_dir/union/5fold/individual/`echo $mybase`_5fold_peaks.bed
    fi
done
# concatenate bed files
mkdir -p $output_dir/union/concat_bed
cat $output_dir/union/5fold/individual/*.bed | sort -k1,1 -k2,2n > $output_dir/union/concat_bed/concat_v2.bed 
# take union of bed regions
bedtools merge -i $output_dir/union/concat_bed/concat_v2.bed > $output_dir/union/union_v2.bed
# compile all bam links
mkdir -p $output_dir/union/bam/
mkdir -p $output_dir/union/counts/
for d in `pwd`/$output_dir/align/* ; do
    mybase=`basename $d`
    if [[ -d $d ]]; then
        fname=`echo $mybase`_unique.sorted.bam
        ln -s $d/$fname $output_dir/union/bam/$fname
    fi
done
# calc read counts for peaks
for f in $output_dir/union/bam/*.bam ; do
    mybase=`basename $f _unique.sorted.dedup.bam`
    bedtools coverage -a $output_dir/union/union_v2.bed -b $f -counts > $output_dir/union/counts/`echo $mybase`_v2.txt
done
# merge read counts
python scripts/merge_counts.py -i $output_dir/union/counts/ -o $output_dir/union/peak_union_counts_v2.txt
# perform differential peak analysis
Rscript scripts/diff_peak_all.R -i $output_dir/union/peak_union_counts_v2.txt -o $output_dir/union/diff_peak_rosa.txt
# split peaks into different bins
mkdir -p $output_dir/union/tmp_diff_bed
python scripts/split_diff_bed.py -i $output_dir/union/diff_peak_rosa.txt -o $output_dir/union/tmp_diff_bed
# link to bigwigs
mkdir -p $output_dir/union/bw
for d in `pwd`/$output_dir/peaks/* ; do
    mybase=`basename $d`
    if [[ -d $d ]]; then
        fname=`echo $mybase`_treat_pileup.bw
        ln -s $d/$fname $output_dir/union/bw/$fname
    fi
done
# create mered bigwig files
mkdir -p $output_dir/union/merged_bw
# cop1 sgrna1
bigWigMerge $output_dir/peaks/COP1_KO_sgRNA1.rep1/COP1_KO_sgRNA1.rep1_treat_pileup.bw $output_dir/peaks/COP1_KO_sgRNA1.rep2/COP1_KO_sgRNA1.rep2_treat_pileup.bw $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup.bdg
cat $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup.bdg | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4/2}' > $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup_average.bdg
wigToBigWig $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup_average.bdg data/mm10.chrom.sizes $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup_average.bw 
# cop1 sgrna2
bigWigMerge $output_dir/peaks/COP1_KO_sgRNA2.rep1/COP1_KO_sgRNA2.rep1_treat_pileup.bw $output_dir/peaks/COP1_KO_sgRNA2.rep2/COP1_KO_sgRNA2.rep2_treat_pileup.bw $output_dir/union/merged_bw/COP1_KO_sgRNA2.treat_pileup.bdg
cat $output_dir/union/merged_bw/COP1_KO_sgRNA2.treat_pileup.bdg | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4/2}' > $output_dir/union/merged_bw/COP1_KO_sgRNA2.treat_pileup_average.bdg
wigToBigWig $output_dir/union/merged_bw/COP1_KO_sgRNA2.treat_pileup_average.bdg data/mm10.chrom.sizes $output_dir/union/merged_bw/COP1_KO_sgRNA2.treat_pileup_average.bw 
# all cop1
bigWigMerge $output_dir/peaks/COP1_KO_sgRNA1.rep1/COP1_KO_sgRNA1.rep1_treat_pileup.bw $output_dir/peaks/COP1_KO_sgRNA1.rep2/COP1_KO_sgRNA1.rep2_treat_pileup.bw $output_dir/peaks/COP1_KO_sgRNA2.rep1/COP1_KO_sgRNA2.rep1_treat_pileup.bw $output_dir/peaks/COP1_KO_sgRNA2.rep2/COP1_KO_sgRNA2.rep2_treat_pileup.bw $output_dir/union/merged_bw/COP1_KO.treat_pileup.bdg
cat $output_dir/union/merged_bw/COP1_KO.treat_pileup.bdg | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4/4}' > $output_dir/union/merged_bw/COP1_KO.treat_pileup_average.bdg
wigToBigWig $output_dir/union/merged_bw/COP1_KO.treat_pileup_average.bdg data/mm10.chrom.sizes $output_dir/union/merged_bw/COP1_KO.treat_pileup_average.bw 
# WT2
bigWigMerge $output_dir/peaks/WT2.rep1/WT2.rep1_treat_pileup.bw $output_dir/peaks/WT2.rep2/WT2.rep2_treat_pileup.bw $output_dir/union/merged_bw/WT2.treat_pileup.bdg
cat $output_dir/union/merged_bw/WT2.treat_pileup.bdg | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4/2}' > $output_dir/union/merged_bw/WT2.treat_pileup_average.bdg
wigToBigWig $output_dir/union/merged_bw/WT2.treat_pileup_average.bdg data/mm10.chrom.sizes $output_dir/union/merged_bw/WT2.treat_pileup_average.bw 
# ROSA
bigWigMerge $output_dir/peaks/ROSA.rep1/ROSA.rep1_treat_pileup.bw $output_dir/peaks/ROSA.rep2/ROSA.rep2_treat_pileup.bw $output_dir/union/merged_bw/ROSA.treat_pileup.bdg
cat $output_dir/union/merged_bw/ROSA.treat_pileup.bdg | awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4/2}' > $output_dir/union/merged_bw/ROSA.treat_pileup_average.bdg
wigToBigWig $output_dir/union/merged_bw/ROSA.treat_pileup_average.bdg data/mm10.chrom.sizes $output_dir/union/merged_bw/ROSA.treat_pileup_average.bw 

# creat heatmap of peaks split by up/down
computeMatrix reference-point -S $output_dir/union/merged_bw/WT2.treat_pileup_average.bw $output_dir/union/merged_bw/ROSA.treat_pileup_average.bw $output_dir/union/merged_bw/COP1_KO_sgRNA1.treat_pileup_average.bw -R $output_dir/union/tmp_diff_bed/up_cop1_ko_vs_rosa.bed $output_dir/union/tmp_diff_bed/down_cop1_ko_vs_rosa.bed -a 1000 -b 1000 -p 8 --referencePoint center --missingDataAsZero --samplesLabel WT ROSA COP1_sgRNA1 -o $output_dir/union/mat_diff_group_up_down.txt.gz
plotHeatmap -m $output_dir/union/mat_diff_group_up_down.txt.gz --colorMap Blues Blues Oranges --yMax 4 4 4 --heatmapHeight 12 -out $output_dir/union/Heatmap_up_down.pdf
