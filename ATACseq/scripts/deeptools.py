import os
import argparse
from time import gmtime, strftime

def parse_arguments():
    parser=argparse.ArgumentParser(description="This is a description of input args")
    parser.add_argument("-bw","--bigwig",dest="bigwig",help="Bigwig files")
    parser.add_argument("-n","--core",dest="core",default='s',help="Number of cores used to run the program")
    parser.add_argument("-o","--outdir",dest="outdir",default='res',help="Direct of outpit file. Default: res")
    return parser.parse_args()

def main(opts):
    if not os.path.exists(opts.outdir+'/deep_tools/'):
        os.makedirs(opts.outdir+'/deep_tools/')
    outdir=opts.outdir
    bg_file=opts.bigwig
    bg_file=(' ').join([i for i in bg_file.split(',')])
    peak_file='%s/merge_count/merge.sort.unique.peak.bed' % outdir
    n_core=opts.core
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    ##### use deeptools process bigwig data https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html
    ### substract signal of each peak from bigwig data
    os.system('multiBigwigSummary BED-file -b %s -o %s/deep_tools/merge.npz --BED %s --outRawCounts %s/deep_tools/merge_count.txt -p %s ' % (bg_file,outdir,peak_file,outdir,n_core))
    ### Heatmap
    # up down and common
    print('Plot heatmap of all peaks')
    os.system('computeMatrix reference-point  -a 2000 -b 2000 -bs 20 -p 20 --missingDataAsZero -o %s/deep_tools/all_matrix --outFileNameMatrix %s/deep_tools/all_matrix.txt --outFileSortedRegions %s/deep_tools/all_matrix.SortedRegions.bed -R %s/diff_peak/com_peak_center.bed %s/diff_peak/down_peak_center.bed %s/diff_peak/up_peak_center.bed -S %s ' % (outdir,outdir,outdir,outdir,outdir,outdir,bg_file))
    os.system('plotHeatmap -m %s/deep_tools/all_matrix  -out %s/deep_tools/all_heatmap.pdf --plotTitle all --colorList white,red --heatmapHeight 4 --heatmapWidth 4 --refPointLabel Peak --yMin 0 --yMax 2 --boxAroundHeatmaps no' % (outdir,outdir))
    print('Plot heatmap of up peaks')
    os.system('computeMatrix reference-point  -a 2000 -b 2000 -bs 20 -p 20 --missingDataAsZero -o %s/deep_tools/up_matrix --outFileNameMatrix %s/deep_tools/up_matrix.txt --outFileSortedRegions %s/deep_tools/up_matrix.SortedRegions.bed -R %s/diff_peak/up_peak_center.bed -S %s ' % (outdir,outdir,outdir,outdir,bg_file))
    os.system('plotHeatmap -m %s/deep_tools/up_matrix  -out %s/deep_tools/up_heatmap.pdf --plotTitle up --colorList white,red --heatmapHeight 4 --heatmapWidth 4 --refPointLabel Peak --yMin 0 --yMax 2 --boxAroundHeatmaps no' % (outdir,outdir))
    print('Plot heatmap of down peaks')
    os.system('computeMatrix reference-point  -a 2000 -b 2000 -bs 20 -p 20 --missingDataAsZero -o %s/deep_tools/down_matrix --outFileNameMatrix %s/deep_tools/down_matrix.txt --outFileSortedRegions %s/deep_tools/down_matrix.SortedRegions.bed -R %s/diff_peak/down_peak_center.bed -S %s ' % (outdir,outdir,outdir,outdir,bg_file))
    os.system('plotHeatmap -m %s/deep_tools/down_matrix  -out %s/deep_tools/down_heatmap.pdf --plotTitle down --colorList white,red --heatmapHeight 4 --heatmapWidth 4 --refPointLabel Peak --yMin 0 --yMax 2 --boxAroundHeatmaps no' % (outdir,outdir))
    print('Plot heatmap of common peaks')
    os.system('computeMatrix reference-point  -a 2000 -b 2000 -bs 20 -p 20 --missingDataAsZero -o %s/deep_tools/commmon_matrix --outFileNameMatrix %s/deep_tools/commmon_matrix.txt --outFileSortedRegions %s/deep_tools/commmon_matrix.SortedRegions.bed -R %s/diff_peak/com_peak_center.bed -S %s ' % (outdir,outdir,outdir,outdir,bg_file))
    os.system('plotHeatmap -m %s/deep_tools/commmon_matrix  -out %s/deep_tools/common_heatmap.pdf --plotTitle common --colorList white,red --heatmapHeight 10 --heatmapWidth 4 --refPointLabel Peak --yMin 0 --yMax 2 --boxAroundHeatmaps no' % (outdir,outdir))
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

