import os
import argparse
from time import gmtime, strftime
import multiprocessing
'''
===================================================================
This tool is used to merge all peaks and count reads of each peak. This tool is developed by Binbin Wang (wangbinbintj@gmail.com), 
with all rights reserved.
===================================================================
'''
def parse_arguments():
	parser=argparse.ArgumentParser(description="This is a description of input args")
	parser.add_argument("-p","--peaks",dest="peaks",help="Paek files of all samples")
	#parser.add_argument("-bw","--bigwig",dest="bigwig",default='s',help="bigwig files")
	parser.add_argument("-n","--core",dest="core",default='s',help="Number of cores used to run the program")
	parser.add_argument("-o","--outdir",dest="outdir",default='res',help="Direct of outpit file. Default: res")
	return parser.parse_args()

def main(opts):
	print('Step 1: merge peaks and calculate read count \n')
	peak_file=opts.peaks
	n_core=opts.core
	outdir=opts.outdir

	peak_file=[i for i in peak_file.split(',')]
	peak_file=(' ').join(peak_file)

	if not os.path.exists(outdir+'/merge_count/'):
	    os.makedirs(outdir+'/merge_count/')
	##### merge peak
	# record running time
	print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
	# cat peaks
	os.system('cat %s > %s/merge_count/merge.peak.bed' % (peak_file,outdir))

	#sort the peak file first:
	os.system('sort -k1,1 -k2,2n %s/merge_count/merge.peak.bed > %s/merge_count/merge.sort.peak.bed' % (outdir,outdir))

	# merge peaks http://bedtools.readthedocs.io/en/latest/content/tools/merge.html
	os.system('bedtools merge -i %s/merge_count/merge.sort.peak.bed > %s/merge_count/merge.sort.unique.peak.bed' % (outdir,outdir))
	print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
