import os
import argparse
from time import gmtime, strftime
import multiprocessing
'''
===================================================================
This tool is used to convert merged bed file to gtf file and count reads of each peak using Htseq-count. This tool is developed by Binbin Wang (wangbinbintj@gmail.com), 
with all rights reserved.
===================================================================
'''
def parse_arguments():
	parser=argparse.ArgumentParser(description="This is a description of input args")
	parser.add_argument("-b","--bed",dest="bed",help="Merged peak bed file")
	parser.add_argument("-ba","--bam",dest="bam",help="bam files")
	#parser.add_argument("-g","--gtf",dest="gtf",default='s',help="gtf files")
	parser.add_argument("-n","--core",dest="core",default='1',help="Number of cores used to run the program")
	parser.add_argument("-o","--outdir",dest="outdir",default='res',help="Direct of outpit file. Default: res")
	return parser.parse_args()
# generate gtf file
def generate_gtf(bed_file,outdir):
	f=open(bed_file,'r')
	tmp_gtf=[]
	i=0
	for line in f:
		if line.strip()!='':
			tmp=line.strip().split('\t')
			tmp_gtf.append(('\t').join([tmp[0],'tmp','exon',tmp[1],tmp[2],'.','.','.','gene_id "%s"'  % ('peak'+str(i))]))
		i+=1
	f.close()
	f=open('%s/count/merge_peak.gtf' % outdir ,'w')
	for j in tmp_gtf:
		f.write(j+'\n')
	f.close()
# count reads using Htseq
def count_reads(bam_file):
	prefix=('.').join(bam_file.split('/')[-1].split('.')[0:-1])
	print('htseq-count --stranded=no -f bam %s %s/count/merge_peak.gtf > %s/count/%s.count.txt' % (bam_file,outdir,outdir,prefix))
	os.system('htseq-count --stranded=no -f bam %s %s/count/merge_peak.gtf > %s/count/%s.count.txt' % (bam_file,outdir,outdir,prefix))

def main(opts):
	if not os.path.exists(opts.outdir+'/count/'):
	    os.makedirs(opts.outdir+'/count/')
	generate_gtf(opts.bed,opts.outdir)
	bam_files=opts.bam
	#bam_files=[i for i in bam_files.split(',')]
	for i in bam_files.split(','):
		count_reads(i)
	# count reads with multi processes
	#gtf_file='%s/count/merge_peak.gtf' % opts.outdir
	#pool=multiprocessing.Pool(processes=int(opts.core))
	#result=pool.map_async(count_reads,bam_files)
	#pool.close()
	#pool.join()

if __name__ == '__main__':
    opts = parse_arguments()
    outdir=opts.outdir
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    main(opts)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))





