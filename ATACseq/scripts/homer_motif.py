import os
import argparse
from time import gmtime, strftime

def parse_arguments():
    parser=argparse.ArgumentParser(description="This is a description of input args")
    #parser.add_argument("-hp","--homer_path",dest="homer_path",help="path of homer")
    #parser.add_argument("-hl","--homer_lib_path",dest="homer_lib_path",help="lib path of homer")
    parser.add_argument("-r","--ref_genome",dest="ref_genome",help="refrence genome")
    parser.add_argument("-o","--outdir",dest="outdir",default='res',help="Direct of outpit file. Default: res")
    return parser.parse_args()

def main(opts):
    if not os.path.exists(opts.outdir+'/homer/'):
        os.makedirs(opts.outdir+'/homer/')
    if not os.path.exists(opts.outdir+'/homer/up_peak/'):
        os.makedirs(opts.outdir+'/homer/up_peak/')
    if not os.path.exists(opts.outdir+'/homer/down_peak/'):
        os.makedirs(opts.outdir+'/homer/down_peak/')

    outdir=opts.outdir
    ref_genome=opts.ref_genome
    ### export homer and lib path
    #os.system('export PATH="%s"' % opts.homer_path)
    #os.system('export PATH="%s"' % opts.homer_lib_path)
    #os.environ['PATH'] = opts.homer_path
    #os.environ['PATH'] = opts.homer_lib_path
    ### run homer with up and down peaks
    print('Motif finding using homer of up peaks')
    os.system('findMotifsGenome.pl %s/diff_peak/up_peak.bed %s %s/homer/up_peak/ -preparse' % (outdir,ref_genome,outdir))

    print('Motif finding using homer of down peaks')
    os.system('findMotifsGenome.pl %s/diff_peak/down_peak.bed %s %s/homer/down_peak/ -preparse' % (outdir,ref_genome,outdir))
    ### scatter plot of enriched 
    #os.system('Rscript motif_scatter.R -f %s/homer/up_peak/knownResults.txt -p up_peak -o %s' % (outdir,outdir))
    #os.system('Rscript motif_scatter.R -f %s/homer/down_peak/knownResults.txt -p down_peak -o %s' % (outdir,outdir))


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

