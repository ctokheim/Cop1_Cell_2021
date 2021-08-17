for f in files/deseq2/Condition/*.txt ; do
    mybase=`basename $f .txt`
    mkdir -p files/gsea/$mybase
    Rscript scripts/mygsea.R -m $f -k data/h.all.v7.4.entrez.gmt -o files/gsea/$mybase
done
