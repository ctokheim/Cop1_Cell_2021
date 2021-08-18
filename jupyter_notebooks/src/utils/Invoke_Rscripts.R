library(optparse)
library(ggplot2)
library(RColorBrewer)
print(getwd())
source('src/immunesig/utils.R')
#parse input parameters
option_list = list(
  make_option(c("-f_m", "--file_mouse"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f_h", "--file_human"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
	make_option(c("-s", "--species"), type="character", default="mm10", 
              help="input species [default= %default]", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

#check input parameters
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

#read the file in
df = read.table(opt$file, header=TRUE)
if (dim(df)[2] >1){
    df$numeric= which(sapply(df, class)=="numeric")
}
print('input gene  list is: ')
print(head(df))
print(opt$species)

#ICB: Genetech cohort study

#write.table(df_out, file=opt$out, row.names=FALSE)