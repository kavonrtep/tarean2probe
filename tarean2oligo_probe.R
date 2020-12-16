#!/usr/bin/env Rscript
suppressPackageStartupMessages({
library(optparse, quiet = TRUE)
})

get_oligo = function(S, ssc=2, formamide=50, optimum=47, min_length=40){
  SS = strsplit(S,"")[[1]]
  LL = seq_along(SS)
  gc=cumsum(SS == "C" | SS=="G")/LL
  T = 81.5+16.6*log10(ssc*0.165)+41*gc-675/LL-0.65*formamide
  T[1:(min_length-1)] = 0
  above_optimum = which(T>optimum)
  if (length(above_optimum) == 0){
    oligo_length = nchar(S)
  }else{
    oligo_length = min(which(T>optimum))
  }
  oligo_tm = T[oligo_length]
  out = substring(S,1,oligo_length)
  names(out) = round(oligo_tm,2)
  return(out)
}


option_list = list(
  make_option(c('-i', '--input_sequences'),action='store',type='character',help='path to fasta file with sequences oriented by TAREAN',default=NULL),
  make_option(c('-g', '--ggmin'),action='store',type='character',help='path to TAREAN file ggmint.RData',default=NULL),
  make_option(c('-d', '--tarean_dir'),action='store',type='character',help='output directory',default=NULL),
  make_option(c('-A', '--use_kmer_file'),action='store_true',type='logical',help='use kmer file instead of ggmin file',default=FALSE),
  make_option(c('-k', '--kmer_length'),action='store',type='numeric',help='what k-mer length to use, only k-mers calculated in ggmin.RData files can be used (typically 11,15,19,23 )',default=15),
  make_option(c('-N', '--number_of_iteration'),action='store',type='numeric',help='Number of iteration - correspond to number of oligos returned',default=NULL),
  make_option(c('-p', '--kmer_prop'),action='store',type='numeric',help='Proportion of most frequent kmers to use, 1 correspond to all kmers',default=1),
  make_option(c('-t', '--oligo_tm'),action='store',type='numeric',help='Required Tm of oligo, default is 47 decgrees of celsia',default=47),
  make_option(c('-o', '--output'), action='store',type='character',help='output base file name',default=NULL)
)

description = paste (strwrap("Script for desing of set of oligo probes from TAREAN results "), collapse ="\n")
epilogue = paste (strwrap("Specify either INPUT_SEQUENCES  and GGMIN file or TAREAN_DIR directory"), collapse ="\n")
parser=OptionParser(
  option_list=option_list,
  epilogue=epilogue,
  description=description
)

opt = parse_args(parser, args=commandArgs(TRUE))

if( !is.null(opt$tarean_dir) & is.null(opt$input_sequences) & is.null(opt$ggmin)){
  s_path = paste0(opt$tarean_dir,"/reads_oriented.fas")
  ggmin_path = paste0(opt$tarean_dir,"/ggmin.RData")
}else if (is.null(opt$tarean_dir) & !is.null(opt$input_sequences) & !is.null(opt$ggmin)){
  s_path = opt$input_sequences
  ggmin_path = opt$ggmin
}else {
  stop("Either TAREAN DIR or GGMIN and INPUT_SEQUENCES must be specified")
}

suppressPackageStartupMessages({
  library(igraph, quiet = TRUE)
  library(Biostrings, quiet = TRUE)
})




s = readDNAStringSet(s_path)
KM = as.character(opt$kmer_length)
## get kmers either from ggmin of from kmer table, -kmer table include all kmers
if (opt$use_kmer_file){
  kmfs = dir(path = opt$tarean_dir, pattern="reads_oriented[.]fas_[0-9]{1,2}[.]kmers")
  KMavail = gsub("[.]kmers$","", gsub("^reads_oriented[.]fas_", "", kmfs))
  if (KM %in% KMavail){
    kmers = read.table(paste0(opt$tarean_dir,"/", kmfs[KMavail %in% KM]), as.is=TRUE, sep="\t")
  }else{
    stop(paste ("kmer length",  KM , "not available\n",
                "possible values: ", KMavail, collapse = " "), "\n")
  }
}else{
  load(ggmin_path)
  if (! KM %in% names(ggmin)){
    stop(paste ("kmer length",  KM , "not available\n",
                "possible values: ", paste(names(ggmin), collapse = " "), "\n")
         )
  }
  kmers = data.frame(seq = V(ggmin[[KM]])$name, freq=V(ggmin[[KM]])$count, stringsAsFactors = FALSE)[order(V(ggmin[[KM]])$count, decreasing = TRUE),]
}
cat("number of available k-mers :", nrow(kmers), "\n")


colnames(kmers) = c("seq", "freq")


Nx_thr= opt$kmer_prop
Niter = 50
kmers_ori = kmers
L = nchar(kmers$seq[1])

schar = sapply(unname(as.character(s)), USE.NAMES=FALSE, get_oligo)



sel_reads_index = numeric()
read_score = numeric()
best_kmers_list = list()
coverage = numeric()
uniqueness = numeric()
counts_list = list()
cat("kmer length            : ", L, "\n")
cat("Number of kmers        : ", nrow(kmers), "\n")
cat("Number of kmers to use : ", min(which(cumsum(kmers$freq)/sum(kmers$freq)>=Nx_thr)), "\n")
cat("Number of reads        : ", length(s), "\n")

Nx = min(which(cumsum(kmers$freq)/sum(kmers$freq)>=Nx_thr))

for (j in 1:Niter){
  counts = integer(length(s))
  Nx_iter = min(Nx, nrow(kmers))
  if (Nx_iter == 0){
    ## no more kmers
    cat("all kmers used\n")
    Niter = j - 1
    break
  }
  cat("Iteration ",j,"\n")
  for (i in 1:Nx_iter){
    positive = grep(kmers$seq[i], schar, fixed=TRUE)
    counts[positive] = counts[positive] + kmers$freq[i]
  }
  ## normalize by oligo length
  counts = counts/nchar(schar)
  counts_list[[j]] = counts
  best_read = s[which.max(counts)]
  W = nchar(best_read)
  best_kmers = substring(as.character(best_read), 1:(W - L +1), L:W)
  km2rm = kmers$seq %in% best_kmers
  ## remove used kmers
  kmers = kmers[!km2rm,]
  sel_reads_index[j] = which.max(counts)
  read_score[j] = max(counts)
  uniqueness[j]  = 100 * (1 - sum(best_kmers %in% unlist(best_kmers_list))/length(best_kmers))
  best_kmers_list[[j]] = best_kmers
  coverage[j] = sum(kmers_ori$freq[kmers_ori$seq %in% best_kmers])/sum(kmers_ori[,2])
}

cumulative_normalized_coverage = sapply(seq_along(best_kmers_list),
                                        function(x){
                                          uk =unique(unlist(best_kmers_list[1:x]))
                                          sum(kmers_ori$freq[kmers_ori$seq %in% uk])/sum(kmers_ori[,2])
                                        }
                                        )




s_out = DNAStringSet(schar[sel_reads_index])


unique_best_kmers = unique(unlist(best_kmers_list))
cat("Number of unique kmers in oligos: ", length(unique_best_kmers),"\n")

kmer_table = matrix(0,nrow=length(unique_best_kmers), ncol=Niter, dimnames = list(unique_best_kmers, names(s)[sel_reads_index]))
for (i in seq_along(best_kmers_list)){
  kmer_table[best_kmers_list[[i]], i] = 1
}

output_df = data.frame(
  read_id = names(s[sel_reads_index]),
  read_sequence = as.character(s[sel_reads_index]),
  oligo_sequence = s_out,
  oligo_tm = names(s_out),
  oligo_length = nchar(s_out),
  read_kmer_uniqueness = uniqueness,
  coverage = coverage,
  cumulative_normalized_coverage = cumulative_normalized_coverage,
  score = read_score
)

write.table(file = paste0(opt$output,"_oligos.csv"), output_df, row.names = FALSE, sep="\t")
write.table(file = paste0(opt$output,"_kmer_ocurence.csv"), kmer_table, row.names = FALSE, sep="\t")
names(s_out) = paste(names(s[sel_reads_index]), names(s_out))
writeXStringSet(s_out, paste0(opt$output,"_oligos.fasta"))
