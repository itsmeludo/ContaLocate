#!/usr/bin/Rscript

### Author: Ludovic V. Mallet, PhD
### 2016.03.22
### licence: GPLv3
### Version: Alpha.1
### Garanteed with misatkes. <- Including this one.
### Important notes:


### dependencies:
library(gtools)
library(getopt) #install.packages("getopt")   #maybe launch it as admin or root if you want it for other users too.
# Kount.py in the exec PATH or where this script is launched



### Standard parameters
win_step=500 #bp: sliding step
win_size=5000 #bp:window size for composition computation
dist="KL"
# 
# genome_fasta="/home/ludovic/taff/projects/inra/sequences/TH12_prefix.fna"
# host_sample_fasta="/home/ludovic/taff/projects/inra/sequences/TH12_prefix.fna"
# conta_sample_fasta="TH12_positions_Burk.GFF.gdna.fasta"

# spec <- matrix(c(
#         'genome'         , 'i', 1, "character", "file from fastq-stats -x (required)",
#         'host_learn'     , 'r', 2, "character", "input gc content file (optional)",
#         'conta_learn'    , 'c', 1, "character", "output filename (optional)",
#         'win_step'       , 't', 2, "int", "output filename (optional)",
#         'win_size'       , 'w', 2, "int", "output filename (optional)",
#         'dist'           , 'd', 2, "character", "KL, JSD or Eucl. the divergence function used. respectivelly Kullback-Leibler, Jensen_Shannon and Euclidean. Default: KL(optional)",
#         'help'           , 'h', 0, "logical",   "this help"
# ),ncol=5,byrow=T)

spec <- matrix(c(
        'genome'         , 'i', 1, "character", "file from fastq-stats -x (required)",
        'host_learn'     , 'r', 2, "character", "input gc content file (optional)",
        'conta_learn'    , 'c', 1, "character", "output filename (optional)",
        'win_step'       , 't', 2, "int", "output filename (optional)",
        'win_size'       , 'w', 2, "int", "output filename (optional)",
        'dist'           , 'd', 2, "character", "KL, JSD or Eucl"
        
),ncol=5,byrow=T)


opt = getopt(spec);
# [[""]]
if (!is.null(opt[["help"]]) || is.null(opt[["genome"]])) {
    cat(paste(getopt(spec, usage=T),"\n"));
}
cat(paste(getopt(spec, usage=T),"\n"));
genome_fasta = opt[["genome"]]
conta_sample_fasta = opt[["conta"]]
host_sample_fasta = ifelse(is.null(opt[["host_learn"]]), test <- "" , test <-opt[["host"]])
dist = ifelse(is.null(opt[["dist"]]), test <- dist , test <-opt[["dist"]])
win_step = ifelse(is.null(opt[["win_step"]]), test <- win_step , test <-opt[["win_step"]])
win_size = ifelse(is.null(opt[["win_size"]]), test <- win_size , test <-opt[["win_size"]])


### compute profiles:
data=list()
if ( is.null(opt[["host_learn"]])) {
  system(paste("ionice -c2 -n3 Kount.py -i ",genome_fasta ," -w ",win_size," -t ",win_step," -d  ",dist," -n 0.5",sep=""))
  data[["host"]]=read.delim(file=paste(basename(genome_fasta),".mcp_windows_vs_whole_",dist,".dist",sep=""), header=F)
}else{
  system(paste("ionice -c2 -n3 Kount.py -i ",genome_fasta ," -w ",win_size," -t ",win_step," -c ",conta_sample_fasta," -r ",host_sample_fasta," -d ",dist," -n 0.5",sep=""))
  data[["host"]]=read.delim(file=paste(basename(genome_fasta),".mcp_hostwindows_vs_host_",basename(host_sample_fasta),"_",dist,".dist",sep=""), header=F)
}
data[["conta"]]=read.delim(file=paste(basename(genome_fasta),".mcp_hostwindows_vs_conta_",basename(conta_sample_fasta),"_",dist,".dist",sep=""), header=F)


### Ask the trusty human to set the thresholds:

threshold_conta=0
repeat{
  plot(density(data[["conta"]][,4],na.rm=TRUE),xlim=c(0,5000),lwd=2)
  abline(v=threshold_conta,col="red")
  new_threshold= ask("Give a different threshold value for the contaminant threshold. Give the same value to confirm it.\n")
   new_threshold <- as.numeric(new_threshold)
   
  if(new_threshold==threshold_conta){
    break
  }
  threshold_conta=new_threshold
}


threshold_host=0
repeat{
  plot(density(data[["host"]][,4],na.rm=TRUE),xlim=c(0,5000),lwd=2)
  abline(v=threshold_host,col="red")
  new_threshold= ask("Give a different threshold value for the host threshold. Give the same value to confirm it.\n")
   new_threshold <- as.numeric(new_threshold)
   
  if(new_threshold==threshold_host){
    break
  }
  threshold_host=new_threshold
}



### Perform the split over the double threshold

data[["Select_conta"]]=which((data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)>=1)
data[["windows_conta"]]=data[["conta"]][data[["Select_conta"]],]
# write(file=paste(f,"fenetres_Burk.txt",sep=""),as.character(data[["Fenetres_Burk"]]))


# a=cbind((data[["conta"]][,4]<= threshold_conta ),(data[["host"]][,4]>= threshold_host),((data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)),(data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)>=1)
# b=data[["conta"]][data[["Select_conta"]],]



### write a GFF file of the positions of the targeted species
# regroup_struct=data[["Select_conta"]]
write(x="##gff-version 2", file = paste(basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
start_index=data[["Select_conta"]][1]
for(i in seq(length(data[["Select_conta"]]))){
#   if(ifelse(is.na(data[["Select_conta"]][i+1]), test <- 0 , test <- data[["Select_conta"]][i+1]) == data[["Select_conta"]][i]+1 | data[["conta"]][i,1] == ifelse(is.na(data[["conta"]][i+1,1]), test <- -1 , test <- data[["conta"]][i+1,1])){ #2 windows have a consecutive index number, so they are physically touching. data[["conta"]][i,1] describe the contig, and regions are only spanning one contig max. we assert that ssuccessive indexes are in the same contig to group them. If everything of this if fine, we will then group the windows in a region, by etending the end of it.
  if(ifelse(is.na(data[["Select_conta"]][i+1]), test <- 0 , test <- data[["Select_conta"]][i+1]) == data[["Select_conta"]][i]+1 ){ #2 windows have a consecutive index number, so they are physically touching. data[["conta"]][i,1] describe the contig, and regions are only spanning one contig max. we assert that ssuccessive indexes are in the same contig to group them. If everything of this if fine, we will then group the windows in a region, by etending the end of it.
    end_index=data[["Select_conta"]][i+1]
  }else{ #well, the window index i+1 is not touching the window i, so i is the last of its island, and the island can be written. i+1 is the start of a new island.
    end_index=data[["Select_conta"]][i]
    line=paste(sep="\t", as.character(data[["conta"]][start_index,1]),"SignatureGohtam\tregion",data[["conta"]][start_index,2],data[["conta"]][end_index,c(3)])#,c(paste(sep="/", data[["conta"]][start_index:end_index,4])))
#     print(line)
    write(x=line,append=TRUE, file = paste(basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
    start_index=data[["Select_conta"]][i+1]
  }
}


### Booom Done

# image.save(file="test.Rdata")





