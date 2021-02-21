library(Biostrings)
library(data.table)
genomes=readBStringSet("D:/Covid_Research/Sars_Cov_Spike_RNA_Human/ncbi_dataset/data/cds.fna")
##data updated as of Feb 10, 2020

#meta=fread("D:/Covid_Research/meta2.csv") #not needed

summary(width(genomes))

widest=which(max(width(genomes))==width(genomes))
PA=pairwiseAlignment(pattern=genomes,subject=genomes[widest[1]],
                     type="global")
aligned_genomes=alignedPattern(PA)
aligned_genomes=stackStrings(aligned_genomes,Lpadding.letter="-",
                             Rpadding.letter="-",from=1, to=max(width(aligned_genomes)))
summary(width(aligned_genomes))

library(Matrix)
#convert to binary
bin_mat=NULL
for(j in (1:length(aligned_genomes))){
  g1=as.character(aligned_genomes[j])
  tot_bin=NULL
  for(i in (1:width(g1))){
    sub_it=substr(g1,i,i)
    
    if(sub_it=="A"){bin_it=c(1,0,0,0)}
    if(sub_it=="T"){bin_it=c(0,1,0,0)}
    if(sub_it=="C"){bin_it=c(0,0,1,0)}
    if(sub_it=="G"){bin_it=c(0,0,0,1)}
    if(sub_it=="-"){bin_it=c(0,0,0,0)}
    
    tot_bin=c(tot_bin,bin_it)
  }
  tot_bin=Matrix(tot_bin,nrow=1)
  bin_mat=rbind(bin_mat,tot_bin)
  
}


library(irlba)
col_means=apply(bin_mat,2,mean)
#col_sds=apply(bin_mat,2,sd) #many columns are constant

svd_out=irlba(bin_mat,nv=50,center=col_means)

svd_out$d


library(ggplot2)

scores=svd_out$u[,c(1:3)]

min1=which(min(scores[,1])==scores[,1])

names(aligned_genomes[min1])

max1=which(max(scores[,1])==scores[,1])

names(aligned_genomes[max1])

min2=which(min(scores[,2])==scores[,2])

names(aligned_genomes[min2])

max2=which(max(scores[,2])==scores[,2])

names(aligned_genomes[max2])

plot(scores)


clus_out=kmeans(scores,4,iter.max = 100,nstart=50)
clus_out$size
clus_out$tot.withinss/clus_out$totss

d=as.data.frame(names(aligned_genomes))
d$clus=clus_out$cluster
fwrite(d,"D:/Covid_Research/sars_cov_cluster_out.csv")

