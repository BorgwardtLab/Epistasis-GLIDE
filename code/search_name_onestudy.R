
search.name<-function(name_vector_set1,
		      name_vector_set2,
                      PartId,
                      vector_index1,
		      vector_index2,
                      nb_chunk_set1=num_chunk1,
		      nb_chunk_set2=num_chunk2,
                      snp_chunksize_set1=chunksize1,
		      snp_chunksize_set2=chunksize2){
  snp_len_set1 <- dim(name_vector_set1)[1]
  snp_len_set2 <- dim(name_vector_set2)[1]
  
  snp1_indexstart <- (PartId[,1])*snp_chunksize_set1 
  snp2_indexstart <- (PartId[,2])*snp_chunksize_set2 

#  if (PartId[,1]==nb_chunk_set1 && PartId[,2]==nb_chunk_set2)
#    {snp_chunksize_set1<-(snp_len_set1-snp1_indexstart)+1;
#     snp_chunksize_set2<-(snp_len_set2-snp2_indexstart)+1;}
  snp1_relindex <- (vector_index1+1);     
  snp2_relindex <- (vector_index2+1);  	
  return(data.frame(SNP1 = name_vector_set1[snp1_indexstart + snp1_relindex,1],
		    SNP2 = name_vector_set2[snp2_indexstart + snp2_relindex,1]))

}

args = commandArgs(TRUE)
filename=args[1]
Set1=args[2]
Set2=args[3]
chunksize1=as.numeric(args[4])
chunksize2=as.numeric(args[5])
nsubjects=as.numeric(args[6])
Resultsnames=as.character(args[7])
BlockSize=as.numeric(args[8])
as.data.frame(read.table(file=filename),stringsAsFactors=FALSE)->resmatrix


name.vector_set1<-read.table(Set1,stringsAsFactors=FALSE,header=FALSE)[1]
name.vector_set2<-read.table(Set2,stringsAsFactors=FALSE,header=FALSE)[1]
PartID<-resmatrix[,c(1,2)]
vector_index1<-resmatrix[,3]*BlockSize+resmatrix[,5]
vector_index2<-resmatrix[,4]*BlockSize+resmatrix[,6]
num_chunk1<-ceiling(dim(name.vector_set1)[1]/chunksize1)
num_chunk2<-ceiling(dim(name.vector_set2)[1]/chunksize2)


search.name(name.vector_set1,name.vector_set2,PartID,vector_index1,vector_index2)->resmatrix[,11:12]
2*pt(-abs(resmatrix[,7]),df=nsubjects-4)->resmatrix[,13]
2*pt(-abs(resmatrix[,8]),df=nsubjects-4)->resmatrix[,14]
2*pt(-abs(resmatrix[,9]),df=nsubjects-4)->resmatrix[,15]
2*pt(-abs(resmatrix[,10]),df=nsubjects-4)->resmatrix[,16]

#resmatrix<-resmatrix[-(which(as.character(resmatrix[,9])==as.character(resmatrix[,10]))),]
colnames(resmatrix)<-c('P1','P2','bidx','bidy','tidx','tidy','Tint','TSnp1','TSnp2','TSnp1n2','Snp1','Snp2','Pint','PSnp1','PSnp2','PSnp1n2')
write.table(resmatrix,file=Resultsnames,quote=FALSE,row.names=FALSE,col.names=T)