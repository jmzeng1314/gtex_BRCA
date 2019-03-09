rm(list = ls())  
options(stringsAsFactors = F)
# The current release is V7 including 11,688 samples, 53 tissues and 714 donors.

if(F){
  options(stringsAsFactors = F)
  GTEx=read.table('~/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz'
                  ,header = T,sep = '\t',skip = 2)
  GTEx[1:4,1:4]
  h=head(GTEx)
  save(h,file = 'GTEx_head.Rdata')
}
if(F){
  load('~/Desktop/GTEx_all.Rdata')
  a[1:4,1:4]
  colnames(a)
  # SMTS	Tissue Type, area from which the tissue sample was taken.   
  # SMTSD	Tissue Type, more specific detail of tissue type
  b=read.table('GTEx_v7_Annotations_SampleAttributesDS.txt',
               header = T,sep = '\t',quote = '')
  table(b$SMTS) 
  breat_gtex=a[,gsub('[.]','-',colnames(a)) %in% b[b$SMTS=='Breast',1]]
  rownames(breat_gtex)=a[,1]
  dat=breat_gtex
  ids=a[,1:2]
  head(ids)
  colnames(ids)=c('probe_id','symbol')
  dat=dat[ids$probe_id,]
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  dat[1:4,1:4]  
  breat_gtex=dat
  save(breat_gtex,file = 'breat_gtex_counts.Rdata')
}
load(file = 'breat_gtex_counts.Rdata')
dat=log2(edgeR::cpm(breat_gtex)+1)
dat[1:4,1:4]
if(T){
  ddata=t(dat)
  ddata[1:4,1:4]
  s=colnames(ddata);head(s)
  library(org.Hs.eg.db)
  s2g=toTable(org.Hs.egSYMBOL)
  g=s2g[match(s,s2g$symbol),1];head(g)
  #  probe Gene.symbol Gene.ID
  dannot=data.frame(probe=s,
                    "Gene.Symbol" =s, 
                    "EntrezGene.ID"=g)
  ddata=ddata[,!is.na(dannot$EntrezGene.ID)]
  dannot=dannot[!is.na(dannot$EntrezGene.ID),] 
  head(dannot)
  library(genefu)
  # c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow")
  
  s<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                         annot=dannot,do.mapping=TRUE)
  table(s$subtype)
  tmp=as.data.frame(s$subtype)
  subtypes=as.character(s$subtype)
}

library(genefu)
pam50genes=pam50$centroids.map[c(1,3)]
pam50genes[pam50genes$probe=='CDCA1',1]='NUF2'
pam50genes[pam50genes$probe=='KNTC2',1]='NDC80'
pam50genes[pam50genes$probe=='ORC6L',1]='ORC6'
x=dat
x=x[pam50genes$probe[pam50genes$probe %in% rownames(x)] ,]
tmp=data.frame(subtypes=subtypes)
rownames(tmp)=colnames(x)
library(pheatmap)


pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_raw.png') 


x=t(scale(t(x)))
x[x>1.6]=1.6
x[x< -1.6]= -1.6

pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_scale.png') 






