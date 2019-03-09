### ---------------
###
### Create: Jianming Zeng
### Date: 2019-03-09 19:20:31
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-09  First version
###
### ---------------

options(stringsAsFactors = F)
GTEx=read.table('~/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz'
             ,header = T,sep = '\t',skip = 2)
GTEx[1:4,1:4]
h=head(GTEx)
save(h,file = 'GTEx_head.Rdata')
save(GTEx,file = 'GTEx_all.Rdata')

