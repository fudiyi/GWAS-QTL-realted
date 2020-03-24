# QTL分析-Multiple QTL Mapping

#### 相关资料：https://rdrr.io/cran/qtl/man/MQM.html

#### 相关文件：A guide to QTL mappingwith R.pdf; MQM-tour 附件已上传Gihub [GWAS-QTL-realted](https://github.com/fudiyi/GWAS-QTL-realted)

#### 原始数据 qtl_data.csv：可通过“python-处理hmp文件用于QTL分析”得到，保持格式与其一致即可

```R
# R 
# FDY-2020-3-24

> setwd("E:/QTL/fdy")
> library(qtl)
> geno = read.cross("csv",  file="qtl_data.csv",na.strings = "-",genotypes = c("AA", "AB","BB"), alleles = c("A", "B"), estimate.map = FALSE, BC.gen=2, F.gen=3) #读取基因型文件，群体为BC2F3，若子代基因型与亲本一致则修改基因型为“AA/BB”（male:AA;female:BB），若为杂合则修改为“AB”，若为测到则修改为“-”
#result
--Read the following data:
	 131  individuals
	 1011  markers
	 1  phenotypes
 --Cross type: bcsft 
Warning message:
In summary.cross(cross) :
  Some markers at the same position on chr 1,2,3,4,5,6,7,8,9,10; use jittermap()
> geno <- jittermap(geno) #文件中会有遗传距离一样的marker，使用jittermap()
> summary(geno)
#result
	BC(2)F(3) cross

    No. individuals:    131 

    No. phenotypes:     1 
    Percent phenotyped: 100 

    No. chromosomes:    10 
        Autosomes:      1 2 3 4 5 6 7 8 9 10 

    Total markers:      1011 
    No. markers:        129 96 81 89 122 106 87 187 39 75 
    Percent genotyped:  97.8 
    Genotypes (%):      AA:42.5  AB:12.3  BB:45.3  not BB:0.0  not AA:0.0 
> geno <- calc.genoprob(geno, step=0.5, off.end=0, error.prob=0.0001, map.function="kosambi", stepwidth="fixed") #calculate QTL conditional probability
> geno <- convert2riself(geno) #convert data to RIL format (not documented, but see rqtl.org/faq)
> operm <- scanone(geno, n.perm=1000, pheno.col="shoushang",method="em")
> operm[1:5]
> plot(operm)
> summary(operm, alpha=c(0.20, 0.05))
#result
LOD thresholds (1000 permutations)
     lod
20% 3.37
5%  4.16

################## QTL mapper, including scanone and mqm ##################

#进行mqm分析，自定义一个函数来进行相关判断
mapper= function(geno, pheno.col, threshold ) {

##################### primary single-QTL scanning ############################ 
out<- scanone (geno,pheno.col=pheno.col,method="hk")
#out<- scanone (geno,pheno.col=pheno.col,method="hk",addcovar=covar)
ylimit= max(as.data.frame(out[,-(1:2)]))
ylimit

######### plot scanone results #######
pdf(paste(pheno.col,"_scanone.pdf",sep=""), width = 6, height=4)
plot(out, lodcolumn=1, lty=1, col=2,ylim=c(0,ceiling(ylimit/10)*10),ylab="LOD score" )
#legend("topright", pheno.col, lty=1, col=2, cex=0.6)
title(main = pheno.col)
dev.off()

##### scanone results summary ####
single.summary=summary(out,threshold=threshold,lodcolumn=1) 
single.summary

######### output scanone results ########
write.csv(out, file=paste(pheno.col,"_scanone.csv",sep=""), quote=FALSE,row.names=T)
    
######### output scanone summarized results #########
write.csv(single.summary, file=paste(pheno.col,"_scanone_summary.csv",sep=""), quote=FALSE,row.names=T)
    
##################### end of single-QTL mapping ############################



######################### beginning of mqm mapping #####################
##### test if QTL mapped by scanone #####
    
if(dim(single.summary)[1] !=0 ){
looper=1
##### make initial QTL object for fitting multiple QTL model #####
chr= single.summary$chr
pos= single.summary$pos
qtl= makeqtl(geno, chr=chr, pos=pos, what="prob")
qtl

while(looper ==1){

##### construct formula of multiple QTL model #####
term= rownames(summary(qtl))
formula1= as.formula(paste("y ~ ", paste(term, collapse= "+")))
formula1

##### refine QTL model #####
qtl= refineqtl(geno, pheno.col=pheno.col, qtl=qtl, formula= formula1, method= "hk")
qtl

##### scan whole genome to check if there are additional significant QTL #####
out.aq <- addqtl(geno, pheno.col=pheno.col, qtl=qtl, formula=formula1, method="hk")
single.summary= summary(out.aq,threshold=threshold,lodcolumn=1 )
single.summary

if(dim(single.summary)[1] !=0 ){

#### if there does more significant QTLs are detected, add into previous QTL object ####
rqtl.added= addtoqtl(geno, qtl, single.summary$chr, single.summary$pos)
qtl= rqtl.added

looper=1

}else{
looper=0
}

}

##### plot final multiple QTL mapping results #####
pdf(paste(pheno.col,"_mqm.pdf",sep=""), width = 6, height=4)
plotLodProfile(qtl,main=pheno.col,ylab="LOD",cex=0.4,lwd=1)
dev.off()

##### output multiple QTL fitting results #####
out.mqm= fitqtl(geno, pheno.col=pheno.col, qtl=qtl, formula= formula1, method= "hk", dropone=T, get.ests=T)
summary(out.mqm)

write.table(summary(out.mqm)[[1]], file= paste(pheno.col,"_mqm_full.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)
write.table(summary(out.mqm)[[2]], file= paste(pheno.col,"_mqm_dropone.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)
write.table(summary(out.mqm)[[3]], file= paste(pheno.col,"_mqm_effect.txt",sep=""), quote=FALSE, sep='\t', row.names=TRUE, col.names=T)

######### Here is QTL scanning results with other QTL fixed as background #########
names(attributes(qtl))                                                                       
lodprof <- attr(qtl, "lodprofile")  # results saved in the component lodprofile                                                                                    
lodprof[[1]] # first QTL scanning results                                                                                                           
lodint(qtl, qtl.index=1)  ##### support interval for the first QTL                                                                                                                  ######## output scanning results for each QTL, saved in lodprof.txt #######
out.lodprof<-capture.output(lodprof)
cat(out.lodprof,file=paste(pheno.col,"_mqm_lodprof.txt",sep=""),sep="\n",append=TRUE)

######## output QTL support interval for each QTL, saved in qtlint.txt #########
out.int=NULL #The purpose here is to empty the possible previously existing file
cat(out.int,file=paste(pheno.col,"_mqm_qtlint.txt",sep=""),sep="\n")
for (i in 1:length(lodprof)){                                                                                                       
qtlint= lodint(qtl, qtl.index=i)
out.int<-capture.output(qtlint)
cat(out.int,file=paste(pheno.col,"_mqm_qtlint.txt",sep=""),sep="\n",append=TRUE)
}

}
}

######################### end of qtl mapper #####################

################################ analyze your traits #############################

##################### for singele trait, input your trait name ###################
> mapper(geno=geno, pheno.col="phenotype", threshold=2)

########### for multiple traits, use a loop #############
# input your traits list
pheno.list=names(geno$pheno)
# loop through your traits
for(i in 1:length(pheno.list)){
mapper(geno=geno, pheno.col=pheno.list[i], threshold=2)
}
################################ end of mapping ###############################
```

