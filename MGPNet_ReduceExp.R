# Ertugrul Dalgic, PhD, 2023-24
# Reduce gene exp files acc. to gene id with the max value
# Reduce tumor exp file based on normal exp file

# download 
#download.file(url='https://www.proteinatlas.org/download/rna_tissue_hpa.tsv.zip',
#              destfile='/truba/home/edalgic/MAPP/rna_tissue_hpa.tsv.zip')
#download.file(url='https://www.proteinatlas.org/download/rna_cancer_sample.tsv.zip',
#              destfile='/truba/home/edalgic/MAPP/rna_cancer_sample.tsv.zip')
#unzip('rna_tissue_hpa.tsv.zip')
#unzip('rna_cancer_sample.tsv.zip')

# use parallel library
library(doParallel)
ncores = detectCores()
clcores = makeCluster(ncores-1)
registerDoParallel(clcores)

###############################
# input gene regulatory network
###############################
library(dorothea)
library(decoupleR)
ps1 <- decoupleR::get_collectri(split_complexes = TRUE)
givenintf = 'EdwinWang_HumanSignalingNet_v6.csv'
ops = read.csv(givenintf,header = FALSE,colClasses = c('character'))
# retain only regulatory ints (Pos or Neg)
ps2 = ops[ops$V5=='Pos'|ops$V5=='Neg',]
# get all the network genes of both datasets
px1 = c(ps1$source,ps1$target)
px2 = c(ps2$V2,ps2$V4)
px = union(px1,px2)


#########################
# Normal RNA expression
#########################
# input hpa rna expression data and compare cors to network values
# reduce rna data only to the proteins in the network
hpaexpdtf = 'rna_tissue_hpa.tsv'
hredtx = read.delim(hpaexpdtf)

# Reduce by network genes
hredtx = hredtx[hredtx$Gene.name %in% px,]

# select by the max TPM
# for repeated values for the same gene in same tissue 
# there could be more ensembl id genes for a gene name; ignore
hregs = unique(hredtx$Gene.name)
hredt <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(hgi$Tissue)){
        hgit = hgi[hgi$Tissue==ti,]
        if (nrow(hgit)==1){
            hx = rbind(hx,hgit)}
        if (nrow(hgit)>1){
            # select max TPM row
            # there might be more than one max just select one of them
            hx = rbind(hx,hgit[which(hgit$TPM==max(hgit$TPM))[1],])}}
    cbind(hx)}
hredtrna = as.data.frame(hredt)
names(hredtrna)=names(hredtx)
write.table(hredtrna,file="Reduced_Normal_RNA_exp_selbyMAX.txt",sep='\t',
            row.names=FALSE,quote=FALSE)


# calculate mean and stdev of multiple tpm values per gene in same tissue 
hregs = unique(hredtx$Gene.name)
hredt <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(as.character(hgi$Tissue))){
        hgit = hgi[hgi$Tissue==ti,]
        # replace with mean value also calculate std dec
        hx = rbind(hx,c(hgit[1,c(1,2,3)],mean(hgit$TPM),sd(hgit$TPM)))}
    cbind(hx)}

hredtrna = as.data.frame(hredt)
names(hredtrna)= c('Gene','Gene.name','Tissue','Mean.exp','Std.exp')
# replace NA with 0 (in std column)
hredtrna[is.na(hredtrna)]=0
# convert all values to character to write next
hredtrna = apply(hredtrna,2,as.character)
write.table(hredtrna,file="Reduced_Normal_RNA_exp.txt",sep='\t',
            row.names=FALSE,quote=FALSE)


##################################
# Tumor RNA expression MEDIAN data
##################################
# input TCGA median data (precalculated in a different script) which has gene names
# calculate mean exp value per gene (for more than one Ensembl ID per gene symbol)
hpaexpdtf = 'TCGA_Cancer_RNA_Median_exp.txt'
hredtxtu = read.delim(hpaexpdtf)
hregs = unique(as.character(hredtxtu$Gene.name))
hredtrnatu <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtxtu[hredtxtu$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(as.character(hgi$Cancer))){
        hgit = hgi[hgi$Cancer==ti,]
        # replace with mean value also calculate std dec
        hx = rbind(hx,c(hgit[1,c(1,2,3)],mean(hgit$Median.exp),sd(hgit$Median.exp)))}
    cbind(hx)}
names(hredtrnatu)=c('Gene','Gene.name','Tissue','Mean.exp','Std.exp')
# replace NA with 0 (in std column)
hredtrnatu[is.na(hredtrnatu)]=0
# convert all values to character to write next
hredtrnatu = apply(hredtrnatu,2,as.character)
write.table(hredtrnatu,file="Reduced_Cancer_RNA_exp.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

#######################################
# Tumor RNA expression original dataset
#######################################
# NOT USED
# input hpa tcga cancer rna expression data and compare cors to network values
# reduce rna data only to the genes in the normal dataset and in the network
# and add Gene names
hpaexpdtf = 'rna_cancer_sample.tsv'
ohredtx = read.delim(hpaexpdtf)
ohredtx = ohredtx[ohredtx$Gene %in% hredtx$Gene,]

oh <- foreach(i=1:nrow(ohredtx),.combine=rbind) %dopar% {
    igname = unique(hredtx[hredtx$Gene==ohredtx[i,1],]$Gene.name)
    cbind(ohredtx[i,],igname)}
oh = as.data.frame(oh)
colnames(oh) = c(colnames(ohredtx),'Gene.name')


# calculate mean exp value per gene
# for repeated values for the same gene in same cancer type 
# there could be more ensembl id genes and multiple samples for a gene name; ignore
hregs = unique(as.character(oh$Gene.name))
ohmean <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = oh[oh$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(as.character(hgi$Cancer))){
        hgit = hgi[hgi$Cancer==ti,]
        # replace with mean value also calculate std dec
        hx = rbind(hx,c(hgit[1,c(1,5,3)],mean(hgit$FPKM),sd(hgit$FPKM)))}
    cbind(hx)}
names(ohmean)=c('Gene','Gene.name','Tissue','Mean.exp','Std.exp')
write.table(ohmean,file="Reduced_Cancer_RNA_exp.txt",sep='\t',
            row.names=FALSE,quote=FALSE)


###########################
# Normal Protein Expression
###########################

ohredtx = read.delim('normal_tissue.tsv')
# clear NA rows
hredtx = ohredtx[ohredtx$Tissue!='N/A',]
hredtx = hredtx[hredtx$Level!='N/A',]
# reduce normal tissue protein data only to the proteins in the network
hredtx = hredtx[hredtx$Gene.name %in% px,]
# convert Level categories to numbers: Not det. 1, Low 2, Medium 3, High 4
hredtx[hredtx$Level=='Not detected',]$Level=1
hredtx[hredtx$Level=='Low',]$Level=2
hredtx[hredtx$Level=='Medium',]$Level=3
hredtx[hredtx$Level=='High',]$Level=4
# ignore other Levels such as ascending ... etc
hredtx = hredtx[hredtx$Level %in% c(1,2,3,4),]
hredtx$Level=as.numeric(hredtx$Level)

# select by the max Level
# for repeated values for the same gene in same tissue and same cell type
# there could be more ensembl id genes for a gene name; ignore
hregs = unique(hredtx$Gene.name)
hredt <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(hgi$Tissue)){
        hgit = hgi[hgi$Tissue==ti,]
        if (nrow(hgit)==1){
            hx = rbind(hx,hgit)}
        if (nrow(hgit)>1){
            for (ci in unique(hgit$Cell.type)){
                hgitc = hgit[hgit$Cell.type==ci,]
                if (nrow(hgitc)==1){
                    hx = rbind(hx,hgitc)}
                if (nrow(hgitc)>1){
                    # select max TPM row
                    # there might be more than one max just select one of them
                    hx = rbind(hx,hgitc[which(hgitc$Level==max(hgitc$Level))[1],])}}}}
    cbind(hx)}
hredtprot = as.data.frame(hredt)
names(hredtprot)=names(hredtx)
write.table(hredtprot,file="Reduced_Normal_Protein_exp_selbyMAX.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

# calculate mean and std for repeated values per gene for each tissue
hregs = unique(hredtx$Gene.name)
hredt <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(hgi$Tissue)){
        hgit = hgi[hgi$Tissue==ti,]
        # replace with mean value also calculate std dec
        hx = rbind(hx,c(hgit[1,c(1,2,3)],mean(hgit$Level),sd(hgit$Level)))}
    cbind(hx)}

hredtprot = as.data.frame(hredt)
names(hredtprot)= c('Gene','Gene.name','Tissue','Mean.exp','Std.exp')
# replace NA with 0 (in std column)
hredtprot[is.na(hredtprot)]=0
# convert all values to character to write next
hredtprot = apply(hredtprot,2,as.character)
write.table(hredtprot,file="Reduced_Normal_Protein_exp.txt",sep='\t',
            row.names=FALSE,quote=FALSE)


############################
# Cancer Protein Expression
############################

###########################
# Cancer protein expression
# input hpa tumor protein expression data
# reduce tumor tissue protein data only to the proteins in the network
# clean dataset
tu = read.delim('pathology.tsv')
tu = na.omit(tu[,1:7])
# remove rows with less than 6 values
sus = apply(tu[,4:7],1,'sum')
tu = tu[which(sus>=6),]
# convert to a single value
# high = 4, medium = 3, low = 2, not.detected = 1
tu = cbind(tu,'Overall'=(4*tu$High+3*tu$Medium+2*tu$Low+1*tu$Not.detected)/
               (tu$High+tu$Medium+tu$Low+tu$Not.detected))
hredtx = tu[tu$Gene.name %in% px,]
hredtx$Overall=as.numeric(hredtx$Overall)

# select the max value 
# for repeated values for the same gene in same tissue 
# there could be more ensembl id genes for a gene name; ignore
hregs = unique(hredtx$Gene.name)
hredtprottu <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(hgi$Cancer)){
        hgit = hgi[hgi$Cancer==ti,]
        if (nrow(hgit)==1){
            hx = rbind(hx,hgit)}
        if (nrow(hgit)>1){
            # select max Overall row
            # there might be more than one max just select one of them
            hx = rbind(hx,hgit[which(hgit$Overall==max(hgit$Overall))[1],])}}
    cbind(hx)}
write.table(hredtprottu,file="Reduced_Tumor_Protein_exp_selbyMAX.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

# calculate mean and std for repeated values per gene for each tissue
hregs = unique(hredtx$Gene.name)
hredt <- foreach(hi=1:length(hregs),.combine=rbind) %dopar% {
    hgi = hredtx[hredtx$Gene.name==hregs[hi],]
    hx=c()
    for (ti in unique(hgi$Cancer)){
        hgit = hgi[hgi$Cancer==ti,]
        # replace with mean value also calculate std dec
        hx = rbind(hx,c(hgit[1,c(1,2,3)],mean(hgit$Overall),sd(hgit$Overall)))}
    cbind(hx)}

hredtprottu = as.data.frame(hredt)
names(hredtprottu)= c('Gene','Gene.name','Tissue','Mean.exp','Std.exp')
# replace NA with 0 (in std column)
hredtprottu[is.na(hredtprottu)]=0
# convert all values to character to write next
hredtprottu = apply(hredtprottu,2,as.character)
write.table(hredtprottu,file="Reduced_Tumor_Protein_exp.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

stopImplicitCluster()


