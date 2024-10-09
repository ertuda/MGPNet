# Ertugrul Dalgic, PhD, 2023-24
# Analyze pairs of proteins with mutual relationship
# Edwin Wang protein signaling network
# function to get pairs
# PP: pospos: pair with mutual positive direct effect
# NN: negneg: pair with mutual negative direct effect
# PN: posneg: pair with mutual effect but one is pos and the other is neg

###########
# FUNCTIONS
###########

# use parallel library
library(doParallel)
ncores = detectCores()
clcores = makeCluster(ncores-1)
registerDoParallel(clcores)

# function to collect mutual pairs
allpas <- function(givps){
    # mutual pairs PP, NN, PN
    # there is no order for PP and NN
    ppt = c()
    nnt = c()
    # for PN order: g1 g2 pos, g2 g1 neg
    pnt = c()
    # pos or neg pairs with no counter interaction
    posnot = c()
    negnot = c()
    # self pairs with pos or neg effect
    sfpost = c()
    sfnegt = c()
    # check all pairs
    for (i in 1:(nrow(givps)-1)){
        g1 = givps[i,2]
        g2 = givps[i,4]
        sign12 = givps[i,5]
        # compare with the rest (givepss) to get mutual relationships
        givpss = givps[(i+1):nrow(givps),]
        if (g1 != g2){
            if (sign12=='Pos'){
                if(nrow(givpss[givpss$V2==g2&givpss$V4==g1,])==0){
                    posnot = rbind(posnot,c(g1,g2))}
                else if(nrow(givpss[givpss$V2==g2&givpss$V4==g1&givpss$V5=='Pos',])>0){
                    ppt = rbind(ppt,c(g1,g2))}
                else if(nrow(givpss[givpss$V2==g2&givpss$V4==g1&givpss$V5=='Neg',])>0){
                    pnt = rbind(pnt,c(g1,g2))}}
            if (sign12=='Neg'){
                if(nrow(givpss[givpss$V2==g2&givpss$V4==g1,])==0){
                    negnot = rbind(negnot,c(g1,g2))}
                else if(nrow(givpss[givpss$V2==g2&givpss$V4==g1&givpss$V5=='Pos',])>0){
                    pnt = rbind(pnt,c(g2,g1))}
                else if(nrow(givpss[givpss$V2==g2&givpss$V4==g1&givpss$V5=='Neg',])>0){
                    nnt = rbind(nnt,c(g1,g2))}}}
        if (g1 == g2){
            if (sign12=='Pos'){
                sfpost = rbind(sfpost,c(g1,g2))}
            if (sign12=='Neg'){
                sfnegt = rbind(sfnegt,c(g1,g2))}}}
    
    # last pair might not be added if not mutual so check it
    lg1 = givps[nrow(givps),2]
    lg2 = givps[nrow(givps),4]
    lsign12 = givps[nrow(givps),5]
    # compare with the rest
    givpss = givps[1:(nrow(givps)-1),]
    if (lg1 != lg2){
        # avoid PN,PP,NN
        if(nrow(givpss[givpss$V2==lg2&givpss$V4==lg1,])==0){ 
            if (lsign12=='Pos'){ 
                posnot = rbind(posnot,c(lg1,lg2))}
            if (lsign12=='Neg'){
                negnot = rbind(negnot,c(lg1,lg2))}}}
    if (lg1 == lg2){
        if (lsign12=='Pos'){
            sfpost = rbind(sfpost,c(g1,g2))}
        if (lsign12=='Neg'){
            sfnegt = rbind(sfnegt,c(g1,g2))}}
    
    retlistt = list(ppt,nnt,pnt,posnot,negnot,sfpost,sfnegt)
    names(retlistt) = c('PP','NN','PN','PO','NE','SP','SN')
    return(retlistt)}

# function to get correlation values of given list of pairs in the given exp datasets
corpas <- function(tpa,tpai1,tpai2,txp,txpcomx){
    # return a table with pairs with their exp cor vals
    # get the given pair with indices for tpa; i1-2: genes/prots
    # get the given exp data 
    # txpcomx is threshold for exp vals
    txps = unique(txp$Gene.name)
    cpt <- foreach(j=1:nrow(tpa),.combine=rbind) %dopar% {
        #for (j in 1:nrow(tpa)){
        p1 = tpa[j,tpai1]
        p2 = tpa[j,tpai2]
        cptx = c(p1,p2,NA,NA)
        if (p1 %in% txps & p2 %in% txps){
            txp1 = txp[txp$Gene.name == p1,]
            txp2 = txp[txp$Gene.name == p2,]
            comtis = intersect(txp1$Tissue,txp2$Tissue)
            # compare two lists if there is enough vals at least txpcomx threshold
            if (length(comtis)>=txpcomx){
                exv1 = as.numeric(txp1[which(txp1$Tissue %in% comtis),]$Mean.exp)
                exv2 = as.numeric(txp2[which(txp2$Tissue %in% comtis),]$Mean.exp)
                cox12 = cor.test(exv1,exv2,method='spearman')
                cptx = c(p1,p2,cox12$estimate,cox12$p.value)}}}
    # remove NA correlation rows
    cpt = cpt[!is.na(cpt[,3]),]
    # add FDR adjusted p values
    cpt = cbind(cpt,p.adjust(as.numeric(cpt[,4]),method = 'fdr'))
    return(cpt)}

# function to calculate degree values of directed network
degpas <- function(degps,degset,pi1,pi2,pi3){
    # return a table of genes with degree vals
    degtb = c()
    # avoid self interactions
    degps2 = c()
    for (i in 1:nrow(degps)){
        if (degps[i,pi1] != degps[i,pi2]){
            degps2 = rbind(degps2,degps[i,])}}
    degps2 = as.data.frame(degps2)
    # all degree combinations: directed/undirected, pos/neg
    # get the given pair with indices; pi1-2: genes/prots,pi3:pos/neg
    deall = unique(c(degps[,pi1], degps[,pi2]))
    # all network or connections of a subset with all network are taken
    if (length(degset) == 1){
        if (degset == 'all'){
            degset = deall}}
    for (pii in degset){
        pi = as.character(pii) 
        netin = degps2[degps2[,pi2]==pi,]
        netinpos = netin[netin[,pi3]=='Pos',]
        netinneg = netin[netin[,pi3]=='Neg',]
        netout = degps2[degps2[,pi1]==pi,]
        netoutpos = netout[netout[,pi3]=='Pos',]
        netoutneg = netout[netout[,pi3]=='Neg',]
        # output numbers as a row:
        # total,in,out,pos,neg,inpos,inneg,outpos,outneg
        degtb = rbind(degtb,c(pi,(nrow(netin)+nrow(netout)),
                nrow(netin),nrow(netout),
                (nrow(netinpos)+nrow(netoutpos)),
                (nrow(netinneg)+nrow(netoutneg)),
           nrow(netinpos),nrow(netinneg),nrow(netoutpos),nrow(netoutneg)))}
    degtb = as.data.frame(degtb)
    colnames(degtb) = c('gene','total_deg','in_deg','out_deg','pos_deg',
                        'neg_deg','inpos_deg','inneg_deg','outpos_deg',
                        'outneg_deg')
    return(degtb)}

# function to calculate interactions between mutual pairs i.e PP-PN ints, etc.
intpas <- function(gpas){
    # use allpas out to calculate ints among and btw PP,PN,NN
    # table of common genes with all pair combinations obtained my multiplication
    # i.e, if gene a belongs to 2 PP and 3 PN then there are 3*2 PP-PN ints via gene a
  
    # vectorized PP, PN, NN member genes
    gpp = c(gpas[['PP']][,1],gpas[['PP']][,2])
    gpn = c(gpas[['PN']][,1],gpas[['PN']][,2])
    gnn = c(gpas[['NN']][,1],gpas[['NN']][,2])
  
    # PP-PN 
    cmpppn = intersect(gpp,gpn)
    tpppn = table(gpp[gpp %in% cmpppn])
    tpnpp = table(gpn[gpn %in% cmpppn])
    # table of genes with number of ints
    intpppn = c()
    for (j in cmpppn){
        intpppn = rbind(intpppn,c(j,'PP-PN',tpppn[[j]]*tpnpp[[j]]))}
    
    # PP-NN 
    cmppnn = intersect(gpp,gnn)
    tppnn = table(gpp[gpp %in% cmppnn])
    tnnpp = table(gnn[gnn %in% cmppnn])
    # table of genes with number of ints
    intppnn = c()
    for (j in cmppnn){
      intppnn = rbind(intppnn,c(j,'PP-NN',tppnn[[j]]*tnnpp[[j]]))}
    
    # PN-NN 
    cmpnnn = intersect(gpn,gnn)
    tpnnn = table(gpn[gpn %in% cmpnnn])
    tnnpn = table(gnn[gnn %in% cmpnnn])
    # table of genes with number of ints
    intpnnn = c()
    for (j in cmpnnn){
      intpnnn = rbind(intpnnn,c(j,'PN-NN',tpnnn[[j]]*tnnpn[[j]]))}

    # self ints PP-PP ints , etc.  
    # PP-PP
    tpppp = table(gpp)
    intpppp = c()
    for (j in unique(gpp)){
        intpppp = rbind(intpppp,c(j,'PP-PP',tpppp[[j]]*tpppp[[j]]))}
    
    # PP-PP
    tpnpn = table(gpn)
    intpnpn = c()
    for (j in unique(gpn)){
        intpnpn = rbind(intpnpn,c(j,'PN-PN',tpnpn[[j]]*tpnpn[[j]]))}
    
    # PP-PP
    tnnnn = table(gnn)
    intnnnn = c()
    for (j in unique(gnn)){
        intnnnn = rbind(intnnnn,c(j,'NN-NN',tnnnn[[j]]*tnnnn[[j]]))}
    
    intmps = as.data.frame(rbind(intpppn,intppnn,intpnnn,
                                 intpppp,intpnpn,intnnnn))
    colnames(intmps) = c('Gene','Int_type','No')
    return(intmps)}
    
#############################
# SIGNED INTERACTION NETWORK
#############################

# input protein signaling network
givenintf = 'EdwinWang_HumanSignalingNet_v6.csv'
ops = read.csv(givenintf,header = FALSE,colClasses = c('character'))

# retain only regulatory ints (Pos or Neg)
ps = ops[ops$V5=='Pos'|ops$V5=='Neg',]
callpas = allpas(ps)

#################
# RNA EXPRESSION
#################

# input hpa normal rna expression data and hpa tcga cancer rna expression data
# use only reduced rna dataset where there is a single value per gene and tissue
hredtrna = read.delim('Reduced_Normal_RNA_exp.txt')
hredtrnatu = read.delim("Reduced_Cancer_RNA_exp.txt")
# remove highly varied genes based on cov (coeff of var. = std/mean*100)
rnavarthr = 1500000000000000000000000000000000000
hredtrna$Cov = hredtrna$Std.exp/hredtrna$Mean.exp*100
hredtrnatu$Cov = hredtrnatu$Std.exp/hredtrnatu$Mean.exp*100
#hredtrna = hredtrna[hredtrna$Cov<rnavarthr,]
#ignore cov for 0 mean genes
hredtrnatu0mean = hredtrnatu[hredtrnatu$Mean.exp==0,]
hredtrnatu = hredtrnatu[hredtrnatu$Cov<rnavarthr,]
hredtrnatu = rbind(hredtrnatu[,1:5],hredtrnatu0mean[,1:5])
#remove NAs
hredtrna = na.omit(hredtrna)
hredtrnatu = na.omit(hredtrnatu)
# use the same set of genes for normal and Cancer
commongs = intersect(unique(hredtrna$Gene.name),unique(hredtrnatu$Gene.name))
hredtrna = hredtrna[hredtrna$Gene.name %in% commongs,]
hredtrnatu = hredtrnatu[hredtrnatu$Gene.name %in% commongs,]

# For correlation analysis require at least a certain number of values per gene
minnovs = 5

# NORMAL CORRELATIONS
# exclude self pairs
# for pairs excluding PP,PN,NN
corponerna = corpas(rbind(callpas[['PO']],callpas[['NE']]),1,2,hredtrna,minnovs)
# for negatives excluding NN,PN
cornerna = corpas(callpas[['NE']],1,2,hredtrna,minnovs)
# for PP
corpprna = corpas(callpas[['PP']],1,2,hredtrna,minnovs)
# for positives excluding PP
corpopnrna = corpas(rbind(callpas[['PO']],callpas[['PN']]),1,2,hredtrna,minnovs)
# for NN
cornnrna = corpas(callpas[['NN']],1,2,hredtrna,minnovs)
# for negatives excluding NN
cornepnrna = corpas(rbind(callpas[['NE']],callpas[['PN']]),1,2,hredtrna,minnovs)
# for PN
corpnrna = corpas(callpas[['PN']],1,2,hredtrna,minnovs)
# for pos or negatives excluding PN
corpoppnennrna = corpas(rbind(callpas[['PO']],callpas[['PP']],callpas[['NE']],callpas[['NN']]),1,2,hredtrna,minnovs)

# ALL
# plot RNA Correlation distributions
svg(file='PSN_normal_rna_cordist.svg',width=4.5, height=4)
plot(density(as.numeric(corpprna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
lines(density(as.numeric(cornnrna[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
lines(density(as.numeric(corpnrna[,3]),na.rm=TRUE),col='blue4',lwd = 1.5)
lines(density(as.numeric(corponerna[,3]),na.rm=TRUE),col='green3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','NN','PN','NO'),col = c('purple','darkorange3','blue4','green3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# boxplot RNA Correlation distributions
svg(file='PSN_normal_rna_cordist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(corpprna[,3]),as.numeric(cornnrna[,3]),as.numeric(corpnrna[,3]),
            as.numeric(corponerna[,3]))
names(dali)=c('PP','NN','PN','NO')
boxplot(dali,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Normal RNA Cor',ylab="RNA Correlation",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

# PP
# plot comparison of other pos correlations to mutual PP correlations
svg(file='PSN_normal_rna_PPdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpopnrna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpprna[,3]),na.rm=TRUE),col='purple',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PN','PP'),col = c('gray2','purple'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# NN
# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_normal_rna_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnrna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(cornnrna[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NE_PN','NN'),col = c('gray2','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# PN
# plot comparison of other pos or neg correlations to mutual PN correlations
svg(file='PSN_normal_rna_PNdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpoppnennrna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpnrna[,3]),na.rm=TRUE),col='lightblue3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PP_NE_NN','PN'),col = c('gray2','lightblue3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# pairwise significance tests
ks.test(as.numeric(corpprna[,3]),as.numeric(corpopnrna[,3]))
ks.test(as.numeric(cornnrna[,3]),as.numeric(cornepnrna[,3]))
ks.test(as.numeric(corpnrna[,3]),as.numeric(corpoppnennrna[,3]))
wilcox.test(as.numeric(corpprna[,3]),as.numeric(corpopnrna[,3]))
wilcox.test(as.numeric(cornnrna[,3]),as.numeric(cornepnrna[,3]))
wilcox.test(as.numeric(corpnrna[,3]),as.numeric(corpoppnennrna[,3]))
ks.test(as.numeric(cornnrna[,3]),as.numeric(corpprna[,3]))
ks.test(as.numeric(cornnrna[,3]),as.numeric(corpnrna[,3]))
ks.test(as.numeric(cornnrna[,3]),as.numeric(corponerna[,3]))
wilcox.test(as.numeric(cornnrna[,3]),as.numeric(corpprna[,3]))
wilcox.test(as.numeric(cornnrna[,3]),as.numeric(corpnrna[,3]))
wilcox.test(as.numeric(cornnrna[,3]),as.numeric(corponerna[,3]))



# Cancer CORRELATIONS
# exclude self pairs
# for pairs excluding PP,PN,NN
corponernatu = corpas(rbind(callpas[['PO']],callpas[['NE']]),1,2,hredtrnatu,minnovs)
# for negatives excluding NN,PN
cornernatu = corpas(callpas[['NE']],1,2,hredtrnatu,minnovs)
# for PP
corpprnatu = corpas(callpas[['PP']],1,2,hredtrnatu,minnovs)
# for positives excluding PP
corpopnrnatu = corpas(rbind(callpas[['PO']],callpas[['PN']]),1,2,hredtrnatu,minnovs)
# for NN
cornnrnatu = corpas(callpas[['NN']],1,2,hredtrnatu,minnovs)
# for negatives excluding NN
cornepnrnatu = corpas(rbind(callpas[['NE']],callpas[['PN']]),1,2,hredtrnatu,minnovs)
# for PN
corpnrnatu = corpas(callpas[['PN']],1,2,hredtrnatu,minnovs)
# for pos or negatives excluding PN
corpoppnennrnatu = corpas(rbind(callpas[['PO']],callpas[['PP']],callpas[['NE']],callpas[['NN']]),1,2,hredtrnatu,minnovs)

# ALL
# plot RNA Correlation distributions
svg(file='PSN_Cancer_rna_cordist.svg',width=4.5, height=4)
plot(density(as.numeric(corpprnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
lines(density(as.numeric(cornnrnatu[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
lines(density(as.numeric(corpnrnatu[,3]),na.rm=TRUE),col='blue4',lwd = 1.5)
lines(density(as.numeric(corponernatu[,3]),na.rm=TRUE),col='green3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','NN','PN','NO'),col = c('purple','darkorange3','blue4','green3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# boxplot RNA Correlation distributions
svg(file='PSN_Cancer_rna_cordist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(corpprnatu[,3]),as.numeric(cornnrnatu[,3]),as.numeric(corpnrnatu[,3]),
            as.numeric(corponernatu[,3]))
names(dali)=c('PP','NN','PN','NO')
boxplot(dali,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Cancer RNA Cor',ylab="RNA Correlation",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

# PP
# plot comparison of other pos correlations to mutual PP correlations
svg(file='PSN_Cancer_rna_PPdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpopnrnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpprnatu[,3]),na.rm=TRUE),col='purple',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PN','PP'),col = c('gray2','purple'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# NN
# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_Cancer_rna_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnrnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(cornnrnatu[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NE_PN','NN'),col = c('gray2','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# PN
# plot comparison of other pos or neg correlations to mutual PN correlations
svg(file='PSN_Cancer_rna_PNdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpoppnennrnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpnrnatu[,3]),na.rm=TRUE),col='lightblue3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PP_NE_NN','PN'),col = c('gray2','lightblue3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# pairwise significance tests
ks.test(as.numeric(corpprnatu[,3]),as.numeric(corpopnrnatu[,3]))
ks.test(as.numeric(cornnrnatu[,3]),as.numeric(cornepnrnatu[,3]))
ks.test(as.numeric(corpnrnatu[,3]),as.numeric(corpoppnennrnatu[,3]))
wilcox.test(as.numeric(corpprnatu[,3]),as.numeric(corpopnrnatu[,3]))
wilcox.test(as.numeric(cornnrnatu[,3]),as.numeric(cornepnrnatu[,3]))
wilcox.test(as.numeric(corpnrnatu[,3]),as.numeric(corpoppnennrnatu[,3]))
ks.test(as.numeric(cornnrnatu[,3]),as.numeric(corpprnatu[,3]))
ks.test(as.numeric(cornnrnatu[,3]),as.numeric(corpnrnatu[,3]))
ks.test(as.numeric(cornnrnatu[,3]),as.numeric(corponernatu[,3]))
wilcox.test(as.numeric(cornnrnatu[,3]),as.numeric(corpprnatu[,3]))
wilcox.test(as.numeric(cornnrnatu[,3]),as.numeric(corpnrnatu[,3]))
wilcox.test(as.numeric(cornnrnatu[,3]),as.numeric(corponernatu[,3]))

wilcox.test(as.numeric(cornnrnatu[,3]),as.numeric(cornnrna[,3]),alternative = 'less')
wilcox.test(as.numeric(corpprnatu[,3]),as.numeric(corpprna[,3]),alternative = 'less')
wilcox.test(as.numeric(corpnrnatu[,3]),as.numeric(corpnrna[,3]),alternative = 'less')
wilcox.test(as.numeric(corponernatu[,3]),as.numeric(corponerna[,3]),alternative = 'less')

####################
# PROTEIN EXPRESSION
####################

hredtprot = read.delim("Reduced_Normal_Protein_exp.txt")
hredtprottu = read.delim("Reduced_Cancer_Protein_exp.txt")
# remove highly varied genes based on cov.
#protvarthr = 150
#hredtprot$Cov = hredtprot$Std.exp/hredtprot$Mean.exp*100
#hredtprottu$Cov = hredtprottu$Std.exp/hredtprottu$Mean.exp*100
#hredtprot = hredtprot[hredtprot$Std.exp<protvarthr,]
#hredtprottu = hredtprottu[hredtprottu$Std.exp<protvarthr,]
# there is no NA row in protein dataset
# use the same set of genes for normal and Cancer
commonps = intersect(unique(hredtprot$Gene.name),unique(hredtprottu$Gene.name))
hredtprot = hredtprot[hredtprot$Gene.name %in% commonps,]
hredtprottu = hredtprottu[hredtprottu$Gene.name %in% commonps,]


# NORMAL CORRELATIONS
# calculate correlations using spearman
# exclude self pairs
# for pairs excluding PP,PN,NN
corponeprot = corpas(rbind(callpas[['PO']],callpas[['NE']]),1,2,hredtprot,minnovs)
# for negatives excluding NN,PN
corneprot = corpas(callpas[['NE']],1,2,hredtprot,minnovs)
# for PP
corppprot = corpas(callpas[['PP']],1,2,hredtprot,minnovs)
# for positives excluding PP
corpopnprot = corpas(rbind(callpas[['PO']],callpas[['PN']]),1,2,hredtprot,minnovs)
# for NN
cornnprot = corpas(callpas[['NN']],1,2,hredtprot,minnovs)
# for negatives excluding NN
cornepnprot = corpas(rbind(callpas[['NE']],callpas[['PN']]),1,2,hredtprot,minnovs)
# for PN
corpnprot = corpas(callpas[['PN']],1,2,hredtprot,minnovs)
# for pos or negatives excluding PN
corpoppnennprot = corpas(rbind(callpas[['PO']],callpas[['PP']],callpas[['NE']],callpas[['NN']]),1,2,hredtprot,minnovs)

# ALL
# plot Protein Correlation distributions
svg(file='PSN_normal_prot_cordist.svg',width=4.5, height=4)
plot(density(as.numeric(corppprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
lines(density(as.numeric(cornnprot[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
lines(density(as.numeric(corpnprot[,3]),na.rm=TRUE),col='blue4',lwd = 1.5)
lines(density(as.numeric(corponeprot[,3]),na.rm=TRUE),col='green3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','NN','PN','NO'),col = c('purple','darkorange3','blue4','green3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# boxplot Protein Correlation distributions
svg(file='PSN_normal_prot_cordist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(corppprot[,3]),as.numeric(cornnprot[,3]),as.numeric(corpnprot[,3]),
            as.numeric(corponeprot[,3]))
names(dali)=c('PP','NN','PN','NO')
boxplot(dali,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Normal Protein Cor',ylab="Protein Correlation",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

# PP
# plot comparison of other pos correlations to mutual PP correlations
svg(file='PSN_normal_prot_PPdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpopnprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corppprot[,3]),na.rm=TRUE),col='purple',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PN','PP'),col = c('gray2','purple'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# NN
# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_normal_prot_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(cornnprot[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NE_PN','NN'),col = c('gray2','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PN
# plot comparison of other pos or neg correlations to mutual PN correlations
svg(file='PSN_normal_prot_PNdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpoppnennprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpnprot[,3]),na.rm=TRUE),col='lightblue3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PP_NE_NN','PN'),col = c('gray2','lightblue3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# pairwise significance tests
ks.test(as.numeric(corppprot[,3]),as.numeric(corpopnprot[,3]))
ks.test(as.numeric(cornnprot[,3]),as.numeric(cornepnprot[,3]))
ks.test(as.numeric(corpnprot[,3]),as.numeric(corpoppnennprot[,3]))
wilcox.test(as.numeric(corppprot[,3]),as.numeric(corpopnprot[,3]))
wilcox.test(as.numeric(cornnprot[,3]),as.numeric(cornepnprot[,3]),alternative='less')
wilcox.test(as.numeric(corpnprot[,3]),as.numeric(corpoppnennprot[,3]))

ks.test(as.numeric(cornnprot[,3]),as.numeric(corppprot[,3]))
ks.test(as.numeric(cornnprot[,3]),as.numeric(corpnprot[,3]))
ks.test(as.numeric(cornnprot[,3]),as.numeric(corponeprot[,3]))
wilcox.test(as.numeric(cornnprot[,3]),as.numeric(corppprot[,3]),alternative='less')
wilcox.test(as.numeric(cornnprot[,3]),as.numeric(corpnprot[,3]),alternative='less')
wilcox.test(as.numeric(cornnprot[,3]),as.numeric(corponeprot[,3]),alternative='less')



# Cancer CORRELATIONS
# calculate correlations using spearman (require at least 5 values per gene)
# exclude self pairs
# for pairs excluding PP,PN,NN
corponeprottu = corpas(rbind(callpas[['PO']],callpas[['NE']]),1,2,hredtprottu,minnovs)
# for negatives excluding NN,PN
corneprottu = corpas(callpas[['NE']],1,2,hredtprottu,minnovs)
# for PP
corppprottu = corpas(callpas[['PP']],1,2,hredtprottu,minnovs)
# for positives excluding PP
corpopnprottu = corpas(rbind(callpas[['PO']],callpas[['PN']]),1,2,hredtprottu,minnovs)
# for NN
cornnprottu = corpas(callpas[['NN']],1,2,hredtprottu,minnovs)
# for negatives excluding NN
cornepnprottu = corpas(rbind(callpas[['NE']],callpas[['PN']]),1,2,hredtprottu,minnovs)
# for PN
corpnprottu = corpas(callpas[['PN']],1,2,hredtprottu,minnovs)
# for pos or negatives excluding PN
corpoppnennprottu = corpas(rbind(callpas[['PO']],callpas[['PP']],callpas[['NE']],callpas[['NN']]),1,2,hredtprottu,minnovs)

# ALL
# plot Protein Correlation distributions
svg(file='PSN_Cancer_prot_cordist.svg',width=4.5, height=4)
plot(density(as.numeric(corppprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
lines(density(as.numeric(cornnprottu[,3]),na.rm=TRUE),col='darkorange3',lwd = 1.5)
lines(density(as.numeric(corpnprottu[,3]),na.rm=TRUE),col='blue4',lwd = 1.5)
lines(density(as.numeric(corponeprottu[,3]),na.rm=TRUE),col='green3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','NN','PN','NO'),col = c('purple','darkorange3','blue4','green3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# boxplot Protein Correlation distributions
svg(file='PSN_Cancer_prot_cordist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(corppprottu[,3]),as.numeric(cornnprottu[,3]),as.numeric(corpnprottu[,3]),
            as.numeric(corponeprottu[,3]))
names(dali)=c('PP','NN','PN','NO')
boxplot(dali,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Cancer Protein Cor',ylab="Protein Correlation",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

# PP
# plot comparison of other pos correlations to mutual PP correlations
svg(file='PSN_Cancer_prot_PPdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpopnprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corppprottu[,3]),na.rm=TRUE),col='purple',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PN','PP'),col = c('gray2','purple'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Pos Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# NN
# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_prot_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 2,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue')
lines(density(as.numeric(cornnprottu[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
lines(density(as.numeric(cornepnprot[,3]),na.rm=TRUE),col='blue3',lwd = 2,lty=2)
lines(density(as.numeric(cornnprot[,3]),na.rm=TRUE),col='darkorange3',lwd = 2,lty=2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,lty=c(1,1,2,2),
       legend=c('Cancer Neg','Cancer NN','Normal Neg','Normal NN'),col = c('blue','darkorange3','blue3','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_Cancer_prot_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 2,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(cornnprottu[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,
       legend=c('NE_PN','NN'),col = c('gray2','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# plot comparison of other neg correlations to mutual NN correlations
svg(file='PSN_normal prot_NNdist.svg',width=4.5, height=4)
plot(density(as.numeric(cornepnprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 2,lty=2,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue3')
lines(density(as.numeric(cornnprot[,3]),na.rm=TRUE),col='darkorange3',lwd = 2,lty=2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,lty=c(2,2),
       legend=c('Normal Neg','Normal NN'),col = c('blue3','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# plot ecdf mutual NN correlations
svg(file='PSN_cdf_prot_NNdist.svg',width=4.5, height=4)
plot(ecdf(as.numeric(cornepnprottu[,3])),ylim=c(0,1.3),lwd = 2,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue')
lines(ecdf(as.numeric(cornnprottu[,3])),col='darkorange3')
lines(ecdf(as.numeric(cornepnprot[,3])),col='blue3',lwd = 2,lty=2)
lines(ecdf(as.numeric(cornnprot[,3])),col='darkorange3',lty=2)
legend(x = -1,y = 1.2,bty = 'n',lwd = 1.5,lty=c(1,1,2,2),
       legend=c('Cancer Neg','Cancer NN','Normal Neg','Normal NN'),col = c('blue','darkorange','blue3','darkorange3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PN
# plot comparison of other pos or neg correlations to mutual PN correlations
svg(file='PSN_Cancer_prot_PNdist.svg',width=4.5, height=4)
plot(density(as.numeric(corpoppnennprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='gray2')
lines(density(as.numeric(corpnprottu[,3]),na.rm=TRUE),col='lightblue3',lwd = 1.5)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PO_PP_NE_NN','PN'),col = c('gray2','lightblue3'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# pairwise significance tests
ks.test(as.numeric(corppprottu[,3]),as.numeric(corpopnprottu[,3]))
ks.test(as.numeric(cornnprottu[,3]),as.numeric(cornepnprottu[,3]))
ks.test(as.numeric(corpnprottu[,3]),as.numeric(corpoppnennprottu[,3]))
wilcox.test(as.numeric(corppprottu[,3]),as.numeric(corpopnprottu[,3]))
wilcox.test(as.numeric(cornnprottu[,3]),as.numeric(cornepnprottu[,3]),alternative='less')
wilcox.test(as.numeric(corpnprottu[,3]),as.numeric(corpoppnennprottu[,3]))


ks.test(as.numeric(cornnprottu[,3]),as.numeric(corppprottu[,3]))
ks.test(as.numeric(cornnprottu[,3]),as.numeric(corpnprottu[,3]))
ks.test(as.numeric(cornnprottu[,3]),as.numeric(corponeprottu[,3]))
wilcox.test(as.numeric(cornnprottu[,3]),as.numeric(corppprottu[,3]),alternative='less')
wilcox.test(as.numeric(cornnprottu[,3]),as.numeric(corpnprottu[,3]),alternative='less')
wilcox.test(as.numeric(cornnprottu[,3]),as.numeric(corponeprottu[,3]),alternative='less')


wilcox.test(as.numeric(cornnprottu[,3]),as.numeric(cornnprot[,3]),alternative = 'less')
wilcox.test(as.numeric(corppprottu[,3]),as.numeric(corppprot[,3]),alternative = 'less')
wilcox.test(as.numeric(corpnprottu[,3]),as.numeric(corpnprot[,3]),alternative = 'less')
wilcox.test(as.numeric(corponeprottu[,3]),as.numeric(corponeprot[,3]),alternative = 'less')


####################################################
# random sample sets for similar size to mutual lists
####################################################
rno = 1000

#####
# PP
####
# for sets of positive pairs with similar size to PP 
pospps = ps[ps$mor==1,c(1,2)]

# normal RNA expression
rcoppvrna = c()
for (rii in 1:rno){
    rappv = pospps[sample(1:nrow(pospps),nrow(callpas[['PP']])),]
    randcors = corpas(rappv,1,2,hredtrna,minnovs)
    rcoppvrna = rbind(rcoppvrna,cbind(randcors,rep(rii,nrow(randcors))))}
rcoppvrna = as.data.frame(rcoppvrna)
names(rcoppvrna)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcoppvrna,file="PSN_randomcors_normal_rna_PP.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcoppvrna = read.delim('PSN_randomcors_normal_rna_PP.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcoppvrna[rcoppvrna$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpprna[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpprna[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_rna_PPdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpprna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
for (i in unique(rcoppvrna$rno)){
    lines(density(as.numeric(rcoppvrna[rcoppvrna$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcoppvrna[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpprna[,3]),na.rm=TRUE),col='purple',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','R'),col = c('purple','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer RNA expression
rcoppvrnatu = c()
for (rii in 1:rno){
    rappv = pospps[sample(1:nrow(pospps),nrow(callpas[['PP']])),]
    randcors = corpas(rappv,1,2,hredtrnatu,minnovs)
    rcoppvrnatu = rbind(rcoppvrnatu,cbind(randcors,rep(rii,nrow(randcors))))}
rcoppvrnatu = as.data.frame(rcoppvrnatu)
names(rcoppvrnatu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcoppvrnatu,file="PSN_randomcors_Cancer_rna_PP.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcoppvrnatu = read.delim('PSN_randomcors_Cancer_rna_PP.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcoppvrnatu[rcoppvrnatu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpprnatu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpprnatu[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_rna_PPdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpprnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
for (i in unique(rcoppvrnatu$rno)){
    lines(density(as.numeric(rcoppvrnatu[rcoppvrnatu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcoppvrnatu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpprnatu[,3]),na.rm=TRUE),col='purple',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','R'),col = c('purple','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# normal protein expression
rcoppvprot = c()
for (rii in 1:rno){
    rappv = pospps[sample(1:nrow(pospps),nrow(callpas[['PP']])),]
    randcors = corpas(rappv,1,2,hredtprot,minnovs)
    rcoppvprot = rbind(rcoppvprot,cbind(randcors,rep(rii,nrow(randcors))))}
rcoppvprot = as.data.frame(rcoppvprot)
names(rcoppvprot)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcoppvprot,file="PSN_randomcors_normal_prot_PP.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcoppvprot = read.delim('PSN_randomcors_normal_prot_PP.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcoppvprot[rcoppvprot$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corppprot[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corppprot[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_prot_PPdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corppprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
for (i in unique(rcoppvprot$rno)){
    lines(density(as.numeric(rcoppvprot[rcoppvprot$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcoppvprot[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corppprot[,3]),na.rm=TRUE),col='purple',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','R'),col = c('purple','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer protein expression
rcoppvprottu = c()
for (rii in 1:rno){
    rappv = pospps[sample(1:nrow(pospps),nrow(callpas[['PP']])),]
    randcors = corpas(rappv,1,2,hredtprottu,minnovs)
    rcoppvprottu = rbind(rcoppvprottu,cbind(randcors,rep(rii,nrow(randcors))))}
rcoppvprottu = as.data.frame(rcoppvprottu)
names(rcoppvprottu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcoppvprottu,file="PSN_randomcors_Cancer_prot_PP.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcoppvprottu = read.delim('PSN_randomcors_Cancer_prot_PP.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcoppvprottu[rcoppvprottu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corppprottu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corppprottu[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_prot_PPdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corppprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='purple')
for (i in unique(rcoppvprottu$rno)){
    lines(density(as.numeric(rcoppvprottu[rcoppvprottu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcoppvprottu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corppprottu[,3]),na.rm=TRUE),col='purple',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PP','R'),col = c('purple','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

####
# NN
###
# for sets of negative pairs with similar size to NN
negpps = ps[ps$mor==-1,c(1,2)]

# normal RNA expression
rconnvrna = c()
for (rii in 1:rno){
    rappv = negpps[sample(1:nrow(negpps),nrow(callpas[['NN']])),]
    randcors = corpas(rappv,1,2,hredtrna,minnovs)
    rconnvrna = rbind(rconnvrna,cbind(randcors,rep(rii,nrow(randcors))))}
rconnvrna = as.data.frame(rconnvrna)
names(rconnvrna)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rconnvrna,file="PSN_randomcors_normal_rna_NN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rconnvrna = read.delim('PSN_randomcors_normal_rna_NN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rconnvrna[rconnvrna$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(cornnrna[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(cornnrna[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_rna_NNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(cornnrna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='darkorange3')
for (i in unique(rconnvrna$rno)){
    lines(density(as.numeric(rconnvrna[rconnvrna$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rconnvrna[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(cornnrna[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','R'),col = c('darkorange3','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer RNA expression
rconnvrnatu = c()
for (rii in 1:rno){
    rappv = negpps[sample(1:nrow(negpps),nrow(callpas[['NN']])),]
    randcors = corpas(rappv,1,2,hredtrnatu,minnovs)
    rconnvrnatu = rbind(rconnvrnatu,cbind(randcors,rep(rii,nrow(randcors))))}
rconnvrnatu = as.data.frame(rconnvrnatu)
names(rconnvrnatu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rconnvrnatu,file="PSN_randomcors_Cancer_rna_NN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rconnvrnatu = read.delim('PSN_randomcors_Cancer_rna_NN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rconnvrnatu[rconnvrnatu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(cornnrnatu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(cornnrnatu[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_rna_NNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(cornnrnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='darkorange3')
for (i in unique(rconnvrnatu$rno)){
    lines(density(as.numeric(rconnvrnatu[rconnvrnatu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rconnvrnatu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(cornnrnatu[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','R'),col = c('darkorange3','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# normal protein expression
rconnvprot = c()
for (rii in 1:rno){
    rappv = negpps[sample(1:nrow(negpps),nrow(callpas[['NN']])),]
    randcors = corpas(rappv,1,2,hredtprot,minnovs)
    rconnvprot = rbind(rconnvprot,cbind(randcors,rep(rii,nrow(randcors))))}
rconnvprot = as.data.frame(rconnvprot)
names(rconnvprot)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rconnvprot,file="PSN_randomcors_normal_prot_NN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rconnvprot = read.delim('PSN_randomcors_normal_prot_NN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rconnvprot[rconnvprot$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(cornnprot[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(cornnprot[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_prot_NNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(cornnprot[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='darkorange3')
for (i in unique(rconnvprot$rno)){
    lines(density(as.numeric(rconnvprot[rconnvprot$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rconnvprot[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(cornnprot[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','R'),col = c('darkorange3','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer protein expression
rconnvprottu = c()
for (rii in 1:rno){
    rappv = negpps[sample(1:nrow(negpps),nrow(callpas[['NN']])),]
    randcors = corpas(rappv,1,2,hredtprottu,minnovs)
    rconnvprottu = rbind(rconnvprottu,cbind(randcors,rep(rii,nrow(randcors))))}
rconnvprottu = as.data.frame(rconnvprottu)
names(rconnvprottu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rconnvprottu,file="PSN_randomcors_Cancer_prot_NN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rconnvprottu = read.delim('PSN_randomcors_Cancer_prot_NN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rconnvprottu[rconnvprottu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(cornnprottu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(cornnprottu[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_prot_NNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(cornnprottu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='darkorange3')
for (i in unique(rconnvprottu$rno)){
    lines(density(as.numeric(rconnvprottu[rconnvprottu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rconnvprottu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(cornnprottu[,3]),na.rm=TRUE),col='darkorange3',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','R'),col = c('darkorange3','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

#####
# PN
####
# for sets of positive pairs with similar size to PN
posnegpps = ps[,c(1,2)]

# normal RNA expression
rcopnvrna = c()
for (rii in 1:rno){
    rappv = posnegpps[sample(1:nrow(posnegpps),nrow(callpas[['PN']])),]
    randcors = corpas(rappv,1,2,hredtrna,minnovs)
    rcopnvrna = rbind(rcopnvrna,cbind(randcors,rep(rii,nrow(randcors))))}
rcopnvrna = as.data.frame(rcopnvrna)
names(rcopnvrna)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcopnvrna,file="PSN_randomcors_normal_rna_PN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcopnvrna = read.delim('PSN_randomcors_normal_rna_PN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcopnvrna[rcopnvrna$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpnrna[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpnrna[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_rna_PNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpnrna[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
for (i in unique(rcopnvrna$rno)){
    lines(density(as.numeric(rcopnvrna[rcopnvrna$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcopnvrna[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpnrna[,3]),na.rm=TRUE),col='blue2',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PN','R'),col = c('blue2','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Normal RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer RNA expression
rcopnvrnatu = c()
for (rii in 1:rno){
    rappv = posnegpps[sample(1:nrow(posnegpps),nrow(callpas[['PN']])),]
    randcors = corpas(rappv,1,2,hredtrnatu,minnovs)
    rcopnvrnatu = rbind(rcopnvrnatu,cbind(randcors,rep(rii,nrow(randcors))))}
rcopnvrnatu = as.data.frame(rcopnvrnatu)
names(rcopnvrnatu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcopnvrnatu,file="PSN_randomcors_Cancer_rna_PN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcopnvrnatu = read.delim('PSN_randomcors_Cancer_rna_PN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcopnvrnatu[rcopnvrnatu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpnrnatu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpnrnatu[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_rna_PNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpnrnatu[,3]),na.rm=TRUE),ylim=c(0,2.5),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
for (i in unique(rcopnvrnatu$rno)){
    lines(density(as.numeric(rcopnvrnatu[rcopnvrnatu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcopnvrnatu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpnrnatu[,3]),na.rm=TRUE),col='blue2',lwd = 2)
legend(x = -1,y = 2.5,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PN','R'),col = c('blue2','gray'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Cancer RNA Cor',ylab ="Frequency Density",
      xlab="RNA Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# normal protein expression
rcopnvprot = c()
for (rii in 1:rno){
    rappv = posnegpps[sample(1:nrow(posnegpps),nrow(callpas[['PN']])),]
    randcors = corpas(rappv,1,2,hredtprot,minnovs)
    rcopnvprot = rbind(rcopnvprot,cbind(randcors,rep(rii,nrow(randcors))))}
rcopnvprot = as.data.frame(rcopnvprot)
names(rcopnvprot)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcopnvprot,file="PSN_randomcors_normal_prot_PN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcopnvprot = read.delim('PSN_randomcors_normal_prot_PN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcopnvprot[rcopnvprot$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpnprot[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpnprot[,3]))))+1)/rno
print(pvalue)

# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_normal_prot_PNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpnprot[,3]),na.rm=TRUE),ylim=c(0,3),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
for (i in unique(rcopnvprot$rno)){
    lines(density(as.numeric(rcopnvprot[rcopnvprot$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcopnvprot[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpnprot[,3]),na.rm=TRUE),col='blue2',lwd = 2)
legend(x = -1,y = 3,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PN','R'),col = c('blue2','gray2'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Normal Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# Cancer protein expression
rcopnvprottu = c()
for (rii in 1:rno){
    rappv = posnegpps[sample(1:nrow(posnegpps),nrow(callpas[['PN']])),]
    randcors = corpas(rappv,1,2,hredtprottu,minnovs)
    rcopnvprottu = rbind(rcopnvprottu,cbind(randcors,rep(rii,nrow(randcors))))}
rcopnvprottu = as.data.frame(rcopnvprottu)
names(rcopnvprottu)=c('G1','G2','rho','p.value','FDR','rno')
write.table(rcopnvprottu,file="PSN_randomcors_Cancer_prot_PN.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

rcopnvprottu = read.delim('PSN_randomcors_Cancer_prot_PN.txt')
rhos = c()
for (ri in 1:rno){
    rhos = c(rhos,median(rcopnvprottu[rcopnvprottu$rno==ri,]$rho))}
pvalue = (length(which(rhos<=median(as.numeric(corpnprottu[,3]))))+1)/rno
print(pvalue)
pvalue = (length(which(rhos>=median(as.numeric(corpnprottu[,3]))))+1)/rno
print(pvalue)


# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_Cancer_prot_PNdist_randsets.svg',width=4.5, height=4)
plot(density(as.numeric(corpnprottu[,3]),na.rm=TRUE),ylim=c(0,3),lwd = 1.5,
     xlim=c(-1,1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
for (i in unique(rcopnvprottu$rno)){
    lines(density(as.numeric(rcopnvprottu[rcopnvprottu$rno==i,][,3]),na.rm=TRUE),col='gray',lwd = 1.5)}
lines(density(as.numeric(rcopnvprottu[,3]),na.rm=TRUE),col='gray2',lwd = 2)
lines(density(as.numeric(corpnprottu[,3]),na.rm=TRUE),col='blue2',lwd = 2)
legend(x = -1,y = 3,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('PN','R'),col = c('blue2','gray2'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN Cancer Protein Cor',ylab ="Frequency Density",
      xlab="Protein Correlation",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()




############################
# differential correlations 
############################

# Cancer - normal expression of each pair
# if Cancer dataset also has the same pair 
# some mismatches due to NA cor values which will be ignored

# PP RNA
corpprnadif = c()
for (ci in 1:nrow(corpprna)){
    norvs = corpprna[ci,]
    # if Cancer dataset also has the same pair 
    tumvs = corpprnatu[corpprnatu[,1]==norvs[1]&corpprnatu[,2]==norvs[2],]
    if (length(tumvs)>0){
        corpprnadif = rbind(corpprnadif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corpprnadif = cbind(corpprnadif,as.numeric(corpprnadif[,6])-as.numeric(corpprnadif[,3]))
corpprnadif = as.data.frame(corpprnadif)
colnames(corpprnadif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corpprnadif,file='PSN_PP_RNA_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# NN RNA
cornnrnadif = c()
for (ci in 1:nrow(cornnrna)){
  norvs = cornnrna[ci,]
  # if Cancer dataset also has the same pair 
  tumvs = cornnrnatu[cornnrnatu[,1]==norvs[1]&cornnrnatu[,2]==norvs[2],]
  if (length(tumvs)>0){
    cornnrnadif = rbind(cornnrnadif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
cornnrnadif = cbind(cornnrnadif,as.numeric(cornnrnadif[,6])-as.numeric(cornnrnadif[,3]))
cornnrnadif = as.data.frame(cornnrnadif)
colnames(cornnrnadif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(cornnrnadif,file='PSN_NN_RNA_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# PN RNA
corpnrnadif = c()
for (ci in 1:nrow(corpnrna)){
  norvs = corpnrna[ci,]
  # if Cancer dataset also has the same pair 
  tumvs = corpnrnatu[corpnrnatu[,1]==norvs[1]&corpnrnatu[,2]==norvs[2],]
  if (length(tumvs)>0){
    corpnrnadif = rbind(corpnrnadif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corpnrnadif = cbind(corpnrnadif,as.numeric(corpnrnadif[,6])-as.numeric(corpnrnadif[,3]))
corpnrnadif = as.data.frame(corpnrnadif)
colnames(corpnrnadif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corpnrnadif,file='PSN_PN_RNA_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)


# NO RNA
corponernadif = c()
for (ci in 1:nrow(corponerna)){
    norvs = corponerna[ci,]
    # if Cancer dataset also has the same pair 
    tumvs = corponernatu[corponernatu[,1]==norvs[1]&corponernatu[,2]==norvs[2],]
    if (length(tumvs)>0){
        corponernadif = rbind(corponernadif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corponernadif = cbind(corponernadif,as.numeric(corponernadif[,6])-as.numeric(corponernadif[,3]))
corponernadif = as.data.frame(corponernadif)
colnames(corponernadif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                             'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corponernadif,file='PSN_NO_RNA_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# boxplot RNA Correlation distributions
svg(file='PSN_RNA_diffcors_boxplot.svg',width=4.5, height=4)
diffa = list(as.numeric(corpprnadif$Diff_cor),as.numeric(cornnrnadif$Diff_cor),
             as.numeric(corpnrnadif$Diff_cor),as.numeric(corponernadif$Diff_cor))
names(diffa)=c('PP','NN','PN','NO')
boxplot(diffa,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN RNA Diff Cor',ylab="Cancer-Normal RNA Cor",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

ks.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corpprnadif$Diff_cor))
ks.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corpnrnadif$Diff_cor))
ks.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corponernadif$Diff_cor))

wilcox.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corpprnadif$Diff_cor),alternative ='less')
wilcox.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corpnrnadif$Diff_cor),alternative ='less')
wilcox.test(as.numeric(as.numeric(cornnrnadif$Diff_cor)),as.numeric(corponernadif$Diff_cor),alternative ='less')



# PP Protein
corppprotdif = c()
for (ci in 1:nrow(corppprot)){
  norvs = corppprot[ci,]
  # if Cancer dataset also has the same pair 
  tumvs = corppprottu[corppprottu[,1]==norvs[1]&corppprottu[,2]==norvs[2],]
  if (length(tumvs)>0){
    corppprotdif = rbind(corppprotdif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corppprotdif = cbind(corppprotdif,as.numeric(corppprotdif[,6])-as.numeric(corppprotdif[,3]))
corppprotdif = as.data.frame(corppprotdif)
colnames(corppprotdif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corppprotdif,file='PSN_PP_Protein_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# NN Protein
cornnprotdif = c()
for (ci in 1:nrow(cornnprot)){
  norvs = cornnprot[ci,]
  # if Cancer dataset also has the same pair 
  tumvs = cornnprottu[cornnprottu[,1]==norvs[1]&cornnprottu[,2]==norvs[2],]
  if (length(tumvs)>0){
    cornnprotdif = rbind(cornnprotdif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
cornnprotdif = cbind(cornnprotdif,as.numeric(cornnprotdif[,6])-as.numeric(cornnprotdif[,3]))
cornnprotdif = as.data.frame(cornnprotdif)
colnames(cornnprotdif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(cornnprotdif,file='PSN_NN_Protein_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# PN Protein
corpnprotdif = c()
for (ci in 1:nrow(corpnprot)){
  norvs = corpnprot[ci,]
  # if Cancer dataset also has the same pair 
  tumvs = corpnprottu[corpnprottu[,1]==norvs[1]&corpnprottu[,2]==norvs[2],]
  if (length(tumvs)>0){
    corpnprotdif = rbind(corpnprotdif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corpnprotdif = cbind(corpnprotdif,as.numeric(corpnprotdif[,6])-as.numeric(corpnprotdif[,3]))
corpnprotdif = as.data.frame(corpnprotdif)
colnames(corpnprotdif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                          'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corpnprotdif,file='PSN_PN_Protein_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# NO Protein
corponeprotdif = c()
for (ci in 1:nrow(corponeprot)){
    norvs = corponeprot[ci,]
    # if Cancer dataset also has the same pair 
    tumvs = corponeprottu[corponeprottu[,1]==norvs[1]&corponeprottu[,2]==norvs[2],]
    if (length(tumvs)>0){
        corponeprotdif = rbind(corponeprotdif,c(norvs,tumvs[3:5]))}}
# add cor diffeence as last column
corponeprotdif = cbind(corponeprotdif,as.numeric(corponeprotdif[,6])-as.numeric(corponeprotdif[,3]))
corponeprotdif = as.data.frame(corponeprotdif)
colnames(corponeprotdif) = c('Gene1','Gene2','Normal_Cor','Normal_pvalue','Normal_fdr',
                           'Cancer_Cor','Cancer_pvalue','Cancer_fdr','Diff_cor')
write.table(corponeprotdif,file='PSN_NO_Protein_NormalvsCancer_Cor_Values.txt',
            sep='\t',row.names = FALSE)

# boxplot Protein Correlation distributions
svg(file='PSN_Protein_diffcors_boxplot.svg',width=4.5, height=4)
diffa = list(as.numeric(corppprotdif$Diff_cor),as.numeric(cornnprotdif$Diff_cor),
             as.numeric(corpnprotdif$Diff_cor),as.numeric(corponeprotdif$Diff_cor))
names(diffa)=c('PP','NN','PN','NO')
boxplot(diffa,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Protein Diff Cor',ylab="Cancer-Normal Protein Cor",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE)
dev.off()

ks.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corppprotdif$Diff_cor))
ks.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corpnprotdif$Diff_cor))
ks.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corponeprotdif$Diff_cor))

wilcox.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corppprotdif$Diff_cor),alternative='greater')
wilcox.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corpnprotdif$Diff_cor),alternative='greater')
wilcox.test(as.numeric(as.numeric(cornnprotdif$Diff_cor)),as.numeric(corponeprotdif$Diff_cor),alternative='greater')


# write a table for all the directed pairs in the parent network
corallnet = c()
ckd = data.frame()
ckd = rbind(ckd,c('none','none'))
names(ckd) = c('A','B')
caPP = as.data.frame(callpas[['PP']])
caPN = as.data.frame(callpas[['PN']])
caNN = as.data.frame(callpas[['NN']])
trowrnaPP = c()
trowrnaPN = c()
trowrnaNN = c()
for (pi in 1:nrow(ps)){
  psp1 = as.character(ps[pi,2])
  psp2 = as.character(ps[pi,4])
  # continue if the pair(both directions) was not checked before
  if (nrow(ckd[ckd$A==psp1&ckd$B==psp2,])==0){
    ckd = rbind(ckd,c(psp1,psp2))
    ckd = rbind(ckd,c(psp2,psp1))
    trowrnaPP = c()
    trowrnaPN = c()
    trowrnaNN = c()
    # PP
    if (nrow(caPP[caPP$V1==psp1&caPP$V2==psp2,])>0|
        nrow(caPP[caPP$V1==psp2&caPP$V2==psp1,])>0){
      if (nrow(corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,])>0 &
          nrow(corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,])>0){
        trowrnaPP = c(psp1,psp2,'PSN','PP',
                      corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,3:9],
                      corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,3:9])}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,])>0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,])>0){
        trowrnaPP = c(psp2,psp1,'PSN','PP',
                      corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,3:9],
                      corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,3:9])}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,])>0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,])==0){
        trowrnaPP = c(psp1,psp2,'PSN','PP',
                      corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,])==0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,])>0){
        trowrnaPP = c(psp1,psp2,'PSN','PP',
                      'NA','NA','NA','NA','NA','NA','NA',
                      corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,3:9])}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,])>0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,])==0){
        trowrnaPP = c(psp2,psp1,'PSN','PP',
                      corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,])==0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,])>0){
        trowrnaPP = c(psp2,psp1,'PSN','PP',
                      'NA','NA','NA','NA','NA','NA','NA',
                      corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,3:9])}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp2&corpprnadif$Gene2==psp1,])==0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp2&corppprotdif$Gene2==psp1,])==0){
        trowrnaPP = c(psp2,psp1,'PSN','PP',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpprnadif[corpprnadif$Gene1==psp1&corpprnadif$Gene2==psp2,])==0 &
               nrow(corppprotdif[corppprotdif$Gene1==psp1&corppprotdif$Gene2==psp2,])==0){
        trowrnaPP = c(psp1,psp2,'PSN','PP',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      corallnet = rbind(corallnet,trowrnaPP)}
    
    #PN p1 p2 (direction is pos neg for PN)
    if (nrow(caPN[caPN$V1==psp1&caPN$V2==psp2,])>0|
        nrow(caPN[caPN$V1==psp2&caPN$V2==psp1,])>0){
      if (nrow(corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,])>0 &
          nrow(corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,])>0){
        trowrnaPN = c(psp1,psp2,'PSN','PN',
                      corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,3:9],
                      corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,3:9])}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,])>0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,])>0){
        trowrnaPN = c(psp2,psp1,'PSN','PN',
                      corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,3:9],
                      corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,3:9])}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,])>0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,])==0){
        trowrnaPN = c(psp1,psp2,'PSN','PN',
                      corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,])==0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,])>0){
        trowrnaPN = c(psp1,psp2,'PSN','PN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,3:9])}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,])>0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,])==0){
        trowrnaPN = c(psp2,psp1,'PSN','PN',
                      corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,])==0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,])>0){
        trowrnaPN = c(psp2,psp1,'PSN','PN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,3:9])}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp2&corpnrnadif$Gene2==psp1,])==0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp2&corpnprotdif$Gene2==psp1,])==0){
        trowrnaPN = c(psp2,psp1,'PSN','PN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(corpnrnadif[corpnrnadif$Gene1==psp1&corpnrnadif$Gene2==psp2,])==0 &
               nrow(corpnprotdif[corpnprotdif$Gene1==psp1&corpnprotdif$Gene2==psp2,])==0){
        trowrnaPN = c(psp1,psp2,'PSN','PN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      corallnet = rbind(corallnet,trowrnaPN)}
    
    
    #NN
    if (nrow(caNN[caNN$V1==psp1&caNN$V2==psp2,])>0|
        nrow(caNN[caNN$V1==psp2&caNN$V2==psp1,])>0){
      if (nrow(cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,])>0 &
          nrow(cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,])>0){
        trowrnaNN = c(psp1,psp2,'PSN','NN',
                      cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,3:9],
                      cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,3:9])}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,])>0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,])>0){
        trowrnaNN = c(psp2,psp1,'PSN','NN',
                      cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,3:9],
                      cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,3:9])}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,])>0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,])==0){
        trowrnaNN = c(psp1,psp2,'PSN','NN',
                      cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,])==0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,])>0){
        trowrnaNN = c(psp1,psp2,'PSN','NN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,3:9])}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,])>0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,])==0){
        trowrnaNN = c(psp2,psp1,'PSN','NN',
                      cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,3:9],
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,])==0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,])>0){
        trowrnaNN = c(psp2,psp1,'PSN','NN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,3:9])}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp2&cornnrnadif$Gene2==psp1,])==0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp2&cornnprotdif$Gene2==psp1,])==0){
        trowrnaNN = c(psp2,psp1,'PSN','NN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      else if (nrow(cornnrnadif[cornnrnadif$Gene1==psp1&cornnrnadif$Gene2==psp2,])==0 &
               nrow(cornnprotdif[cornnprotdif$Gene1==psp1&cornnprotdif$Gene2==psp2,])==0){
        trowrnaNN = c(psp1,psp2,'PSN','NN',
                      'NA','NA','NA','NA','NA','NA','NA',
                      'NA','NA','NA','NA','NA','NA','NA')}
      corallnet = rbind(corallnet,trowrnaNN)}
    
    
    if (length(trowrnaPP)==0&length(trowrnaPN)==0&length(trowrnaNN)==0){
      corallnet = rbind(corallnet,c(psp1,psp2,'PSN','NO',
                                    'NA','NA','NA','NA','NA','NA','NA',
                                    'NA','NA','NA','NA','NA','NA','NA'))}}}
corallnet = as.data.frame(corallnet)
colnames(corallnet) = c('Gene1','Gene2','Network','Mutual_type',
                        'Normal_RNA_Cor','Normal_RNA_pvalue','Normal_RNA_FDR',
                        'Cancer_RNA_Cor','Cancer_RNA_pvalue','Cancer_RNA_FDR','Diff_RNA_cor',
                        'Normal_Protein_Cor','Normal_Protein_pvalue','Normal_Protein_FDR',
                        'Cancer_Protein_Cor','Cancer_Protein_pvalue','Cancer_Protein_FDR','Diff_Protein_cor')

corallnet = apply(corallnet,2,as.character)
write.table(corallnet,file='PSN_All_Cor_Values.txt',sep='\t',row.names = FALSE,quote=FALSE)

corna = corallnet[is.na(corallnet$Diff_RNA_cor)==FALSE,]
cornann = corna[corna$Mutual_type=='NN',]
cornapn = corna[corna$Mutual_type=='PN',]
cornapp = corna[corna$Mutual_type=='PP',]

# plot diff rna cor
svg(file='PSN_difRNA_dist.svg',width=4.5, height=4)
plot(density(as.numeric(cornann$Diff_RNA_cor),na.rm=TRUE),ylim=c(0,2.3),lwd = 1.5,
     xlim=c(-1.3,1.3),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
lines(density(as.numeric(cornapn$Diff_RNA_cor),na.rm=TRUE),col='orange2',lwd = 2)
lines(density(as.numeric(cornapp$Diff_RNA_cor),na.rm=TRUE),col='red2',lwd = 2)
legend(x = -1,y = 2.3,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','PN','PP'),col = c('blue2','orange2','red2'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Dif RNA Cor',ylab ="Frequency Density",
      xlab="Differential RNA Cor",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

#corallnet=read.delim('PSN_All_Cor_Values.txt')
coprot = corallnet[is.na(corallnet$Diff_Protein_cor)==FALSE,]
coprotnn = coprot[coprot$Mutual_type=='NN',]
coprotpn = coprot[coprot$Mutual_type=='PN',]
coprotpp = coprot[coprot$Mutual_type=='PP',]


# plot random set correlation dists vs real set mutual correlations
svg(file='PSN_difprot_dist.svg',width=4.5, height=4)
plot(density(as.numeric(coprotnn$Diff_Protein_cor),na.rm=TRUE),ylim=c(0,2.3),lwd = 1.5,
     xlim=c(-1.3,1.3),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='blue2')
lines(density(as.numeric(coprotpn$Diff_Protein_cor),na.rm=TRUE),col='orange2',lwd = 2)
lines(density(as.numeric(coprotpp$Diff_Protein_cor),na.rm=TRUE),col='red2',lwd = 2)
legend(x = -1,y = 2.3,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('NN','PN','PP'),col = c('blue2','orange2','red2'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN Dif Protein Cor',ylab ="Frequency Density",
      xlab="Differential Protein Cor",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

#######################
# network calculations
#######################

# interactions among PP PN NN
cmts = intpas(callpas)
# unique lists of PP PN NN genes
capp = unique(c(callpas[['PP']][,1],callpas[['PP']][,2]))
capn = unique(c(callpas[['PN']][,1],callpas[['PN']][,2]))
cann = unique(c(callpas[['NN']][,1],callpas[['NN']][,2]))
# unique list of other genes (all genes ecluding PP PN NN)
psls = unique(c(ps$V2,ps$V4))
otherpss = setdiff(psls,unique(c(capp,capn,cann)))
# degree values of PP PN NN genes and others
cdgpp = degpas(ps,capp,2,4,5)
cdgpn = degpas(ps,capn,2,4,5)
cdgnn = degpas(ps,cann,2,4,5)
cdgos = degpas(ps,otherpss,2,4,5)
# limit otherpss to those with min 2 connections
indsmin2 = which(as.numeric(cdgos$total_deg)>=2)
# average stats of real vals
realcvs = c(nrow(callpas[['PP']]),nrow(callpas[['NN']]),nrow(callpas[['PN']]),
      nrow(callpas[['PO']]),nrow(callpas[['NE']]),
      sum(as.numeric(cmts[cmts$Int_type=='PP_PN',]$No)),
      sum(as.numeric(cmts[cmts$Int_type=='PP_NN',]$No)),
      sum(as.numeric(cmts[cmts$Int_type=='PN_NN',]$No)),
      sum(as.numeric(cmts[cmts$Int_type=='PP_PP',]$No)),
      sum(as.numeric(cmts[cmts$Int_type=='PN_PN',]$No)),
      sum(as.numeric(cmts[cmts$Int_type=='NN_NN',]$No)),
      median(as.numeric(cdgpp$total_deg)),median(as.numeric(cdgpp$in_deg)),
      median(as.numeric(cdgpp$out_deg)),median(as.numeric(cdgpp$pos_deg)),
      median(as.numeric(cdgpp$neg_deg)),median(as.numeric(cdgpp$inpos_deg)),
      median(as.numeric(cdgpp$inneg_deg)),median(as.numeric(cdgpp$outpos_deg)),
      median(as.numeric(cdgpp$outneg_deg)),
      median(as.numeric(cdgpn$total_deg)),median(as.numeric(cdgpn$in_deg)),
      median(as.numeric(cdgpn$out_deg)),median(as.numeric(cdgpn$pos_deg)),
      median(as.numeric(cdgpn$neg_deg)),median(as.numeric(cdgpn$inpos_deg)),
      median(as.numeric(cdgpn$inneg_deg)),median(as.numeric(cdgpn$outpos_deg)),
      median(as.numeric(cdgpn$outneg_deg)),
      median(as.numeric(cdgnn$total_deg)),median(as.numeric(cdgnn$in_deg)),
      median(as.numeric(cdgnn$out_deg)),median(as.numeric(cdgnn$pos_deg)),
      median(as.numeric(cdgnn$neg_deg)),median(as.numeric(cdgnn$inpos_deg)),
      median(as.numeric(cdgnn$inneg_deg)),median(as.numeric(cdgnn$outpos_deg)),
      median(as.numeric(cdgnn$outneg_deg)))


# boxplot of degree distributions (nonzero values in log scale)

# total deg
svg(file='PSN_degdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$total_deg),as.numeric(cdgnn$total_deg),
            as.numeric(cdgpn$total_deg),as.numeric(cdgos$total_deg)[indsmin2])
names(dali)=c('PP','NN','PN','NO')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
dali$NO=dali$NO[dali$NO!=0]
boxplot(dali,col=c('purple','darkorange3','blue4','green3'),
        main ='PSN Distribution',ylab="Total Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$total_deg),as.numeric(cdgnn$total_deg),alternative='greater')
wilcox.test(as.numeric(cdgpn$total_deg),as.numeric(cdgnn$total_deg),alternative='greater')
wilcox.test(as.numeric(cdgpp$total_deg),as.numeric(cdgos$total_deg)[indsmin2],alternative='greater')
wilcox.test(as.numeric(cdgpn$total_deg),as.numeric(cdgos$total_deg)[indsmin2],alternative='greater')
wilcox.test(as.numeric(cdgnn$total_deg),as.numeric(cdgos$total_deg)[indsmin2],alternative='greater')

# in deg
svg(file='PSN_indegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$in_deg),as.numeric(cdgnn$in_deg),
            as.numeric(cdgpn$in_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="In Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$in_deg),as.numeric(cdgnn$in_deg),alternative='greater')
wilcox.test(as.numeric(cdgpn$in_deg),as.numeric(cdgnn$in_deg),alternative='greater')

# out deg
svg(file='PSN_outdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$out_deg),as.numeric(cdgnn$out_deg),
            as.numeric(cdgpn$out_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="Out Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$out_deg),as.numeric(cdgnn$out_deg),alternative='less')
wilcox.test(as.numeric(cdgpn$out_deg),as.numeric(cdgnn$out_deg),alternative='less')

# pos deg
svg(file='PSN_posdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$pos_deg),as.numeric(cdgnn$pos_deg),
            as.numeric(cdgpn$pos_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="Pos Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$pos_deg),as.numeric(cdgnn$pos_deg),alternative='greater')
wilcox.test(as.numeric(cdgpn$pos_deg),as.numeric(cdgnn$pos_deg),alternative='greater')

# neg deg
svg(file='PSN_negdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$neg_deg),as.numeric(cdgnn$neg_deg),
            as.numeric(cdgpn$neg_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="Neg Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$neg_deg),as.numeric(cdgnn$neg_deg),alternative='less')
wilcox.test(as.numeric(cdgpn$neg_deg),as.numeric(cdgnn$neg_deg),alternative='less')

# in pos deg
svg(file='PSN_inposdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$inpos_deg),as.numeric(cdgnn$inpos_deg),
            as.numeric(cdgpn$inpos_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="In Pos Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$inpos_deg),as.numeric(cdgnn$inpos_deg),alternative='greater')
wilcox.test(as.numeric(cdgpn$inpos_deg),as.numeric(cdgnn$inpos_deg),alternative='greater')

# in neg deg
svg(file='PSN_innegdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$inneg_deg),as.numeric(cdgnn$inneg_deg),
            as.numeric(cdgpn$inneg_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="In Neg Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$inneg_deg),as.numeric(cdgnn$inneg_deg),alternative='less')
wilcox.test(as.numeric(cdgpn$inneg_deg),as.numeric(cdgnn$inneg_deg),alternative='less')

# out pos deg
svg(file='PSN_outposdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$outpos_deg),as.numeric(cdgnn$outpos_deg),
            as.numeric(cdgpn$outpos_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="Out Pos Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$outpos_deg),as.numeric(cdgnn$outpos_deg),alternative='less')
wilcox.test(as.numeric(cdgpn$outpos_deg),as.numeric(cdgnn$outpos_deg),alternative='less')

# out neg deg
svg(file='PSN_outnegdegdist_boxplot.svg',width=4.5, height=4)
dali = list(as.numeric(cdgpp$outneg_deg),as.numeric(cdgnn$outneg_deg),
            as.numeric(cdgpn$outneg_deg))
names(dali)=c('PP','NN','PN')
dali$PP=dali$PP[dali$PP!=0]
dali$PN=dali$PN[dali$PN!=0]
dali$NN=dali$NN[dali$NN!=0]
boxplot(dali,col=c('purple','darkorange3','blue4'),
        main ='PSN Distribution',ylab="Out Neg Degree",font=2,
        cex.lab=1.2,cex.axis=1.2,cex.main=1.3,notch=TRUE,log='y')
dev.off()

wilcox.test(as.numeric(cdgpp$outneg_deg),as.numeric(cdgnn$outneg_deg),alternative='less')
wilcox.test(as.numeric(cdgpn$outneg_deg),as.numeric(cdgnn$outneg_deg),alternative='less')

######################################
# random sampling network calculations
######################################
# randomize network and calculate mutual pair numbers
# as well as network properties
time1 = Sys.time()
rno = 1000
rpall <- foreach(ri=1:rno,.combine=rbind) %dopar% {
    # randomize only pos neg (shuffle) keeping int pairs same
    rps = cbind(ps[,1:4],V5=sample(ps$V5))
    rallpas = allpas(rps)
    rmts = intpas(rallpas)
    rapp = unique(c(rallpas[['PP']][,1],rallpas[['PP']][,2]))
    rapn = unique(c(rallpas[['PN']][,1],rallpas[['PN']][,2]))
    rann = unique(c(rallpas[['NN']][,1],rallpas[['NN']][,2]))
    rdgpp = degpas(rps,rapp,2,4,5)
    rdgpn = degpas(rps,rapn,2,4,5)
    rdgnn = degpas(rps,rann,2,4,5)
    cbind(nrow(rallpas[['PP']]),nrow(rallpas[['NN']]),nrow(rallpas[['PN']]),
          nrow(rallpas[['PO']]),nrow(rallpas[['NE']]),
          sum(as.numeric(rmts[rmts$Int_type=='PP_PN',]$No)),
          sum(as.numeric(rmts[rmts$Int_type=='PP_NN',]$No)),
          sum(as.numeric(rmts[rmts$Int_type=='PN_NN',]$No)),
          sum(as.numeric(rmts[rmts$Int_type=='PP_PP',]$No)),
          sum(as.numeric(rmts[rmts$Int_type=='PN_PN',]$No)),
          sum(as.numeric(rmts[rmts$Int_type=='NN_NN',]$No)),
          median(as.numeric(rdgpp$total_deg)),median(as.numeric(rdgpp$in_deg)),
          median(as.numeric(rdgpp$out_deg)),median(as.numeric(rdgpp$pos_deg)),
          median(as.numeric(rdgpp$neg_deg)),median(as.numeric(rdgpp$inpos_deg)),
          median(as.numeric(rdgpp$inneg_deg)),median(as.numeric(rdgpp$outpos_deg)),
          median(as.numeric(rdgpp$outneg_deg)),
          median(as.numeric(rdgpn$total_deg)),median(as.numeric(rdgpn$in_deg)),
          median(as.numeric(rdgpn$out_deg)),median(as.numeric(rdgpn$pos_deg)),
          median(as.numeric(rdgpn$neg_deg)),median(as.numeric(rdgpn$inpos_deg)),
          median(as.numeric(rdgpn$inneg_deg)),median(as.numeric(rdgpn$outpos_deg)),
          median(as.numeric(rdgpn$outneg_deg)),
          median(as.numeric(rdgnn$total_deg)),median(as.numeric(rdgnn$in_deg)),
          median(as.numeric(rdgnn$out_deg)),median(as.numeric(rdgnn$pos_deg)),
          median(as.numeric(rdgnn$neg_deg)),median(as.numeric(rdgnn$inpos_deg)),
          median(as.numeric(rdgnn$inneg_deg)),median(as.numeric(rdgnn$outpos_deg)),
          median(as.numeric(rdgnn$outneg_deg)))}
time2 = Sys.time()
print(time2-time1)
# add first row as real values and rest of the rows as random
rpall = rbind(realcvs,rpall)
rpall = as.data.frame(rpall)
colnames(rpall)=c('PP_no','NN_no','PN_no','PO_no','NE_no','PP_PN_sum',
                  'PP_NN_sum','PN_NN_sum','PP_PP_sum','PN_PN_sum','NN_NN_sum',
                  'PP_totdeg','PP_indeg','PP_outdeg','PP_posdeg','PP_negdeg',
                  'PP_inposdeg','PP_innegdeg','PP_outposdeg','PP_outnegdeg',
                  'PN_totdeg','PN_indeg','PN_outdeg','PN_posdeg','PN_negdeg',
                  'PN_inposdeg','PN_innegdeg','PN_outposdeg','PN_outnegdeg',
                  'NN_totdeg','NN_indeg','NN_outdeg','NN_posdeg','NN_negdeg',
                  'NN_inposdeg','NN_innegdeg','NN_outposdeg','NN_outnegdeg')
write.table(rpall,file="PSN_randomvals.txt",sep='\t',
            row.names=FALSE,quote=FALSE)

# plot pairs numbers random vs real
svg(file='PSN_randoms.svg',width=4.5, height=4)
plotb = barplot(c(mean(rpall$PP_no[-1]),rpall$PP_no[1],
                  mean(rpall$NN_no[-1]),rpall$NN_no[1],
                  mean(rpall$PN_no[-1]),rpall$PN_no[1]),
                names.arg=c('PP R','PP','NN R','NN','PN R','PN'),
                ylim=c(0,1200),main='PSN Real vs Random Size',
                ylab='Number of Pairs',col=c('gray','darkred'))
arrows(x0 = plotb,                      
       y0 = c(mean(rpall$PP_no[-1])+sd(rpall$PP_no[-1]),0,mean(rpall$NN_no[-1])+sd(rpall$NN_no[-1]),0,
              mean(rpall$PN_no[-1])+sd(rpall$PN_no[-1]),0),
       y1 = c(mean(rpall$PP_no[-1])-sd(rpall$PP_no[-1]),0,mean(rpall$NN_no[-1])-sd(rpall$NN_no[-1]),0,
              mean(rpall$PN_no[-1])-sd(rpall$PN_no[-1]),0),
       angle = 90,code = 3,length = 0.1)
dev.off()

# PP no
# plot distribution of randoms and compare to real
svg(file='PSN_PPno_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PP_no[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PP_no)-20,max(rpall$PP_no)+20),xaxt='n',yaxt='n')
abline(v=rpall$PP_no[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP',ylab ="Frequency Density",
      xlab="Number of Pairs",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PN no
# plot distribution of randoms and compare to real
svg(file='PSN_PNno_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PN_no[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PN_no)-20,max(rpall$PN_no)+20),xaxt='n',yaxt='n')
abline(v=rpall$PN_no[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN',ylab ="Frequency Density",
      xlab="Number of Pairs",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# NN no
# plot distribution of randoms and compare to real
svg(file='PSN_NNno_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$NN_no[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$NN_no)-20,max(rpall$NN_no)+20),xaxt='n',yaxt='n')
abline(v=rpall$NN_no[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN',ylab ="Frequency Density",
      xlab="Number of Pairs",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# random sampling based p values
length(which(rpall$PP_no<=rpall$PP_no[1]))/length(rpall$PP_no)
length(which(rpall$PN_no<=rpall$PN_no[1]))/length(rpall$PN_no)
length(which(rpall$NN_no<=rpall$NN_no[1]))/length(rpall$NN_no)
length(which(rpall$PO_no<=rpall$PO_no[1]))/length(rpall$PO_no)
length(which(rpall$NE_no<=rpall$NE_no[1]))/length(rpall$NE_no)

length(which(rpall$PP_totdeg<=rpall$PP_totdeg[1]))/length(rpall$PP_totdeg)
length(which(rpall$PP_indeg<=rpall$PP_indeg[1]))/length(rpall$PP_indeg)
length(which(rpall$PP_outdeg<=rpall$PP_outdeg[1]))/length(rpall$PP_outdeg)
length(which(rpall$PP_posdeg<=rpall$PP_posdeg[1]))/length(rpall$PP_posdeg)
length(which(rpall$PP_negdeg<=rpall$PP_negdeg[1]))/length(rpall$PP_negdeg)
length(which(rpall$PP_inposdeg<=rpall$PP_inposdeg[1]))/length(rpall$PP_inposdeg)
length(which(rpall$PP_innegdeg<=rpall$PP_innegdeg[1]))/length(rpall$PP_innegdeg)
length(which(rpall$PP_outposdeg<=rpall$PP_outposdeg[1]))/length(rpall$PP_outposdeg)
length(which(rpall$PP_outnegdeg<=rpall$PP_outnegdeg[1]))/length(rpall$PP_outnegdeg)

length(which(rpall$PN_totdeg<=rpall$PN_totdeg[1]))/length(rpall$PN_totdeg)
length(which(rpall$PN_indeg<=rpall$PN_indeg[1]))/length(rpall$PN_indeg)
length(which(rpall$PN_outdeg<=rpall$PN_outdeg[1]))/length(rpall$PN_outdeg)
length(which(rpall$PN_posdeg<=rpall$PN_posdeg[1]))/length(rpall$PN_posdeg)
length(which(rpall$PN_negdeg<=rpall$PN_negdeg[1]))/length(rpall$PN_negdeg)
length(which(rpall$PN_inposdeg<=rpall$PN_inposdeg[1]))/length(rpall$PN_inposdeg)
length(which(rpall$PN_innegdeg<=rpall$PN_innegdeg[1]))/length(rpall$PN_innegdeg)
length(which(rpall$PN_outposdeg<=rpall$PN_outposdeg[1]))/length(rpall$PN_outposdeg)
length(which(rpall$PN_outnegdeg<=rpall$PN_outnegdeg[1]))/length(rpall$PN_outnegdeg)

length(which(rpall$NN_totdeg<=rpall$NN_totdeg[1]))/length(rpall$NN_totdeg)
length(which(rpall$NN_indeg<=rpall$NN_indeg[1]))/length(rpall$NN_indeg)
length(which(rpall$NN_outdeg<=rpall$NN_outdeg[1]))/length(rpall$NN_outdeg)
length(which(rpall$NN_posdeg<=rpall$NN_posdeg[1]))/length(rpall$NN_posdeg)
length(which(rpall$NN_negdeg<=rpall$NN_negdeg[1]))/length(rpall$NN_negdeg)
length(which(rpall$NN_inposdeg<=rpall$NN_inposdeg[1]))/length(rpall$NN_inposdeg)
length(which(rpall$NN_innegdeg<=rpall$NN_innegdeg[1]))/length(rpall$NN_innegdeg)
length(which(rpall$NN_outposdeg<=rpall$NN_outposdeg[1]))/length(rpall$NN_outposdeg)
length(which(rpall$NN_outnegdeg<=rpall$NN_outnegdeg[1]))/length(rpall$NN_outnegdeg)


# plot PP,PN,NN ints (cooccurence) numbers random vs real
svg(file='PSN_randoms_co.svg',width=4.5, height=4)
plotb = barplot(c(mean(rpall$PP_PP_sum[-1]),rpall$PP_PP_sum[1],
                  mean(rpall$PN_PN_sum[-1]),rpall$PN_PN_sum[1],
                  mean(rpall$NN_NN_sum[-1]),rpall$NN_NN_sum[1],
                  mean(rpall$PP_PN_sum[-1]),rpall$PP_PN_sum[1],
                  mean(rpall$PP_NN_sum[-1]),rpall$PP_NN_sum[1],
                  mean(rpall$PN_NN_sum[-1]),rpall$PN_NN_sum[1]),
                names.arg=c('PP-PP R','PP-PP','PN-PN R','PN-PN',
                            'NN-NN R','NN-NN','PP-PN R','PP-PN',
                            'PP-NN R','PP-NN','PN-NN R','PN-NN'),
                ylim=c(0,20000),main='PSN Real vs Random Mutual Pair ',
                ylab='Number of Co-occurences',col=c('gray','darkred'))
arrows(x0 = plotb,                      
       y0 = c(mean(rpall$PP_PP_sum[-1])+sd(rpall$PP_PP_sum[-1]),0,
              mean(rpall$PN_PN_sum[-1])+sd(rpall$PN_PN_sum[-1]),0,
              mean(rpall$NN_NN_sum[-1])+sd(rpall$NN_NN_sum[-1]),0,
              mean(rpall$PP_PN_sum[-1])+sd(rpall$PP_PN_sum[-1]),0,
              mean(rpall$PP_NN_sum[-1])+sd(rpall$PP_NN_sum[-1]),0,
              mean(rpall$PN_NN_sum[-1])+sd(rpall$PN_NN_sum[-1]),0),
       y1 = c(mean(rpall$PP_PP_sum[-1])-sd(rpall$PP_PP_sum[-1]),0,
              mean(rpall$PN_PN_sum[-1])-sd(rpall$PN_PN_sum[-1]),0,
              mean(rpall$NN_NN_sum[-1])-sd(rpall$NN_NN_sum[-1]),0,
              mean(rpall$PP_PN_sum[-1])-sd(rpall$PP_PN_sum[-1]),0,
              mean(rpall$PP_NN_sum[-1])-sd(rpall$PP_NN_sum[-1]),0,
              mean(rpall$PN_NN_sum[-1])-sd(rpall$PN_NN_sum[-1]),0),
       angle = 90,code = 3,length = 0.1)
dev.off()

# PP-PP no
# plot distribution of randoms and compare to real
svg(file='PSN_PP-PPco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PP_PP_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PP_PP_sum)-20,max(rpall$PP_PP_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$PP_PP_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP-PP',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PN-PN no
# plot distribution of randoms and compare to real
svg(file='PSN_PN-PNco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PN_PN_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PN_PN_sum)-20,max(rpall$PN_PN_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$PN_PN_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN-PN',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# NN-NN no
# plot distribution of randoms and compare to real
svg(file='PSN_NN-NNco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$NN_NN_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$NN_NN_sum)-20,max(rpall$NN_NN_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$NN_NN_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN NN-NN',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PP-PN no
# plot distribution of randoms and compare to real
svg(file='PSN_PP-PNco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PP_PN_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PP_PN_sum)-20,max(rpall$PP_PN_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$PP_PN_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP-PN',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PN-NN no
# plot distribution of randoms and compare to real
svg(file='PSN_PN-NNco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PN_NN_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PN_NN_sum)-20,max(rpall$PN_NN_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$PN_NN_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PN-NN',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# PP-NN no
# plot distribution of randoms and compare to real
svg(file='PSN_PP-NNco_randvsreal.svg',width=4.5, height=4)
plot(density(rpall$PP_NN_sum[-1],na.rm=TRUE),lwd = 1.5,main='',xlab='',ylab='',
     xlim = c(min(rpall$PP_NN_sum)-20,max(rpall$PP_NN_sum)+20),xaxt='n',yaxt='n')
abline(v=rpall$PP_NN_sum[1],lwd=2.5,col='darkred')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='PSN PP-NN',ylab ="Frequency Density",
      xlab="Number of Co-occurences",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()


# random sampling based p values
length(which(rpall$PP_PN_sum<=rpall$PP_PN_sum[1]))/length(rpall$PP_PN_sum)
length(which(rpall$PP_NN_sum<=rpall$PP_NN_sum[1]))/length(rpall$PP_NN_sum)
length(which(rpall$PN_NN_sum<=rpall$PN_NN_sum[1]))/length(rpall$PN_NN_sum)
length(which(rpall$PP_PP_sum<=rpall$PP_PP_sum[1]))/length(rpall$PP_PP_sum)
length(which(rpall$PN_PN_sum<=rpall$PN_PN_sum[1]))/length(rpall$PN_PN_sum)
length(which(rpall$NN_NN_sum<=rpall$NN_NN_sum[1]))/length(rpall$NN_NN_sum)

stopImplicitCluster()