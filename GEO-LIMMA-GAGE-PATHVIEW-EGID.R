options(scipen=999)
#insert the study you want, preferrably small subset
gset <- "GSE84511"
gset_inside <- paste0(gset,"_series_matrix.txt.gz")
#download it from GEO
library(GEOquery)
data.geo <- getGEO(gset)
#pull out all the useful parts but first unwrap data.geo
data.inside <- data.geo$GSE84511_series_matrix.txt.gz #hard coded for now...

#geoquery data obtained this way is already normalized either with RMA or MAS5 or something so we don't fret here


#the actual expression table
expn.raw <- as.matrix(assayData(data.inside)$exprs)
#drop out the row labels
rownames(expn.raw) <- c()
#get list of gene symbols to prepare for collapsing multiple probe to single gene refs
gene.syms <- featureData(data.inside)$`Entrez Gene`


####experimental, revert to just using normal gene.syms if this breaks code#####
##### some probes map to multiple entrez IDs that don't aggregate because they carry different names so we will homogenize
gene.syms<-sapply(gene.syms,function(x) x=str_match(x,"^\\d*")[1,1]) #entries with multiple entrez ID sep like 11111 /// 22222 /// 33333 reduced to 1st ID
names(gene.syms) <- c() #remove the entry names
#####

#get list of probe ids since we dropped them out of expn.raw
#####probe.ids <- featureData(data.inside)$`ID`
#now make the collapsed table and prepare to collapse probes in it
expn.raw <- as.data.frame(expn.raw)
#####collapse <- cbind(gene.syms,probe.ids,expn.raw)
collapse <- cbind(gene.syms,expn.raw)
colnames(collapse)[1] <- "gene"
#####colnames(collapse)[2] <- "probe"




#use aggregate to collapse probes to genes using max function
#####collapse <- aggregate(. ~ gene,data=collapse,FUN=max,na.rm=TRUE)
collapse <- aggregate(. ~ gene,data=collapse,function(x) {x[which.max(abs(x))]})
write.table(collapse,file="collapsed_expression_data.txt",sep = "\t",row.names = FALSE)



#time for limma on the collapsed set...
library(limma)
control.idx <- rep("young",24)
test.idx <- rep("old",24)
groupings <- c(test.idx,control.idx)
f <- factor(groupings,levels=c("old","young"))
design <- model.matrix(~ 0 + f)
colnames(design) <- c("old","young")
#better make a collapsed matrix that's like the example data.matrix...
collapsed.matrix <- collapse
rownames(collapsed.matrix) <- collapse$gene
collapsed.matrix$gene <- c()
collapsed.matrix$probe <- c()
#continuing with limma now
data.fit <- lmFit(collapsed.matrix,design)
contrast.matrix <- makeContrasts(old-young,levels=design)
data.fit.con <- contrasts.fit(data.fit,contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)
#export a ranked list with fold change for GSEA
gsea.out <- cbind(collapse$gene,as.data.frame(data.fit.eb$`coefficients`)) ################ugh
colnames(gsea.out)[1:2] <- c('egid','fc')
write.table(gsea.out,"fc.rnk",sep ="\t")
#do top table to control FDR within DEGs (top 5000, adjust number parameter as needed and BH method is FDR control)
tab <- topTable(data.fit.eb,coef=1,number=21435,adjust.method = "BH")
#select the top genes < 0.05 adjusted p
topgenes <- tab[tab[, "adj.P.Val"] < 0.05,]
write.table(topgenes,"limma_top_genes_p005.txt",sep="\t")
write.table(tab,"limma_top_21435.txt",sep="\t")
write.table(collapsed.matrix,"full_matrix.txt",sep="\t")
#volcano plot
pdf("volcano.pdf")
volcanoplot(data.fit.eb,coef=1,highlight=10)
dev.off() #otherwise graphs don't show


#gage analysis
library(gage) #gage works on entrez gene IDs
mouse.kegg.gset <- kegg.gsets("mmu") #get mouse kegg gene sets (entrez id)
kegg.gs=mouse.kegg.gset$kg.sets[mouse.kegg.gset$sigmet.idx]
gage.yng.idx <- c(25:48) #last 24 samples are young
gage.old.idx <- c(1:24) #first 24 samples are old
collapsed.matrix.firstrowless <- collapsed.matrix[-1,] #remove first row of na probe collapses
collapsed.matrix.firstrowless <- as.matrix(collapsed.matrix.firstrowless) #might need to be a matrix not df for this to work
kegg.out <- gage(collapsed.matrix.firstrowless,kegg.gs,ref = gage.yng.idx,samp = gage.old.idx,same.dir=F,use.fold=T) #kegg by fold change (2 way perturbation)
#kegg.out <- gage(collapsed.matrix.firstrowless,kegg.gs,ref = gage.yng.idx,samp = gage.old.idx,same.dir=F,use.fold=F) #kegg by p val
print("2d perturbed KEGG");head(kegg.out$greater[,1:5],5)
write.table(kegg.out$greater,"KEGG_FC_2D.txt",sep="\t")
kegg.2d.sig <- sigGeneSet(kegg.out,outname="KEGG_sig")
write.table(kegg.2d.sig,"KEGG_sig.txt",sep="\t")

#change 1:10 below for inclusion of more top pathways
#essential unique genes in top 10 pathways
kegg.2d.ess <- esset.grp(kegg.out$greater,collapsed.matrix.firstrowless,gsets=kegg.gs,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="KEGG_ESS_2D1",make.plot=T)
gs=unique(unlist(kegg.gs[rownames(kegg.out$greater)[1:50]]))
kegg.essData=essGene(gs, collapsed.matrix.firstrowless, ref =gage.yng.idx, samp =gage.old.idx)
for (gs in rownames(kegg.out$greater)[1:50]) {
outname = gsub(" |:|/", "_", substr(gs, 10, 100))
geneData(genes=kegg.gs[[gs]], exprs=kegg.essData, ref=gage.yng.idx, samp=gage.old.idx, outname="KEGG_ESS_2D2", txt=T, heatmap=T, Colv=F, Rowv=F, dendrogram="none",limit=3,scatterplot=T)
}



library(GO.db)
go.hs <- go.gsets(species="mouse",id.type="eg") #go analysis
go.bp <- go.hs$go.sets[go.hs$go.subs$BP] #bioproc
go.mf <- go.hs$go.sets[go.hs$go.subs$MF] #molec factor
go.cc <- go.hs$go.sets[go.hs$go.subs$CC] #cell compart

#bioproc
go.out.bp <- gage(collapsed.matrix.firstrowless,go.bp,ref = gage.yng.idx,samp = gage.old.idx,dir=T)
print("GO bioproc // up w/ age"); head(go.out.bp$greater[,c(1,5)],5); print("down w/ age");head(go.out.bp$less[,c(1,5)],5)
write.table(go.out.bp$greater,"GO_BP_FC_UP.txt",sep="\t")
write.table(go.out.bp$less,"GO_BP_FC_DN.txt",sep="\t")
go.bp.sig <- sigGeneSet(go.out.bp,out="GO_BP_SIG")
write.table(go.bp.sig$greater,"GO_BP_FC_UP_SIG.txt",sep="\t")
write.table(go.bp.sig$less,"GO_BP_FC_DN_SIG.txt",sep="\t")

#mol func
go.out.mf <- gage(collapsed.matrix.firstrowless,go.mf,ref = gage.yng.idx,samp = gage.old.idx,dir=T)
print("GO bioproc //up w/ age");head(go.out.mf$greater[,c(1,5)],5); print("down w/ age");head(go.out.mf$less[,c(1,5)],5)
go.out.cc <- gage(collapsed.matrix.firstrowless,go.cc,ref = gage.yng.idx,samp = gage.old.idx,dir=T)
write.table(go.out.mf$greater,"GO_MF_FC_UP.txt",sep="\t")
write.table(go.out.mf$less,"GO_MF_FC_DN.txt",sep="\t")
go.mf.sig <- sigGeneSet(go.out.mf,out="GO_MF_SIG")
write.table(go.mf.sig$greater,"GO_MF_FC_UP_SIG.txt",sep="\t")
write.table(go.mf.sig$less,"GO_MF_FC_DN_SIG.txt",sep="\t")

#cell comp
print("GO bioproc //up w/ age");head(go.out.cc$greater[,c(1,5)],5); print("down w/ age");head(go.out.cc$less[,c(1,5)],5)
write.table(go.out.cc$greater,"GO_CC_FC_UP.txt",sep="\t")
write.table(go.out.cc$less,"GO_CC_FC_DN.txt",sep="\t")
go.cc.sig <- sigGeneSet(go.out.cc,out="GO_CC_SIG")
write.table(go.cc.sig$greater,"GO_CC_FC_UP_SIG.txt",sep="\t")
write.table(go.cc.sig$less,"GO_CC_FC_DN_SIG.txt",sep="\t")


#try to use same esest function on go
#change 1:10 below for inclusion of more top pathways
#bp up
go.bp.ess.up <- esset.grp(go.out.bp$greater,collapsed.matrix.firstrowless,gsets=go.bp,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_BP_UP_ESS",make.plot=T)
#bp down
go.bp.ess.dn <- esset.grp(go.out.bp$less,collapsed.matrix.firstrowless,gsets=go.bp,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_BP_DN_ESS",make.plot=T)
#mf up
go.mf.ess.up <- esset.grp(go.out.mf$greater,collapsed.matrix.firstrowless,gsets=go.mf,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_MF_UP_ESS",make.plot=T)
#mf down
go.mf.ess.dn <- esset.grp(go.out.mf$less,collapsed.matrix.firstrowless,gsets=go.mf,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_MF_DN_ESS",make.plot=T)
#cc up
go.cc.ess.up <- esset.grp(go.out.cc$greater,collapsed.matrix.firstrowless,gsets=go.cc,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_CC_UP_ESS",make.plot=T)
#cc down
go.cc.ess.dn <- esset.grp(go.out.cc$less,collapsed.matrix.firstrowless,gsets=go.cc,ref=gage.yng.idx,samp=gage.old.idx,same.dir=F,outname="GO_CC_DN_ESS",make.plot=T)







#visualize wiring diagram
library(pathview)
collapsed.matrix.firstrowless.d <- collapsed.matrix.firstrowless[,gage.old.idx]-collapsed.matrix.firstrowless[,gage.yng.idx]
path.ids <- c("mmu04630 Jak−STAT signaling pathway","mmu04152 AMPK signaling pathway","mmu04150 mTOR signaling pathway","mmu04060 Cytokine-cytokine receptor interaction","mmu04550 Signaling pathways regulating pluripotency of stem cells")
path.ids2 <- substr(path.ids,1,8) #extract the pathway designation, eg. mmu04630
#native KEGG view
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data=collapsed.matrix.firstrowless.d[,1:2],pathway.id=pid,species="mmu"))





####fast gsea doesn't work right now#####
#fast gsea
library(stringr)
library(fgsea)
data(examplePathways)
my_pathways <- gmtPathways("mGSKB_Entrez.gmt") #use for mmu specific pathways
my_pathways <- examplePathways #doesn't work w mmu genes
gsea.out <- gsea.out[-c(1),, drop=FALSE] #the ,,drop=FALSE is critical to preserve DF struct
rownames(gsea.out) <- c()



gsea.out.test <- gsea.out
gsea.out.num.egid <- sapply(gsea.out$egid,function(x) as.character(x)) #numeric
gsea.out.num.rank <- gsea.out$fc #numeric
names(gsea.out.num.rank) <- gsea.out.num.egid



testfgsea <- fgsea(my_pathways, gsea.out.num.rank, nperm=10000, maxSize=1000,minSize=10)
num_fgsea_hits <- sum(testfgsea[, padj < 0.05])
testfgsea[order(pval), ][1:num_fgsea_hits,]
#fgsea plots
topPathwaysUp <- testfgsea[NES > 0][head(order(padj), n=20), pathway]
topPathwaysDown <- testfgsea[NES < 0][head(order(padj), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf("gsea3.pdf",width=30,height=30)
plotGseaTable(my_pathways[topPathways], gsea.out.num.rank, testfgsea, gseaParam = 0.5)
dev.off()
library(data.table)
fwrite(testfgsea, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

