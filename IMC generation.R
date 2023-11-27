##Generation of IMC classification
##Copy Number Variation (CNV), Somatic Mutation (Mut), Methylation (Meth), mRNA expression (mRNA), miRNA expression (miRNA)
##data were downloaded from TCGA-HNSC project.

library(SNFtool)
library(pheatmap)
comsam = Reduce(intersect, list(colnames(cnv),colnames(mut),colnames(meth),colnames(mrna),colnames(mirna)))

list = list(mut = as.matrix(mut[,comsam]),
            cnv = as.matrix(cnv[,comsam]),
            meth = as.matrix(meht[,comsam]),
            mrna = as.matrix(mrna[,comsam]),
            mirna = as.matrix(mirna[,comsam]))

para.K = 50
para.alpha = 0.6
para.T = 50

W.list <- list()
for (m in 1:length(moic.list)) {
  W.list[[m]] <- affinityMatrix(as.matrix(dist(t(moic.list[[m]]))), K = para.K, sigma = para.alpha)
}

W <- SNF(W.list, K = para.K, t = para.T) 

n.clust <- 2:15
for (n in n.clust) {
  message("--generating heatmap with cluster number of ",n,"...")
  clust <- spectralClustering(W, K = n)
  ind <- sort(as.vector(clust), index.return = TRUE)
  ind <- ind$ix
  diag(W) <- median(as.vector(W))
  indata <- normalize(W)
  indata <- indata + t(indata)
  hm <- pheatmap(indata[ind, ind],
                 color = NMF:::ccRamp(c("white","red"),64),
                 border_color = NA,
                 cluster_cols = F,
                 cluster_rows = F,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 name = "SNF",
                 filename = paste0("snf heatmap with cluster number of ",n,".pdf"))
}

##Line chart
library(SNFtool)
estimateNumberOfClustersGivenGraph(W,2:10)
degs = rowSums(W)

# compute unnormalized Laplacian
degs[degs == 0] = .Machine$double.eps    
D = diag(degs)    
L = D - W
Di = diag(1 / sqrt(degs))
L = Di %*% L %*% Di

# compute the eigenvectors corresponding to the k smallest
# eigs$valuess
eigs = eigen(L)
eigs_order = sort(eigs$values, index.return=T)$ix
eigs$values = eigs$values[eigs_order]
eigs$vectors = eigs$vectors[, eigs_order]
eigengap = abs(diff(eigs$values))
eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])

NUMC <- 2:15
quality = list()
for (c_index in 1:length(NUMC)) {
  ck = NUMC[c_index]
  UU = eigs$vectors[, 1:ck]
  EigenvectorsDiscrete <- .discretisation(UU)[[1]]
  EigenVectors = EigenvectorsDiscrete^2
  
  # MATLAB: sort(EigenVectors,2, 'descend');
  temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
  temp1 <- t(apply(temp1, 1, sort, TRUE))  
  
  quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
    sum(sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
}

estimateNumberOfClustersGivenGraph(W,2:15)
t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)
t2 <- sort(unlist(quality), index.return=TRUE)

K1 = NUMC[t1$ix[1]]
K2 <- NUMC[t2$ix[1]]

plotdata1 <- data.frame(Eigengap = t1$x,
                        Cluster = t1$ix+1)
plotdata2 <- data.frame(Rotation.Cost = t2$x,
                        Cluster = t2$ix+1)
plotdata <- merge(plotdata1,plotdata2,by="Cluster")

plotdata$Rotation.Cost <- plotdata$Rotation.Cost/10000


library(plotrix)
library(latticeExtra)
obj1 <- xyplot(Eigengap~Cluster,plotdata,type="l",lwd=2,
               col = "blue",xlim = 1:16,skip=F,scales = list(tick.number = 15))
obj2 <- xyplot(Rotation.Cost~Cluster,plotdata,
               type="l",lwd=2,col="red",xlim = 1:16,skip=F,scales = list(tick.number=15))
doubleYScale(obj1,obj2,add.ylab2 = TRUE,use.style = FALSE)

n = 3
clust <- spectralClustering(W, K = n)
ind <- sort(as.vector(clust), index.return = TRUE)
ind <- ind$ix
diag(W) <- median(as.vector(W))
indata <- normalize(W)
indata <- indata + t(indata)

sam_order <- row.names(indata[ind, ind])
sam.info = data.frame(
  sample = sam_order,
  cluster = ind
)

##Multi-omics heatmap generation
library(ComplexHeatmap)
library(data.table)
library(circlize)

sample.info <- sample.info[order(sample.info$cluster),]
annCol.tcga <- data.frame(row.names = sample.info$sample,
                          cluster = sample.info$cluster)
annCol.tcga$cluster <- paste("IMC",sample.info$cluster,sep = "")

ha1 <- HeatmapAnnotation(
  "IMC" = annCol.tcga$cluster,
  which = "column",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    "IMC" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 2,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "IMC"
    )), 
  simple_anno_size = unit(0.3,"cm"),
  col = list(
    "IMC" = c(
      "IMC1" = "#00468BFF",
      "IMC2" = "#ED0000FF",
      "IMC3" = "#42B540FF"
    )))

rna.plotdata <- t(scale(t(rna.plotdata)))

ht_mrna <- Heatmap(
  matrix = rna.plotdata,
  name = "mRNA/z-score",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    annCol.tcga$cluster,
    levels = unique(annCol.tcga$cluster)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  col = colorRamp2(seq(-2, 2, by = 4/99), pal_rna),
  show_column_names = FALSE,
  show_row_names = FALSE
)

meth.plotdata <- t(scale(t(meth.plotdata)))
ht_meth <- Heatmap(
  matrix = meth.plotdata,
  name = "Methylation/M-value",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    annCol.tcga$cluster,
    levels = unique(annCol.tcga$cluster)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows = F,
  col = colorRamp2(seq(-8, 8, by = 16/99), pal_protein),
  show_column_names = FALSE,
  show_row_names = FALSE
)


library(gplots)
mirna_pal <- colorpanel(64,low="blue",mid = "white",high="red")
mirna.plotdata <- t(scale(t(mirna)))
mirna.plotdata <- na.omit(mirna.plotdata)

ht_mirna <- Heatmap(
  matrix = mirna.plotdata,
  name = "miRNA/z-score",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    annCol.tcga$cluster,
    levels = unique(annCol.tcga$cluster)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  col = colorRamp2(seq(-2, 2, by = 4/63), mirna_pal),
  show_column_names = FALSE,
  show_row_names = FALSE
)


cnv.plotdata <- cnv[,rownames(annCol.tcga)]
mad <- apply(cnv.plotdata, 1, mad)
cnv.plotdata <- cnv.plotdata[rev(order(mad))[1:5000],]

ht_cnv <- Heatmap(
  matrix = cnv.plotdata,
  name = "CNV",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    annCol.tcga$cluster,
    levels = unique(annCol.tcga$cluster)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  col = colorRamp2(seq(-2, 2, by = 4/63), mirna_pal),
  show_column_names = FALSE,
  show_row_names = FALSE,
  top_annotation = ha1
)

mut.plotdata <- mut[apply(mut, 1, function(x)sum(x==0)<0.98*ncol(mut)),rownames(annCol.tcga)]
ht_mut <- Heatmap(
  matrix = mut.plotdata,
  name = "Mutation",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    annCol.tcga$cluster,
    levels = unique(annCol.tcga$cluster)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  col = c("white","black"),
  show_column_names = FALSE,
  show_row_names = FALSE
)

pdf(
  "multi-omics_ht.pdf",
  width = 17,
  height = 20
)
ht_list <- 
  ht_cnv %v%
  ht_mut %v%
  ht_meth %v%
  ht_mrna %v%
  ht_mirna

draw(
  ht_list,
  column_title = "TCGA", 
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  ht_gap = unit(1, "mm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
