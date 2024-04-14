##data as bulk sample transcriptomic data
##marker_gene as IMC marker genes
indata <- t(scale(t(tpm_mRNA[marker_gene,])))

library(lolR)
# nearest centroid classifier
classifier <- lol.classify.nearestCentroid(t(indata), as.numeric(factor(group))) 

testdata <- t(scale(t(data[marker_gene,])))

group2 <- paste0("Cluster",predict(classifier, t(testdata)))
names(group2) <- colnames(testdata) 
group2 <- sort(group2) 

group2 = data.frame(
  sample = names(group2),
  cluster = group2
)
