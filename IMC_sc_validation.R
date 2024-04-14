library(Seurat)
library(CMScaller)

Idents(GSE164690) = "patient"

test_mat = AverageExpression(GSE164690,assay = "RNA",
                             features = rownames(GSE164690),
                             slot = "scale.data")
test_res = ntp(test_mat,template,distance = "pearson")

