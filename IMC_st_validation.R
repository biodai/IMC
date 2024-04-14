##comb_as is the combined data of GSE208253
##gene_list is the list of IMC marker genes.
gene_list1 = list[[1]]
gene_list2 = list[[2]]
gene_list3 = list[[3]]

comb_as <-AddModuleScore(comb_as,features = gene_list1,name = "IMC1",replace = T)
comb_as <-AddModuleScore(comb_as,features = gene_list2,name = "IMC2")
comb_as <-AddModuleScore(comb_as,features = gene_list3,name = "IMC3")

pdf("IMC1.pdf",width = 5.5,height=6)
print(SCpubr::do_NebulosaPlot(comb_as, features = "IMC1", method = "wkde"))
dev.off()

pdf("IMC2.pdf",width = 5.5,height=6)
print(SCpubr::do_NebulosaPlot(comb_as, features = "IMC2", method = "wkde"))
dev.off()

pdf("IMC3.pdf",width = 5.5,height=6)
print(SCpubr::do_NebulosaPlot(comb_as, features = "IMC3", method = "wkde"))
dev.off()
