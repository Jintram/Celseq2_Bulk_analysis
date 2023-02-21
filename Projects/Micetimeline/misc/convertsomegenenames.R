convertMouseGeneList

genes_of_interest = c('MYL2', 'TNNI3', 'TNNC1', 'MYL3', 'FOXN3', 'ZEB1')

convertMouseGeneList(x = genes_of_interest, the_other_animal = 'mouse')

toString(sort(genes_of_interest))
mouse_genes = convertMouseGeneList_usingJAX(x = genes_of_interest)
toString(mouse_genes)


convertMouseGeneList_usingJAX(c('MYL2', 'TNNI3', 'TNNC1', 'MYL3', 'FOXN3', 'ZEB1','TPM1'))
convertMouseGeneList_usingJAX(c('RPL32','HPRT1'))



