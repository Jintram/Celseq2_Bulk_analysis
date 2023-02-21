
# Use preprocessing_counttables.R to convert umicount output to tables

data_dir = '/Users/m.wehrens/Data/2022_09_micetimeline/test_libraries/counttables/'

################################################################################

library(pals)
library(patchwork)
cols_alphabet1 = pals::alphabet(26); names(col_alphabet)=NULL
cols_alphabet2 = pals::alphabet2(26); names(col_alphabet2)=NULL

# Load some scripts I wrote previously
cfg=list()
cfg$SCRIPT_LOCATION_MW = '/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/'
cfg$SCRIPT_LOCATION_MW_2 = '/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/'
cfg$species = "mouse"
source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_load_libraries.R')

if (F) {
    
    # NOTE THAT THIS CODE WAS ADDED LATER BUT IS REDUNDANT WITH CODE BELOW
    # DON'T EXECUTE
    # Load some resources
    source('/Users/m.wehrens/Documents/git_repos/Resources_bioinf_RNAseq/human_LR_list_v2.R') 
    source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')
    mouse_ligands = conversion_human_to_mouse_symbols[ligand_list_unique]
    
}

################################################################################
# Notes about samples

# Note that I used barcodes from row M (13), col 21, 23, 24
used_BCs = 12*24+c(21,23,24)

# Annotation:
# BC 309 = sample 1 = WT, 2939M1
# BC 311 = sample 2 = Homo, 2795M2
# BC 312 = sample 3 = MYBPC3 mutant Maya (M82)

################################################################################


raw_data_table_mousetest =
    read.table(paste0(data_dir, 'Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.nonRibo_E99_Aligned.out.counts.tsv'), header = 1)

total_counts_percell =
    aggregate(list(count=raw_data_table_mousetest$count), list(cell=raw_data_table_mousetest$cell), sum)
    
df_toplot_totalsraw=total_counts_percell[total_counts_percell$cell %in% c(309, 311, 312),]
df_toplot_totalsraw$samplenames = factor(c('WT.t','MYBPC3.t','MYBPC3.m'), levels=c('WT.t','MYBPC3.t','MYBPC3.m'))
ggplot(df_toplot_totalsraw)+
    geom_bar(aes(x=samplenames, y=count), stat='identity')+theme_bw()+
    xlab(element_blank())+ylab('Total UMI count')+ggtitle('cDNA library from 10 PCR cycles')+give_better_textsize_plot(12)+
    theme(axis.text.x = element_text(angle = 90))


################################################################################

# Now conversion of these count tables 
data_table_mousetest = list()
data_table_mousetest$cycles10 = 
    read.table(paste0(data_dir, 'Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.nonRibo_E99_Aligned.out.counts.table.tsv'))
data_table_mousetest$cycles12 = 
    read.table(paste0(data_dir, 'Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.nonRibo_E99_Aligned.out.counts.table.tsv'))
data_table_mousetest$cycles10_rawc = 
    read.table(paste0(data_dir, 'Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.counts.table-raw.tsv'))
data_table_mousetest$cycles12_rawc = 
    read.table(paste0(data_dir, 'Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.counts.table-raw.tsv'))
# ERCC data
data_table_mousetest$ERCC_cycles10 = 
    read.table(paste0(data_dir, 'ERCC_Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.nonRibo_E99_Aligned.out.counts.table.tsv'))
data_table_mousetest$ERCC_cycles12 = 
    read.table(paste0(data_dir, 'ERCC_Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.nonRibo_E99_Aligned.out.counts.table.tsv'))





# View(data_table_mousetest)

which_dataset = 'cycles12'
which_dataset = 'ERCC_cycles12'

data_table_mousetest_sel = 
    lapply(data_table_mousetest, function(X) {
        X=X[,c('X309','X311')]
        colnames(X) = c('WT_8wk','Hom_8wk')
        #X=X[,c('X309','X311','X312')]
        #colnames(X) = c('WT_8wk','Hom_8wk','Mybpc3Maya')
        X
        })
names(data_table_mousetest_sel) = names(data_table_mousetest)

# Plot general stats
totals_10cycles = apply(data_table_mousetest_sel$cycles10, 2, sum)
totals_12cycles = apply(data_table_mousetest_sel$cycles12, 2, sum)
totals_10cycles
totals_12cycles
# For ERCC, and also calculate percentage of ERCC reads
totals_10cycles_ERCC = apply(data_table_mousetest_sel$ERCC_cycles10, 2, sum)
totals_12cycles_ERCC = apply(data_table_mousetest_sel$ERCC_cycles12, 2, sum)
totals_10cycles_ERCC
totals_12cycles_ERCC
totals_10cycles_ERCC/totals_10cycles
totals_12cycles_ERCC/totals_12cycles


totals_12cycles_ERCC_aboveX = apply(data_table_mousetest_sel$ERCC_cycles12[apply(data_table_mousetest_sel$ERCC_cycles12[,1:2], 1, sum)>5,][,1:2], 2, sum)
expected_ratio  = totals_12cycles_ERCC[1]/totals_12cycles_ERCC[2]
totals_12cycles_ERCC_aboveX[1]/totals_12cycles_ERCC_aboveX[2]
expected_ratio_log10  = log10(.1+totals_12cycles_ERCC[1])/log10(.1+totals_12cycles_ERCC[2])
expected_ratio_log10_adj  = log10(.1+totals_12cycles_ERCC_aboveX[1])/log10(.1+totals_12cycles_ERCC_aboveX[2])
expected_ratio2 = (totals_12cycles_ERCC/totals_12cycles)[1] / (totals_12cycles_ERCC/totals_12cycles)[2] # very similar
expected_ratio3 = (totals_12cycles_ERCC+totals_12cycles)[1] / (totals_12cycles_ERCC+totals_12cycles)[2] # wouldn't make sense

# ERCC read-out
sum(rownames(data_table_mousetest_sel$ERCC_cycles10) %in% rownames(data_table_mousetest_sel$ERCC_cycles12))/dim(data_table_mousetest_sel$ERCC_cycles12)[1]
data_table_mousetest_sel$ERCC_cycles10$gene = row.names(data_table_mousetest_sel$ERCC_cycles10)
data_table_mousetest_sel$ERCC_cycles12$gene = row.names(data_table_mousetest_sel$ERCC_cycles12)
ggplot(data_table_mousetest_sel$ERCC_cycles12, mapping = aes(x=log10(.1+WT_8wk), y=log10(.1+Hom_8wk)))+
    geom_point()+
    geom_abline(slope = 1, intercept = 0)+
    geom_abline(slope=1/expected_ratio_log10_adj, intercept = 0, lty=2)+
    theme_bw()


data_table_mousetest_sel$ERCC_cycles10$cycles=10
data_table_mousetest_sel$ERCC_cycles12$cycles=12
df_ERCC_merged_melted = rbind(data_table_mousetest_sel$ERCC_cycles10, data_table_mousetest_sel$ERCC_cycles12)
cols_alphabet1_=cols_alphabet1
names(cols_alphabet1_)=NULL
ggplot(df_ERCC_merged_melted, mapping = aes(x=log10(.1+WT_8wk), y=log10(.1+Hom_8wk), color=as.factor(cycles)))+
    geom_point()+
    geom_abline(slope = 1, intercept = 0)+
    geom_abline(slope=1/expected_ratio_log10_adj, intercept = 0, lty=2)+
    theme_bw()+scale_color_manual(values=cols_alphabet1_)+theme(legend.position='none')

# dependency on cycles
df_ERCC_merged = merge(data_table_mousetest_sel$ERCC_cycles10, data_table_mousetest_sel$ERCC_cycles12, by.x = 'gene', by.y = 'gene', suffixes = c('_c10','_c12'))
sel_adj = (df_ERCC_merged$WT_8wk_c10+df_ERCC_merged$WT_8wk_c12)>5
expected_ratio_log10_c  = log10(.1+sum(df_ERCC_merged$WT_8wk_c10))/log10(.1+sum(df_ERCC_merged$WT_8wk_c12))
expected_ratio_log10_c_adj  = log10(.1+sum(df_ERCC_merged$WT_8wk_c10[sel_adj]))/log10(.1+sum(df_ERCC_merged$WT_8wk_c12[sel_adj]))
ggplot(df_ERCC_merged, mapping = aes(x=log10(.1+WT_8wk_c10), y=log10(.1+WT_8wk_c12)))+
    geom_point()+
    geom_abline(slope = 1, intercept = 0)+
    geom_abline(slope=1/expected_ratio_log10_c_adj, intercept = 0, lty=2)+
    theme_bw()+scale_color_manual(values=cols_alphabet1_)+theme(legend.position='none')

# Amount of genes that are equal
sum(rownames(data_table_mousetest_sel$cycles10) %in% rownames(data_table_mousetest_sel$cycles12))
shared_genes = rownames(data_table_mousetest_sel$cycles10)[rownames(data_table_mousetest_sel$cycles10) %in% rownames(data_table_mousetest_sel$cycles12)]

# Similarity between datasets by scattering
df_compare_reads = data.frame(
    WT_10C = data_table_mousetest_sel$cycles10[shared_genes,]$Hom_8wk,
    WT_12C = data_table_mousetest_sel$cycles12[shared_genes,]$Hom_8wk,
    gene_name=shared_genes)
ggplot(df_compare_reads, aes(x=log10(WT_10C+.1),y=log10(.1+WT_12C)))+
    geom_smooth(formula = y~x, data=df_compare_reads[(df_compare_reads$WT_10C+df_compare_reads$WT_12C)>4,])+
    geom_point()+
    theme_bw()+geom_abline(slope = 1, intercept = 0)+give_better_textsize_plot(12)

# Read counts vs. UMI counts (later better version)
shared_genes2=rownames(data_table_mousetest_sel$cycles10)[rownames(data_table_mousetest_sel$cycles10) %in% rownames(data_table_mousetest_sel$cycles10_rawc)]
rownames_list_all=lapply(data_table_mousetest_sel[c('cycles10','cycles12')], rownames)
shared_genes_all=Reduce(intersect, rownames_list_all)
df_compare_reads_rawumi = data.frame(
    UMI = c(data_table_mousetest_sel$cycles10[shared_genes_all,]$Hom_8wk,
            data_table_mousetest_sel$cycles12[shared_genes_all,]$Hom_8wk),
    Cnt = c(data_table_mousetest_sel$cycles10_rawc[shared_genes_all,]$Hom_8wk,
            data_table_mousetest_sel$cycles12_rawc[shared_genes_all,]$Hom_8wk),
    gene_name=shared_genes_all,
    cycles=rep(c('10','12'), each=length(shared_genes_all)))
ggplot(df_compare_reads_rawumi, aes(x=Cnt,y=UMI))+
    #geom_smooth(formula = y~x, data=df_compare_reads_rawumi[(df_compare_reads_rawumi$Cnt+df_compare_reads_rawumi$UMI)>4,])+
    geom_point()+
    theme_bw()
ggplot(df_compare_reads_rawumi, aes(x=log10(Cnt+.1),y=log10(.1+UMI), color=cycles))+
    #geom_smooth(formula = y~x, data=df_compare_reads_rawumi[(df_compare_reads_rawumi$Cnt+df_compare_reads_rawumi$UMI)>4,])+
    geom_point(size=1,alpha=.5)+
    geom_smooth()+
    scale_color_manual(values=cols_alphabet2[c(2,3)])+
    #geom_density_2d(bins=3)+
    theme_bw()

# General stats again
df_totals_Cnt = aggregate(x=list(Cnt=df_compare_reads_rawumi$Cnt), by = list(cycles=df_compare_reads_rawumi$cycles), FUN=sum)
df_totals_UMI = aggregate(x=list(UMI=df_compare_reads_rawumi$UMI), by = list(cycles=df_compare_reads_rawumi$cycles), FUN=sum)
p1=ggplot(df_totals_Cnt, aes(x = cycles,y = Cnt, fill=cycles))+
    geom_bar(stat='identity')+scale_fill_manual(values=cols_alphabet2[c(2,3)])+theme_bw()+give_better_textsize_plot(12)+theme(legend.position = 'none')
p2=ggplot(df_totals_UMI, aes(x = cycles,y = UMI, fill=cycles))+
    geom_bar(stat='identity')+scale_fill_manual(values=cols_alphabet2[c(2,3)])+theme_bw()+give_better_textsize_plot(12)+theme(legend.position = 'none')
p1+p2

# General stats again, no mitochondrial
nomito_selection = !grepl(':mt-',df_compare_reads_rawumi$gene_name)
df_totals_Cnt_nomito = aggregate(x=list(Cnt=df_compare_reads_rawumi$Cnt[nomito_selection]), by = list(cycles=df_compare_reads_rawumi$cycles[nomito_selection]), FUN=sum)
df_totals_UMI_nomito = aggregate(x=list(UMI=df_compare_reads_rawumi$UMI[nomito_selection]), by = list(cycles=df_compare_reads_rawumi$cycles[nomito_selection]), FUN=sum)
p1=ggplot(df_totals_Cnt_nomito, aes(x = cycles,y = Cnt, fill=cycles))+
    geom_bar(stat='identity')+scale_fill_manual(values=cols_alphabet2[c(2,3)])+theme_bw()+give_better_textsize_plot(12)+theme(legend.position = 'none')
p2=ggplot(df_totals_UMI_nomito, aes(x = cycles,y = UMI, fill=cycles))+
    geom_bar(stat='identity')+scale_fill_manual(values=cols_alphabet2[c(2,3)])+theme_bw()+give_better_textsize_plot(12)+theme(legend.position = 'none')
p1+p2+plot_annotation(title='Without mito genes')

# 
x=df_compare_reads_rawumi$Cnt
y=df_compare_reads_rawumi$UMI
bins=10^seq(-1,6,.1)
generate_binned_vals = function(x,y,bins) {
    vals=sapply(1:(length(bins)-1), function(idx) {
        idxs_thisbin = (x>=bins[idx])&(x<=bins[idx+1])
        mean(y[idxs_thisbin],na.rm=T)
        }
        )
    sds=sapply(1:(length(bins)-1), function(idx) {
        idxs_thisbin = x>bins[idx]&x<bins[idx+1]
        sd(y[idxs_thisbin],na.rm=T)
        }
        )
    bincenters= (bins[1:(length(bins)-1)]+bins[ 2:(length(bins))])/2
    return(list(vals=vals,bincenters=bincenters, sds=sds))
}

df_binned_10 = data.frame(generate_binned_vals(df_compare_reads_rawumi[df_compare_reads_rawumi$cycles=='10',]$Cnt,
                                               df_compare_reads_rawumi[df_compare_reads_rawumi$cycles=='10',]$UMI,bins))
df_binned_12 = data.frame(generate_binned_vals(df_compare_reads_rawumi[df_compare_reads_rawumi$cycles=='12',]$Cnt,
                                               df_compare_reads_rawumi[df_compare_reads_rawumi$cycles=='12',]$UMI,bins))

ggplot(df_compare_reads_rawumi, aes(x=log10(Cnt+.1),y=log10(.1+UMI), color=cycles))+
    #geom_smooth(formula = y~x, data=df_compare_reads_rawumi[(df_compare_reads_rawumi$Cnt+df_compare_reads_rawumi$UMI)>4,])+
    geom_point(size=1,alpha=.5)+
    geom_point(data=df_binned_10, aes(x=log10(.1+bincenters),y=log10(.1+vals)), color='black')+
    geom_point(data=df_binned_12, aes(x=log10(.1+bincenters),y=log10(.1+vals)), color='black')+
    geom_line(data=df_binned_10, aes(x=log10(.1+bincenters),y=log10(.1+vals)), color='black')+
    geom_line(data=df_binned_12, aes(x=log10(.1+bincenters),y=log10(.1+vals)), color='black')+
    geom_hline(yintercept = log10(.1+5),linetype='dashed')+
    geom_smooth(formula = y ~ poly(x, 5, raw=TRUE))+
    scale_color_manual(values=cols_alphabet2[c(2,3)])+
    #geom_density_2d(bins=3)+
    theme_bw()+xlim(c(0,NA))+ylim(c(0,NA))+give_better_textsize_plot(12)

#pals::pal.bands(cols_alphabet1)

# ggplot(df_compare_reads_rawumi, aes(x=Cnt,color=cycles))+
#     geom_freqpoly(bins=100)+
#     theme_bw()
ggplot(df_compare_reads_rawumi, aes(x=log10(Cnt+.1),color=cycles))+
    #geom_smooth(formula = y~x, data=df_compare_reads_rawumi[(df_compare_reads_rawumi$Cnt+df_compare_reads_rawumi$UMI)>4,])+
    geom_freqpoly()+
    theme_bw()
ggplot(df_compare_reads_rawumi, aes(x=log10(UMI+.1),color=cycles))+
    #geom_smooth(formula = y~x, data=df_compare_reads_rawumi[(df_compare_reads_rawumi$Cnt+df_compare_reads_rawumi$UMI)>4,])+
    geom_freqpoly()+
    theme_bw()


# # Compare maya's Mybpc3 with Thomas' one
# df_compare_reads2 = data.frame(
#     Wk8 = data_table_mousetest_sel$cycles12[shared_genes,]$WT_8wk,
#     Maya = data_table_mousetest_sel$cycles12[shared_genes,]$Mybpc3Maya,
#     gene_name=shared_genes)
# ggplot(df_compare_reads2, aes(x=log10(Wk8+.1),y=log10(.1+Maya)))+
#     geom_smooth(formula = y~x, data=df_compare_reads2[(df_compare_reads$Wk8+df_compare_reads$Maya)>4,])+
#     geom_point()+
#     theme_bw()+geom_abline(slope = 1, intercept = 0)
mean_sum=mean(sapply(data_table_mousetest_sel$cycles12, function(X){sum(X)}))
data_table_mousetest_sel_norm =  data.frame(apply(data_table_mousetest_sel$cycles12,2,function(X){
                                        X/sum(X)*mean_sum
                                        }))

# Now also a table without mitochondrial genes
data_table_mousetest_sel_c12_nomito = data_table_mousetest_sel$cycles12[!grepl(':mt-',rownames(data_table_mousetest_sel$cycles12)),]
mean_sum2=mean(sapply(data_table_mousetest_sel_c12_nomito, function(X){sum(X)}))
data_table_mousetest_sel_c12_nomito_norm =  data.frame(apply(data_table_mousetest_sel_c12_nomito,2,function(X){
                                            X/sum(X)*mean_sum2
                                        }))

data_table_mousetest_s12sn = data_table_mousetest_sel_c12_nomito_norm[apply(data_table_mousetest_sel_c12_nomito_norm, 1, sum)>10,]

data_table_mousetest_s12sn[grepl(':Mybpc3',row.names(data_table_mousetest_s12sn)),]

DE_genes_test = 
    data.frame(FC=data_table_mousetest_s12sn$Hom_8wk/data_table_mousetest_s12sn$WT_8wk, gene=row.names(data_table_mousetest_s12sn))
gene_names_split = str_split(string = DE_genes_test$gene, pattern = ':')
DE_genes_test$sym = sapply(gene_names_split, function(X){X[[2]]})
DE_genes_test$ens = sapply(gene_names_split, function(X){X[[1]]})
row.names(DE_genes_test) = row.names(data_table_mousetest_s12sn)

dim(DE_genes_test)
View(DE_genes_test)

DE_genes_test[grepl('Tnni3',DE_genes_test$gene),]
DE_genes_test[grepl('Mybpc3',DE_genes_test$gene),]
DE_genes_test[grepl('Csrp3',DE_genes_test$gene),]

df_FC_sarco =
    DE_genes_test[DE_genes_test$sym %in% c('Actc1','Myl3','Myl2','Mybpc3','Csrp3','Myh7','Tpm1','Tnni3','Tnnc1','Tnnt2','Myzo2'), ]
ggplot(df_FC_sarco)+
    geom_bar(aes(x=sym, y=log2(FC)),stat='identity')+theme_bw()+give_better_textsize_plot(12)+
    theme(axis.text.x = element_text(angle = 90))

df_qPCR_comparison =
    DE_genes_test[DE_genes_test$sym %in% c('Hprt','Rpl32','Csrp3','Mybpc3','Nppa','Nppb', 'Myl2','Myl3','Tnnc1','Tnni3','Foxn3','Zeb1'), ]
ggplot(df_qPCR_comparison)+
    geom_bar(aes(x=sym, y=log2(FC)),stat='identity')+theme_bw()+give_better_textsize_plot(12)+
    theme(axis.text.x = element_text(angle = 90))


View(DE_genes_test[,c('sym','FC')])

data_table_mousetest_s12sn_ext = data_table_mousetest_s12sn
gene_names_split2 = str_split(string = row.names(data_table_mousetest_s12sn_ext), pattern = ':')
data_table_mousetest_s12sn_ext$sym = sapply(gene_names_split2, function(X){X[[2]]})
data_table_mousetest_s12sn_ext$ens = sapply(gene_names_split2, function(X){X[[1]]})

data_table_mousetest_s12sn_ext[grepl('Flnc',data_table_mousetest_s12sn_ext$sym),]
data_table_mousetest_s12sn_ext[grepl('Psma5',data_table_mousetest_s12sn_ext$sym),]

################################################################################
# Some more specific genes of interest
genes_cl2mod4 = c('CMYA5', 'XIRP2', 'ZNF106', 'MAP4', 'NRAP', 'ANKRD1', 'TTN', 'HIPK2', 'DES', 'TECRL', 'TCAP', 'MYOM1')
genes_cl2mod4_mouse = conversion_human_to_mouse_symbols[genes_cl2mod4]

View(
    DE_genes_test[DE_genes_test$sym %in% genes_cl2mod4_mouse,]
)

toString(genes_cl2mod4_mouse[!(genes_cl2mod4_mouse %in% DE_genes_test$sym)])

# Plus TFs
genes_cl2mod4_relTFs = c('SCL', 'TEAD4', 'TCF7L2', 'GATA1', 'MEF2C')
genes_cl2mod4_relTFs_mouse = conversion_human_to_mouse_symbols[genes_cl2mod4_relTFs]

View(
    DE_genes_test[DE_genes_test$sym %in% genes_cl2mod4_relTFs,]
)
# None are available


################################################################################
################################################################################
################################################################################

# Trying to extract some information from this n=1

################################################################################
# First decide which of these genes I consider "significant" (n=1, but ..)
DE_genes_test$FC_negposfold = DE_genes_test$FC
DE_genes_test$FC_negposfold[DE_genes_test$FC<1] = -1/(DE_genes_test$FC_negposfold[DE_genes_test$FC<1])

ggplot(DE_genes_test, aes(x=FC_negposfold))+
    geom_histogram(bins=100)+theme_bw()+
    geom_vline(xintercept = c(2,3,4,5))

thelimits = calc_limits(DE_genes_test$FC_negposfold, percentile = .05)
thelimits

# DE_genes_test_sel = DE_genes_test[DE_genes_test$FC>3|DE_genes_test$FC<1/3,]
DE_genes_test_sel = DE_genes_test[DE_genes_test$FC>(2.7)|DE_genes_test$FC<1/(2.4),]
dim(DE_genes_test_sel)


# Let's check for TFs from Trrust
DE_genes_test_sel[DE_genes_test_sel$sym %in% TF_mouse_trrust_table$V1,]
# --> 
# Conclusion: Ankrd1 2.8105976, Sfpq   0.3279031, Dbp    2.9511275

# Now let's check for TFs from Lambert
# Sanity check: only 82% of human TFs have mouse orthologs
#sum(human_TFs_ens %in% names(convert_enshuman_to_ensmouse_lookup))/length(human_TFs_ens) 
# Sanity sanity check: ±all genes can be found in my downloaded db
#sum(human_TFs_ens %in% names(conversion_ens_sym_107manual_lookup))/length(human_TFs_ens)
#human_TFs_ens[!(human_TFs_ens %in% names(conversion_ens_sym_107manual_lookup))]
# Now produce the mouse TFs:
mouse_TFs_convertedFromLambert2018 = convert_enshuman_to_ensmouse_lookup[human_TFs_ens[human_TFs_ens %in% names(convert_enshuman_to_ensmouse_lookup)]]
names(mouse_TFs_convertedFromLambert2018)=NULL
mouse_TFs_convertedFromLambert2018

DE_genes_test_sel[DE_genes_test_sel$ens %in% mouse_TFs_convertedFromLambert2018,]
DE_genes_test_sel[DE_genes_test$sym=='Mafk',] # not detected, conversion_human_to_mouse_symbols['MAFK'] --> gives Mafk
# --> 
# Conclusion: Also Dbp, and Nfe2l1 2.951127



# Both methods combined:
myfavoriteTFs=unique(c(rownames(DE_genes_test_sel[DE_genes_test_sel$ens %in% mouse_TFs_convertedFromLambert2018,]), 
                       rownames(DE_genes_test_sel[DE_genes_test_sel$sym %in% TF_mouse_trrust_table$V1,])))
View(DE_genes_test_sel[myfavoriteTFs,c('sym','ens','FC_negposfold')])                     

# Now check for ligands:
ENSCONVERSIONTABLEDIR='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/'
ENSEMBL_VERSION_FORTHISONE=93
source('/Users/m.wehrens/Documents/git_repos/Resources_bioinf_RNAseq/human_LR_list_v2.R')
source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')
ligand_list_unique 
load(paste0(ENSCONVERSIONTABLEDIR,'Rdata/sym_to_ens_conv_table__',ENSEMBL_VERSION_FORTHISONE,'_v2.Rdata'))  # sym_to_ens_conv_table__XX_v2
# Convert ligand list to ens
# sum(ligand_list_unique %in% names(sym_to_ens_conv_table__XX_v2))/length(ligand_list_unique) # Most can be converted to ens
ligand_list_unique[!(ligand_list_unique %in% names(sym_to_ens_conv_table__XX_v2))]
ligand_list_unique_ens = sym_to_ens_conv_table__XX_v2[ligand_list_unique[ligand_list_unique %in% names(sym_to_ens_conv_table__XX_v2)]]
print(paste0('I wasnt able to convert: ', toString(ligand_list_unique[!(ligand_list_unique %in% names(sym_to_ens_conv_table__XX_v2))])))
# Now convert the ligands to mouse ligands
ligand_list_unique_ens_miceHomologues = convert_enshuman_to_ensmouse_lookup[ligand_list_unique_ens[ligand_list_unique_ens %in% names(convert_enshuman_to_ensmouse_lookup)]]
# Now show DE ligands
DE_genes_test_sel[DE_genes_test_sel$ens %in% ligand_list_unique_ens_miceHomologues,]


View(DE_genes_test_sel[DE_genes_test_sel$ens %in% ligand_list_unique_ens_miceHomologues,c('sym','ens','FC_negposfold')])







DE_genes_test_temp = DE_genes_test
rownames(DE_genes_test_temp) = shorthand_cutname(rownames(DE_genes_test_temp))
DE_genes_test_temp$FC_negposfold=round(DE_genes_test_temp$FC_negposfold,2)
View(DE_genes_test_temp[conversion_human_to_mouse_symbols[c('YY1', 'NFE2L1', 'NFIC', 'SRF', 'MAFK', 'PIN1', 'TSC22D1', 'DBP', 'ANKRD1', 'SFPQ', 'JUND', 'CALM3', 'HSP90AA1', 'PKM', 'RTN4', 'PSAP', 'MFGE8', 'PROS1', 'TGM2', 'CCN2', 'CALR', 'GNAS', 'LPL')],])







