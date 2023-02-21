
library(openxlsx)
library(ggplot2)
library(pals)
alph2=pals::alphabet2(); names(alph2)=NULL

myHeartData =
    #openxlsx::read.xlsx('/Users/m.wehrens/Data/2022_09_micetimeline/tables.mice.parameters/data_thomas_202211.xlsx', startRow = T)
    openxlsx::read.xlsx('/Users/m.wehrens/Data/2022_09_micetimeline/tables.mice.parameters/data_thomas_202211_v2-all.xlsx', startRow = T)


myHeartData$HWdivBW = myHeartData$heart_weight/myHeartData$body_weight
    
ggplot(data = myHeartData, mapping= aes(x=age_days, y=HWdivBW, fill=genotype))+
    geom_bar(stat = 'identity', position='dodge')+theme_bw()+
    xlab('Mouse age (wks)')+ylab('Heart weight / TL')

heartBWmean = 
    sapply(split( myHeartData$HWdivBW, paste0(myHeartData$age_days, myHeartData$genotype)), mean, na.rm=T)

ListnormBy =
    paste0(myHeartData$age_days, 'WT')

myHeartData$lbl = paste0(myHeartData$age_days, myHeartData$genotype)
myHeartData$genotype_ = factor(myHeartData$genotype, levels = c('WT','HOM'))
myHeartData$age_days_ = factor(myHeartData$age_days, levels=sort(unique(myHeartData$age_days)))
myHeartData$HWdivBW_Norm=myHeartData$HWdivBW/heartBWmean[ListnormBy]

myHeartData_avg =
    aggregate(x= list(HWdivBW_Norm=myHeartData$HWdivBW_Norm), by=list(age_days=myHeartData$age_days, genotype=myHeartData$genotype), FUN=mean, na.rm=T)
myHeartData_avg$genotype_ = factor(myHeartData_avg$genotype, levels = c('WT','HOM'))
myHeartData_avg$age_days_ = factor(myHeartData_avg$age_days, levels=sort(unique(myHeartData$age_days)))

ggplot(data = myHeartData, mapping= aes(x=age_days_, y=HWdivBW_Norm, fill=genotype_))+
    #geom_point(position='dodge', alpha=.5, mapping=aes(color=genotype_))+
    #geom_boxplot(position=position_dodge(width=.5), width=.25, outlier.colour = 'grey', outlier.size = 2)+
    geom_bar(data=myHeartData_avg, position=position_dodge(width=.5), mapping=aes(group=genotype_), shape=4, stat='identity', width=.5)+
    geom_point(position=position_dodge(width=.5), mapping=aes(group=genotype_))+#, shape=4)+
    #geom_jitter(position=position_dodge(width=.5), mapping=aes(group=genotype_), shape=4)+
    theme_bw()+
    xlab('Mouse age (days)')+ylab('Heart weight / body weight (%WT)')+
    scale_fill_manual(values=c('darkgrey','firebrick'))+theme(legend.position = c(.1,.9))+#'bottom')
    ylim(c(0,1.2*max(myHeartData$HWdivBW_Norm, na.rm=T)))#+give_better_textsize_plot(10)

# Export data for Easy use in Prism


#position=position_dodge(width=.5), width=.25, outlier.colour = 'grey')


export_df_1 = 
    data.frame(genotype=myHeartData$genotype,
               HWBW = myHeartData$HWdivBW_Norm,
               age = myHeartData$age_days)

repnr = max(table(export_df_1[,c('genotype','age')]))

myages=sort(unique(export_df_1$age))
df_WT_perage_ = data.frame( 
    lapply(myages,
           function(age) {
               output=export_df_1[export_df_1$age==age&export_df_1$genotype=='WT',]$HWBW
               c(output,rep(NA, repnr-length(output)))
           }))
colnames(df_WT_perage_) = myages

df_HOM_perage_ = data.frame( 
    lapply(myages,
           function(age) {
               output=export_df_1[export_df_1$age==age&export_df_1$genotype=='HOM',]$HWBW
               c(output,rep(NA, repnr-length(output)))
           }))
colnames(df_HOM_perage_) = myages

export_df_2 = cbind(t(df_WT_perage_), t(df_HOM_perage_))
export_df_final  = data.frame(export_df_2)

openxlsx::write.xlsx(x = export_df_final, file = '/Users/m.wehrens/Data/2022_09_micetimeline/tables.mice.parameters/data_thomas_202211_v2-simplifiedForPrism_byR.xlsx')

