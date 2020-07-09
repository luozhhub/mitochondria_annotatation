
library(circlize)
root='.'
setwd(root)
MAF_0=read.csv("Original_all_site_add_age_add_annotation.csv",header = T,stringsAsFactors = F)
dat=read.csv("transv|i_human_MT_mutation_annotation.csv",stringsAsFactors = F)
pos_whole=read.csv("gene_annotation.MT.csv",stringsAsFactors = F)
pos_prot=read.table("woDloop_mito_anno_file.txt",stringsAsFactors = F,header = T)
pos_prot=pos_prot[,-3]
colnames(pos_prot)<-c('Position','Type')
# base df
pos_DL=pos_whole[which(pos_whole$Type=='D-loop'),]
pos=rbind(pos_prot,pos_DL)
# change pos increasing 
pos=pos[with(pos,order(Position)),]
# set cirle base
pos$Type[which(pos$Position>=16024)]="D-loop'"
fa=unique(pos$Type)
fa=factor(fa,levels =fa)

xlim=matrix(0,nrow=39,ncol=2)
div_bin=c(1)
for(i in 1:(nrow(pos)-1)){
  if(pos$Type[i]!=pos$Type[i+1]){
    div_bin=append(div_bin,pos$Position[i+1])
  }
}
div_bin=div_bin[-c(10,11,12,13,seq(28,117),seq(124,135))]
xlim[,1]=div_bin
xlim[,2]=c(div_bin[-1]-1,16569)
row.names(xlim)<-c(unique(pos$Type))

# analyzing set
trmt=MAF_0[,c(4,11,14,17)]
# fit factor names to prot
for(i in 1:nrow(trmt)){
  for(j in 1:nrow(xlim)){
    if(xlim[j,1]<=trmt$POS[i] &(xlim[j,2]>=trmt$POS[i])){
      trmt$Type[i]=rownames(xlim)[j]
      break
    }
  }
}
ctrl=trmt[which(trmt$Diagnosis=="Control"),]
disease=trmt[which(trmt$Diagnosis!="Control"),]

# circular figure
# set the canvas parameter 
circos.par$cell.padding = c(0.02, 0, 0.02, 0)
circos.par$start.degree=90
circos.par$track.margin=c(mm_h(0.001),mm_h(0.001))
circos.par$gap.degree=0.2
circos.par$canvas.ylim=c(-1.2,1.2)
# set 38 color for each gene
# generate by rand_col(number,luminosity='dark')
col_set=c("#BC890BFF" ,"#0A3C64FF" ,"#DF0793FF", "#024A81FF", "#0A1D84FF" ,"#0B1475FF" ,"#430C81FF" ,"#B56E0DFF", "#390787FF",
          "#DA880CFF", "#B97A06FF" ,"#B807EAFF", "#5B7601FF" ,"#737F06FF", "#C89103FF", "#0A9A72FF","#5F0A9CFF","#510A74FF",
          "#49AE10FF" ,"#540488FF", "#BB017CFF", "#C67E07FF", "#063464FF", "#00911BFF" ,"#CE0BB8FF", "#CA0CBCFF","#0A8B09FF",
          "#030D73FF", "#808600FF", "#A50607FF", "#D101BFFF" ,"#09107BFF" ,"#0C8169FF", "#0A9C81FF", "#C7016AFF","#190D8FFF",
          "#047731FF", "#41077DFF","#BC890BFF")
col_fun = colorRamp2(seq(1,39), col_set)


circos.initialize(factors =fa , xlim=xlim)

circos.track(factors =trmt$Type , x= trmt$POS, 
             ylim=c(0,1),track.height =0.05,bg.col=col_set,bg.border=NA,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 4, CELL_META$sector.index,
                           cex=0.3,niceFacing = T,facing = 'clockwise')
             })
circos.track(factors =trmt$Type , x= trmt$POS,  y=trmt$MLE,
             ylim=c(0,1),track.height =0.9,bg.border=NA,
             panel.fun = function(x, y) {
               circos.points(x, y,pch =16,cex = 0.4,col = col_fun(CELL_META$sector.numeric.index))
             })
# circos.yaxis(side='right',labels.niceFacing = TRUE,labels.cex = 0.3) 
circos.clear()
# circular figure:For Control
circos.par$cell.padding = c(0.02, 0, 0.02, 0)
circos.par$start.degree=90
circos.par$track.margin=c(mm_h(0.001),mm_h(0.001))
circos.par$gap.degree=0.2
circos.par$canvas.ylim=c(-1.2,1.2)
circos.initialize(factors =fa , xlim=xlim)

circos.track(factors =ctrl$Type , x= ctrl$POS, 
             ylim=c(0,1),track.height =0.05,bg.col=col_set,bg.border=NA,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 4, CELL_META$sector.index,
                           cex=0.3,niceFacing = T,facing = 'clockwise')
             })
circos.track(factors =ctrl$Type , x= ctrl$POS,  y=ctrl$MLE,
             ylim=c(0,1),track.height =0.9,bg.border=NA,
             panel.fun = function(x, y) {
               circos.points(x, y,pch =16,cex = 0.4,col = col_fun(CELL_META$sector.numeric.index))
             })
# circos.yaxis(side='right',labels.niceFacing = TRUE,labels.cex = 0.3) 
circos.clear()
# circular figure:For Disease
circos.par$cell.padding = c(0.02, 0, 0.02, 0)
circos.par$start.degree=90
circos.par$track.margin=c(mm_h(0.001),mm_h(0.001))
circos.par$gap.degree=0.2
circos.par$canvas.ylim=c(-1.2,1.2)
circos.initialize(factors =fa , xlim=xlim)

circos.track(factors =disease$Type , x= disease$POS, 
             ylim=c(0,1),track.height =0.05,bg.col=col_set,bg.border=NA,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 4, CELL_META$sector.index,
                           cex=0.3,niceFacing = T,facing = 'clockwise')
             })
circos.track(factors =disease$Type , x= disease$POS,  y=disease$MLE,
             ylim=c(0,1),track.height =0.9,bg.border=NA,
             panel.fun = function(x, y) {
               circos.points(x, y,pch =16,cex = 0.4,col = col_fun(CELL_META$sector.numeric.index))
             })
# circos.yaxis(side='right',labels.niceFacing = TRUE,labels.cex = 0.3) 
circos.clear()



