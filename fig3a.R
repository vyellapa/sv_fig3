library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library(dplyr)
options(scipen=999)
setwd("/Users/yellapav/SV_project/data/")
###### Create bins ################
cyto.bed=read.table("cytoband_b37.bed",sep="\t")
colnames(cyto.bed)=c("chr","start","stop","band","direction")
cyto <- data.frame(chrom=character(), 
                   start=numeric(), 
                   stop=numeric(),
                   arm=character(), 
                   stringsAsFactors=FALSE)


#Create cytoband df with min,max for each arm
cyto.bed$arm=paste(cyto.bed$chr,substr(cyto.bed$band,0,1),sep="")
for(i in unique(cyto.bed$arm)) {subset=cyto.bed[cyto.bed$arm==i,]; 
min=min(subset$start)
max=max(subset$stop)
cyto <- rbind(cyto, data.frame(chrom=unique(subset$chr), start=min, stop=max,arm=i))
}

bin.size=10000
#Create dataframe with bins of bin.size into df called bins.df
bins.df <- data.frame(chrom=character(), 
                      start=numeric(), 
                      stop=numeric(),
                      mid=numeric(),
                      bin.num=character(), 
                      stringsAsFactors=FALSE)
bin=0
sa=0

for(i in unique(cyto$arm)){
  sub=cyto[cyto$arm==i,]
  chrom=sub$chrom
  start=sub$start
  stop=sub$stop
  
  
  starts=seq(start, stop, by=bin.size)
  stops=starts+bin.size
  chroms=rep(as.character(chrom),length(starts))
  arms=rep(as.character(i),length(starts))
  bins.df=rbind(bins.df, data.frame(chrom=chroms, start=starts, stop=stops, mid=starts,bin.num=arms))
}

bins.df$bin.num=paste(bins.df$bin.num,bins.df$start,bins.df$stop,sep="_")




### upload latest SV file
sv_all<- read.delim("delly_mapq_60_all.txt", sep="\t", stringsAsFactors = F)
head(sv_all)
sv <- sv_all[!sv_all$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[!sv$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[sv$sample %in% grep("_1_BM",unique(sv$sample),value=TRUE),]

#A unique key for each SV
sv$key=paste(sv$sample,sv$chrom1,sv$pos1,sv$chrom2,sv$pos2,sv$SVTYPE,sep="_")

#Merge 2 SV break points into a melted df
sv1=sv[,c("sample","key","chrom1","pos1")]
sv2=sv[,c("sample","key","chrom2","pos2")]
colnames(sv2)=colnames(sv1)
sv.bed=rbind(sv1,sv2)

#A wonky SV
sv.bed=sv.bed[!sv.bed$key=="MMRF_1079_2_BM__34523640_2_34526079_DEL",]

#Get inter SV break-point
sv.bed$dist=0
sv.bed$pos1=as.numeric(as.character(sv.bed$pos1))
sv.bed=sv.bed[order(sv.bed$chrom1,sv.bed$pos1),]

for(i in 1:dim(sv.bed)[1]) {
  if(i>1 && (sv.bed[i,c("chrom1")] == sv.bed[(i-1),c("chrom1")]) ) {sv.bed[i,c("dist")]=sv.bed[i,c("pos1")]-sv.bed[i-1,c("pos1")]}
}


gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
sv.bed=sv.bed[!is.na(sv.bed$pos1),]
gr.sv1 = GRanges(seqnames=Rle(as.character(sv.bed$chrom1)), IRanges(sv.bed$pos1, sv.bed$pos1+1), sample=sv.bed$sample, key=sv.bed$key)

overlapGenes <- findOverlaps(gr.bins, gr.sv1)
df.sv = data.frame(sv.bed[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])

#write.table(df.sv,file="~/Desktop/SV_paper/sv_cyto_10kb.txt", append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE,quote=F)



#sv.binned = data.frame(chrom=character(), 
#                       pos=numeric(), 
#                       bin.num=character(),
#                       sample=character(), 
#                       stringsAsFactors=FALSE)


######### Prepare CNV ###########
mb.lim=as.numeric(as.character(1000000))
cnv = read.table("commpass_cnv_new_fm6.txt",sep="\t",header=T,quote='~')
cnv.n = cnv %>% dplyr::mutate(num.markers=1000) %>% dplyr::select("IDA","seqnames","startA","endA","major") %>% 
  dplyr::rename(ID=IDA, chr=seqnames,start=startA, end=endA,seg.mean=major) %>% dplyr::filter(seg.mean!="2") #%>% dplyr::filter(chr==8)

########################################################################
####gb=read.table("~/local/resources/GRCh37.e75.gene_boundaries.bed", header=F, sep="\t")
####gr.bins = GRanges(seqnames=Rle(gb$V1), IRanges(gb$V2, gb$V3))
#gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
#gr.sv1 = GRanges(seqnames=Rle(as.character(cnv.n$chr)), IRanges(cnv.n$start, cnv.n$end))
#overlapGenes <- findOverlaps(gr.bins, gr.sv1)
#see=data.frame(cnv.n[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])
#write.table(see,file="~/Desktop/SV_paper/cnv_see.txt", append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE,quote=F)


#del.df = see %>% dplyr::filter(seg.mean<2) %>% unique()
#amp.df = see %>% dplyr::filter(seg.mean>2) %>% unique()

####################################################
sv.bed.m = sv.bed %>% dplyr::mutate(start=pos1-mb.lim,stop=pos1+mb.lim)
cnv.bed.m = cnv.n %>% dplyr::filter(seg.mean!=2)

cnv.bed.new = data.frame(head(cnv.bed.m,n=1),head(sv.bed.m,n=1))

for(i in unique(sv.bed.m$sample)){
  sv.sub = sv.bed.m[sv.bed.m$sample==i,]
  cnv.sub = cnv.bed.m[cnv.bed.m$ID==i,] 
  
  gr.bins = GRanges(seqnames=Rle(sv.sub$chrom1), IRanges(sv.sub$start, sv.sub$stop))
  gr.sv1 = GRanges(seqnames=Rle(as.character(cnv.sub$chr)), IRanges(cnv.sub$start, cnv.sub$end))
  overlapGenes <- findOverlaps(gr.bins, gr.sv1)
  cnv.bed.new = rbind(cnv.bed.new,data.frame(cnv.sub[subjectHits(overlapGenes),], sv.sub[queryHits(overlapGenes),]))
  
  
}

i="MMRF_1940_1_BM"
sv.sub = sv.bed.m[sv.bed.m$sample==i,]
cnv.sub = cnv.bed.m[cnv.bed.m$ID==i,] 

gr.bins = GRanges(seqnames=Rle(sv.sub$chrom1), IRanges(sv.sub$start, sv.sub$stop))
gr.sv1 = GRanges(seqnames=Rle(as.character(cnv.sub$chr)), IRanges(cnv.sub$start, cnv.sub$end))
overlapGenes <- findOverlaps(gr.bins, gr.sv1)
cnv.bed.new.1 = data.frame(cnv.sub[subjectHits(overlapGenes),], sv.sub[queryHits(overlapGenes),])


(cnv.bed.new) %>% dplyr::mutate( xx=abs(start-pos1),yy=abs(end-pos1)) %>% dplyr::mutate( xx=abs(start-pos1),yy=abs(end-pos1)) %>% dplyr::filter(xx < mb.lim | yy < mb.lim) %>% head()

cnv.bed.new = cnv.bed.new[-1,]
#cnv.bed.new.bkup = cnv.bed.new
cnv.bed.new = unique(cnv.bed.new) %>% dplyr::mutate(start=as.numeric(as.character(start)),pos1=as.numeric(as.character(pos1)),end==as.numeric(as.character(end))) %>% dplyr::mutate( xx=abs(start-pos1),yy=abs(end-pos1)) %>% dplyr::filter(xx < mb.lim | yy < mb.lim) %>% dplyr::select(ID,chr,start,end,seg.mean)
cnv.bed.new = unique(cnv.bed.new)
## Overlap with 10kb bins
gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
gr.sv = GRanges(seqnames=Rle(as.character(cnv.bed.new$chr)), IRanges(cnv.bed.new$start, cnv.bed.new$end))

overlapGenes <- findOverlaps(gr.bins, gr.sv)
df = data.frame(cnv.bed.new[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])
df = df %>% dplyr::mutate(seg.mean=as.numeric(as.character(seg.mean)))
write.table(df,file="~/Desktop/SV_paper/cnv_cyto_svSpecific_10kb.txt", append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE,quote=F)

#df= data.table::fread("~/Desktop/SV_paper/cnv_cyto_svSpecific_10kb.txt",header=T,sep="\t")
head(df)
del.df = df %>% dplyr::filter(seg.mean<2) %>% unique()
amp.df = df %>% dplyr::filter(seg.mean>2) %>% unique()

head(del.df)
#del.df = del.df %>% dplyr::select(ID,chr,start,end,seg.mean,bin.num) %>% dplyr::group_by(bin.num) %>% tally() %>% data.frame()  %>% left_join(bins.df, by =c('bin.num' = 'bin.num')) 
#amp.df = amp.df %>% dplyr::select(ID,chr,start,end,seg.mean,bin.num) %>% dplyr::group_by(bin.num) %>% tally() %>% data.frame()  %>% left_join(bins.df, by =c('bin.num' = 'bin.num')) 
#del.df$chrom=paste0("chr",del.df$chrom)
#amp.df$chrom=paste0("chr",amp.df$chrom)

head(del.df)

head(bins.df)
bins.df_test<- bins.df[bins.df$start > 9000000 &bins.df$stop< 15000000 & bins.df$chrom==16,]
del.df_test<- del.df[del.df$start.1 > 9000000 &del.df$stop< 15000000 & del.df$chrom==16,]
head(del.df_test)
kk<- as.data.frame(table(del.df_test$bin.num))
colnames(kk)[1]<-"bin.num"
require(plyr)
int<- join(bins.df_test, kk, by="bin.num")
int$Freq[is.na(int$Freq)]<-0
plot(int$Freq, type="l")

###### read hotspots ###############################################
h.spots = read.table("manual.hotspots.txt",header=TRUE,sep="\t")
h.spots = (h.spots) %>% dplyr::mutate(start=start.bp-10000000, end=end.bp+10000000)
h.spots$chr = as.character(gsub("chr","",h.spots$chr))


library("Gviz")
plot.chr="chr8"
plot.start = 126806779
plot.end = 131113499
#Enhancer
enhancers = read.table("Merged.MM.primary.K27ac.remove2.5k.bed",sep="\t")
en.myc = enhancers  %>% dplyr::mutate(start=as.numeric(as.character(V2)), end=as.numeric(as.character(V3))) ## %>% dplyr::filter(V1=="chr8" & start>126806779 & end<131113499)
gr.myc = GRanges(seqnames=Rle(en.myc$V1), IRanges(en.myc$start, en.myc$end))
head(enhancers)

#################
####### Promoter
#################
promoters = read.table("promoter_track.txt",sep="\t",header=T,quote='~')
promoters = promoters[grep("promoter",promoters$Annotation),] %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(end))) #%>% dplyr::filter(chr=="chr8" & start>126806779 & end<131113499)
head(promoters)
pr.gr = GRanges(seqnames=Rle(promoters$chr), IRanges(promoters$start, promoters$end))


#SVs

df.sv.hs.tally = (df.sv) %>% dplyr::group_by(bin.num) %>% tally() %>% data.frame()  %>% left_join(bins.df, by =c('bin.num' = 'bin.num')) %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(stop))) 
df.sv.hs.tally$chrom=paste0("chr",df.sv.hs.tally$chrom)
gr = GRanges(seqnames=Rle(df.sv.hs.tally$chrom), IRanges(df.sv.hs.tally$start, df.sv.hs.tally$end), freq=df.sv.hs.tally$n)
dTrack <- DataTrack(gr, name="SV Breakpoints / Bin",background.title="darkblue",type=c("s"),col=c("#007B22"))

#################
####### CNV #####
#################

del.df.n = del.df %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(stop))) ##%>% dplyr::filter(chrom=="chr8" & start>126806779 & end<131113499)
head(del.df.n)

amp.df.n = amp.df %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(stop))) ##%>% dplyr::filter(chrom=="chr8" & start>126806779 & end<131113499)
head(amp.df.n)
all.df = dplyr::full_join(del.df.n,amp.df.n, by="bin.num") %>% dplyr::rename(del=n.x,amp=n.y)
head(all.df)
all.df[is.na(all.df$chrom.x),c("start.x")]=all.df[is.na(all.df$chrom.x),c("start.y")]
all.df[is.na(all.df$chrom.x),c("end.x")]=all.df[is.na(all.df$chrom.x),c("end.y")]
all.df[is.na(all.df$chrom.x),c("del")]=0
all.df[is.na(all.df$chrom.x),c("chrom.x")]=all.df[is.na(all.df$chrom.x),c("chrom.y")]

cnv.gr = GRanges(seqnames=Rle(all.df$chrom.x), IRanges(all.df$start.x, all.df$end.x), Del=all.df$del,Amp=all.df$amp)
cnvTrack <- DataTrack(cnv.gr, name="CNV / Bin",background.title="darkblue",groups=c("Del", "Amp"),type=c("s","a"))





#################
######## ENSEMBL
#################

genes = read.table("~/Desktop/SV_paper/data/genes.txt",header=F,sep="\t")
genes = read.table("genes1.txt",header=F,sep="\t")
genes.vec = as.vector(genes$V1)
#biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", filters=list(hgnc_symbol=c("MYC","PVT1","PCAT1","ASAP1"),transcript_source="havana"),stacking="squish",background.title="darkblue",transcriptAnnotation="symbol")
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", filters=list(hgnc_symbol=genes.vec,transcript_source="havana"),stacking="squish",background.title="darkblue",transcriptAnnotation="symbol")


#chr <- as.character(unique(seqnames(plot.chr)))
#gen <- genome(gr.myc)
atrack <- AnnotationTrack(gr.myc, name="Enhancers",background.title="darkblue",fill="#6836A5",stacking="dense")
ptrack <- AnnotationTrack(pr.gr, name="Promoter",background.title="darkblue",fill="#21AABD",stacking="dense")

gtrack <- GenomeAxisTrack(transcriptAnnotation="symbol")

####### GISTIC ########
gistic = read.table("all_lesions.conf_90.txt",header=T,sep="\t")
gistic = gistic[,c(4,7)]

x=data.frame(matrix(unlist(strsplit(gsub("-",":",as.character(gistic$Peak.Limits)),"[:,-,(]")),byrow=T,ncol=5))
gistic=cbind(gistic,x)
gistic = (gistic) %>% dplyr::select(2:5) %>% dplyr::mutate(X1=as.character(X1),X2=as.numeric(as.character(X2)),X3=as.numeric(as.character(X3)),neg.log.q=-log(Residual.q.values.after.removing.segments.shared.with.higher.peaks))

bins.df.chr = bins.df %>% dplyr::mutate(chrom=paste0("chr",chrom))
gr.bins.chr = GRanges(seqnames=Rle(bins.df.chr$chrom), IRanges(bins.df.chr$start, bins.df.chr$stop))
gistic.gr = GRanges(seqnames=Rle(gistic$X1), IRanges(gistic$X2, gistic$X3), Del=gistic$neg.log.q)
overlapGenes <- findOverlaps(gr.bins.chr, gistic.gr)
gistic.new = data.frame(gistic[subjectHits(overlapGenes),], bins.df.chr[queryHits(overlapGenes),])
gistic.new=unique(gistic.new)





gistic.gr = GRanges(seqnames=Rle(gistic.new$chrom), IRanges(gistic.new$start, gistic.new$stop), Del=gistic.new$neg.log.q)
gistic.gr = GRanges(seqnames=Rle(gistic.new$chrom), IRanges(gistic.new$start, gistic.new$stop))
#gisticTrack <- DataTrack(gistic.gr, name="GISTIC Peak Limits",background.title="darkblue",type=c("s"))
gistrack <- AnnotationTrack(gistic.gr, name="GISTIC Boundary",background.title="darkblue",fill="#240046",stacking="dense")
head(gistic)


pdf("plots-sv2.pdf", useDingbats=FALSE)

for(i in 1:nrow(h.spots)) {
 # i=72
  sub=h.spots[i,]
  plot.chr=paste0("chr",sub$chr)
  plot.start=as.numeric(as.character(sub$start))+7000000
  plot.stop=as.numeric(as.character(sub$end))-7000000
  
  
  #head(df.sv.hs.tally)
  kk <- df.sv.hs.tally[df.sv.hs.tally$chrom == plot.chr & df.sv.hs.tally$start > plot.start & df.sv.hs.tally$end < plot.stop,]
  #sum(kk$n)
  
  ## select SV catalogue withion range hotspot
  sv.bed.k = sv.bed
  sv.bed.k$chrom1 = paste0("chr",sv.bed.k$chrom1)
  h.spots_sam <-sv.bed.k[sv.bed.k$chrom1==plot.chr & sv.bed.k$pos1 > plot.start & sv.bed.k$pos1 < plot.stop,] 
  
  ## select sampleID with SV within the hotspot
  sample_list_bin <- unique(h.spots_sam$sample)
  
  ## select CNV data only for samples with SV within hotspot
  del.df_bin<- del.df[del.df$ID %in% sample_list_bin,]
  amp_df_bin<- amp.df[amp.df$ID %in% sample_list_bin,]
  
  amp_df_bin2<- amp_df_bin[amp_df_bin$chr==sub$chr & amp_df_bin$start > plot.start & amp_df_bin$end< plot.stop,] 
  del.df_bin2<- del.df_bin[del.df_bin$chr==sub$chr & del.df_bin$start > plot.start & del.df_bin$end< plot.stop,] 
  
  
  kk_del<- as.data.frame(table(del.df_bin2$bin.num))
  
  #### reorder data frame according to the position
 # plot(kk_gdel$Freq)
  
  var1 = unique(as.character(kk_gains$Var1))
  kk_gains<- as.data.frame(table(amp_df_bin2$bin.num))
  
  
  cn.kkk = data.frame(unique(c(as.character(kk_del$Var1),as.character(kk_gains$Var1))))
  colnames(cn.kkk) = c("bin.num")
  cn.kkk = cn.kkk %>% left_join(bins.df, by="bin.num") %>% dplyr::mutate(chrom = paste0("chr",chrom))
  
if(nrow(cn.kkk)>0) {  
  if(nrow(kk_gains)>0) {
    cn.kkk = cn.kkk %>% left_join(kk_gains, by=c("bin.num" = "Var1")) %>% dplyr::rename(amp=Freq) 
  } else {
    cn.kkk$amp=0
  }
  
  if(nrow(kk_del)>0) {
    cn.kkk = cn.kkk %>% left_join(kk_del, by=c("bin.num" = "Var1")) %>% dplyr::rename(del=Freq)
  } else {
    cn.kkk$del=0
  }
} else {
  cn.kkk=data.frame(bin.num=1,chrom="chr1",start=1,stop=1,mid=1,amp=1,del=1 )
} 
  cn.kkk$amp[is.na(cn.kkk$amp)]<-0
  cn.kkk$del[is.na(cn.kkk$del)]<-0
  kkk.gr = GRanges(seqnames=Rle(cn.kkk$chrom), IRanges(cn.kkk$start, cn.kkk$stop), Del=cn.kkk$del,Amp=cn.kkk$amp)
  kkkTrack <- DataTrack(kkk.gr, name="CNV / Bin",background.title="darkblue",groups=c("Del", "Amp"),type=c("s"),col=c("brown2","dodgerblue"))
  
  #### reorder data frame according to the position
  #plot(kk_gains$Freq, type="l")
  
  #plot(cn.kkk$amp, type="l")
  
  print(c(plot.chr,plot.start,plot.stop))
  itrack <- IdeogramTrack(genome="hg19", chromosome=plot.chr)
  #plotTracks(list(itrack, gtrack,dTrack,cnvTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL,
  #           chromosome=plot.chr, from=plot.start, to=plot.stop )
  
  
  plotTracks(list(itrack, gtrack,dTrack,kkkTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL,
             chromosome=plot.chr, from=plot.start, to=plot.stop )
 }

dev.off()
itrack <- IdeogramTrack(genome="hg19", chromosome="chr15")
(plotTracks(list(itrack, gtrack,dTrack,cnvTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL, chromosome="chr15", from=22116803, to=28664453 ))

plotTracks(list(gtrack,biomTrack),chromosome="chr1", from=232702064, to=245702064)

dTrack <- DataTrack(gr, name="SV Breakpoints / Bin",background.title="darkblue",type=c("s"),col=c("#007B22"))
plotTracks(dTrack,chromosome="chr15", from=22116803, to=28664453)

ggplot(ddf,aes(seg.len))+geom_histogram()+facet_wrap(.~as.factor(chr),scales = "free")

par(mfrow = c(4, 1))
c1 = plotTracks(list(itrack, gtrack,dTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL, chromosome="chr1", from=0, to=249250621, panel.only = TRUE)

c2 = plotTracks(list(itrack, gtrack,dTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL, chromosome="chr2", from=0, to=243199373)
c3 = plotTracks(list(itrack, gtrack,dTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL, chromosome="chr3", from=0, to=198022430)

c4 = plotTracks(list(itrack, gtrack,dTrack,gistrack,atrack,ptrack,biomTrack),collapseTranscripts="longest",col = NULL, chromosome="chr4", from=0, to=191154276)
