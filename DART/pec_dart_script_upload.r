##################################################################
##################################################################
## Processing and exploration of DArTseq data used in O'Donnell
## et al. (2022) Molecular and morphological analyses support recognition 
## of Prostanthera volucris (Lamiaceae), a new species from the 
## Central Tablelands of New South Wales
##################################################################
##################################################################
## Requires *SNP_2.csv file supplied by DArT and matching sample
## *metadata.csv in a folder specific to the DArT run (within 
## the data directory). 
##################################################################
## Load packages
##################################################################
library(dartR)
library(StAMPP)
library(gplots)
library(ggplot2)
library(plotly)
library(hierfstat)
library(poppr)
library(reshape2)
library(ape)
library(phangorn)
library(magrittr)
library(devtools)
library(RSplitsTree)
library(cluster)
library(LEA)
library(beepr)
library(RColorBrewer)
library(map)
library(dplyr)
library(heatmaply)
library(leaflet)
library(leaflet.minicharts)
library(conStruct)
library(ggmap)
library(ozmaps)
library(sf)
library(ggspatial)
library(ggrepel)
library(cowplot)
library(treeio)
library(viridis)
library(ggtree)
library(stats)
library(phytools)
library(geiger)
library(picante)
library(tidyverse)
library(scico)
library(viridisLite)
##################################################################
## Setup
##################################################################
## Specify main parameters

setwd("")
overwrite=FALSE
indir=getwd()
analysis_name=""

## Calculate other parameters
dartname=list.files(indir,"Report_DPro21-6054_2_moreOrders_SNP_2.csv")
dartfile=file.path(indir,dartname)
metadataname=list.files(indir,"DPro21-6054_coanalysis_metadata3.csv")
metadatafile=file.path(indir,metadataname)
outdir=file.path("output",analysis_name)
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("temp")) dir.create("temp")
if(dir.exists(outdir)){
  if(overwrite){
    print("Overwriting existing directory")
  }else print("Directory already exists")
} else dir.create(outdir)

##################################################################
## Load data
##################################################################
## Read data file. Note that metadata are not read in here, as
##  the order must be the same, but the dartfile is not in a
##  sensible order.

ECgl=gl.read.dart(filename = "Report_DPro21-6054_2_moreOrders_SNP_2.csv", ind.metafile = "DPro21-6054_coanalysis_metadata3.csv",)

## Read metadata file and sort according to sample order in dartfile
EC.metadata=read.csv("DPro21-6054_coanalysis_metadata3.csv", stringsAsFactors = T) #all populations delineated
EC.metadata=EC.metadata[match(ECgl@ind.names,EC.metadata$inds),]

## Add sample metadata to ECgl
ECgl@pop=factor(EC.metadata$pop)
levels(ECgl@pop)
ECgl <- ECgl[order(ECgl$pop),]

##################################################################
## Filter data - first pass
##################################################################
## Explore data quality
gl.report.bases(ECgl)
gl.report.callrate(ECgl)
gl.report.callrate(ECgl,method = "ind")
gl.report.secondaries(ECgl)
gl.report.monomorphs(ECgl)

##filter secondary SNPs
ECgl2 <-gl.filter.secondaries(ECgl)
##filter by read depth
ECgl3 <- gl.filter.rdepth(ECgl2)
##filter by reproducibility (repavg)
ECgl4 <- gl.filter.reproducibility(ECgl3,  threshold=0.98)
##filter by locus callrate
ECgl5 <- gl.filter.callrate(ECgl4, method="loc",threshold=0.95)
##filter by individual callrate
ECgl6 <- gl.filter.callrate(ECgl5, method="ind",threshold=0.8)
ECgl6$ind.names
##check number of inds and loci
nInd(ECgl6)
nLoc(ECgl6)

##################################################################
## Identify clones vs unique individuals
##################################################################

##separate pops
ECgl7 <- gl.filter.callrate(ECgl6, threshold = 1)
clone.pops <- seppop(ECgl6)

###use technical duplicates to determine clone threshold
## first - calculate the proportion of shared alleles between inds from Towac pop
res <- gl.propShared(clone.pops$gilesiiTowac)

## convert shared allele matrix to data frame 
clones_temp <- as.data.frame(as.table(as.matrix(res)))

#identify the threshold from the tech dupes
TDfilter<- filter(clones_temp, Var1%in%c("HZ1","HZ93","JJB3602ee"))
TDfilter <- TDfilter[TDfilter$Var2 %in% c("HZ1.1","HZ93.1","JJB3602ee.1"),]
TDfilter
##threshold is 0.9910912 from JJB3602ee & JJB3602ee.1

clones_temp_2<- clones_temp[clones_temp$Freq>0.9910912,]
clones_temp_2

##remove techdupes from dataset
clones_temp <- subset(clones_temp, !(Var2 %in% c("HZ1.1","HZ93.1","JJB3602ee.1")))
clones_temp_2 <- clones_temp[which(clones_temp$Var1!=clones_temp$Var2),]

# add up sum of clone pairings above threshold (count)
clone_means <- clones_temp_2[-2] %>% group_by(Var1) %>% summarise(mean=mean(Freq),count=sum(Freq>=0.9910912),n=n())
clone_means

##group by number of shared clone pairings and store the number of groups as the likely number of clusters
#this is derived by taking the maximum number of matches from the number sampled + the lowest number

ind_groups <- clone_means %>% group_by(count)

write.csv(ind_groups,file="ind_groups_gl.csv")
ind_groups <- clone_means %>% group_by(count) %>% distinct(count) %>% nrow()

### filter to find individuals with the least amount of clone pairings above the threshold 
##ie most likely to either be unique individuals or individuals with several dupes 
gl_inds <- clone_means[clone_means$count <= 61-(ind_groups),]
gl_inds

#unique inds are "HZ111","JJB3602ee.1","HZ1.1"
#unique walls are "TCW577","RWM7111"


###EVANS CROWN CLONES 
# calculating the proportion of shared alleles

clone.pops

# calculating the proportion of shared alleles
res_EC <- gl.propShared(clone.pops$EvansCrown)

# converting a matrix into columns 
clones_temp <- as.data.frame(as.table(as.matrix(res_EC)))
clones_temp %>% group_by(Var1)
clones_temp
clones_temp_2 <- clones_temp[which(clones_temp$Var1!=clones_temp$Var2),]
# filtering individuals using a threshold of 0.99, add up sum of clone pairings above threshold (count)
clone_means <- clones_temp_2[-2] %>% group_by(Var1) %>% summarise(mean=mean(Freq),count=sum(Freq>=0.9910912),n=n())
clone_means

##group by number of shared clone pairings and store the number of groups as the likely number of clusters

ind_groups <- clone_means %>% group_by(count)

write.csv(ind_groups,file="ind_groups_EC.csv")
ind_groups <- clone_means %>% group_by(count) %>% distinct(count) %>% nrow()

ind_groups

### filter to find individuals with the least amount of clone pairings above the threshold ie most likely to be 
##unique individuals - this is derived by taking the total number of distinct count groups from the total number sampled
gl_inds <- clone_means[clone_means$count <= 36-(ind_groups-10),]
gl_inds


#unique inds are "GMT853c", "RPO30c", "RPO55","RPO56","RPO30f","RPO30g","RPO29a",RPO29i","RPO29b","RPO29j"
##inds selected from remaining groups "GMT853b","RPO30j","RPO29d","RPO28d"

##repeat for remaining pops of phylicifolia - just replace at the res step
##unique phyKNP = "RWM2170961a", "RWN2170961b"
#unique phyDNR = "RWM217099a1","RWM217099b1","RWM217099c1","RWM217099d1","RWM217099e1"
#unique phyAda = "RWM217100a1","RWM217100b1","RWM217100c1","RWM217100d1","RWM217100e1"
#unique phyTinderry = "RPO61","RPO62"
#unique phyCobrunga = "NGW9013","NGW8979"

##make final gl with distinct individuals only

ECgl7 = gl.keep.ind(ECgl6, ind.list = c("HZ111","JJB3602ee.1","HZ1.1","TCW577","RWM7111",
                                        "GMT853c", "RPO30c", "RPO55","RPO56","RPO30f","RPO30g",
                                        "RPO29a","RPO29i","RPO29b","RPO29j","GMT853b","RPO30j",
                                        "RPO29d","RPO28d","RWM2170961a", "RWM2170961b",
                                        "RWM217099a1","RWM217099b1","RWM217099c1","RWM217099d1",
                                        "RWM217099e1","RWM217100a1","RWM217100b1","RWM217100c1",
                                        "RWM217100d1","RWM217100e1","RPO61","RPO62","NGW9013","NGW8979"),
                                        recalc = TRUE, mono.rm = TRUE)

save.image(".Rdata") ##name it what you want

##################################################################
## Ordination - PCA
##################################################################
## Principal coordinates analysis
mypc <- gl.pcoa(ECgl7,nfactors = 6)

## 3D plot
gl.pcoa.plot(mypc, ECgl7,xaxis=1,yaxis=2, zaxis=3, hadjust=1.5, 
             vadjust=1, pop.labels = "pop",axis.label.size=3, pt.size=3)
gl.pcoa.plot(mypc, ECgl7,xaxis=1,yaxis=3, hadjust=1.5, vadjust=1, pop.labels = "pop")

##################################################################
## Neighbour-Network analysis
##################################################################

##create a distance matrix from the dataset
distancegl <- ECgl7
distancegl@pop <- as.factor(ECgl7@ind.names)
dist2splits <- gl.dist.ind(distancegl, method = "euclidian")

##export distance object to nexus to be opened in splitstree
splitstree(dist2splits, nexus.file = "")##name it what you want

##################################################################
## Population-level statistics
##################################################################

##make species genlight

species.gl <- gl.merge.pop(ECgl7, old=c("phyAda","phyCobrunga","phyDNR","phyTinderry","phyKNP"), new="phylicifolia")
species.gl <- gl.merge.pop(species.gl, old=c("gilesiiWalls","gilesiiTowac"), new="gilesii")

##remove pops where N < 5
pop.gl <- gl.drop.pop(ECgl7, pop.list=c("phyCobrunga","phyTinderry","phyKNP"), recalc = TRUE, mono.rm = TRUE)
pop.gl <- gl.merge.pop(pop.gl, old=c("gilesiiWalls","gilesiiTowac"), new="gilesii")

## Compare heterozygosity among species
ho_pop <- gl.report.heterozygosity(species.gl)
pdf(file.path(outdir,paste0(analysis_name,"_ho_species_full.pdf")));gl.report.heterozygosity(species.gl);dev.off()
write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_species_full.txt")))
ho_pop

## Compare heterozygosity among populations
ho_pop <- gl.report.heterozygosity(pop.gl)
pdf(file.path(outdir,paste0(analysis_name,"_ho_pop_full.pdf")));gl.report.heterozygosity(pop.gl);dev.off()
write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_pop_full.txt")))
ho_pop

## Examine private alleles by pairs of populations
gl.report.pa(species.gl)
pa_pw_pop <- gl.report.pa(species.gl)
pa_pw_pop
write.table(pa_pw_pop,file = file.path(outdir,paste0(analysis_name,"_pa_pw_pop.txt")))
global <- gl.basic.stats(ECgl7)
global$overall
write.table(global$overall,file = file.path(outdir,paste0(analysis_name,"_global.txt")))

## Nei's (1972) genetic distance

mygeno <- stamppConvert(ECgl7,type = "genlight")
mydist <- stamppNeisD(mygeno,pop=F)
write.table(mydist,file = file.path(outdir,paste0(analysis_name,"_pop_dist_Nei1972.txt")),sep=",")

## F statistics with STAMPP -  species
pwfst <-stamppFst(species.gl, nboots=1000, percent=95, nclusters=1)
pwfst$Fsts
pwfst$Pvalues
write.table(pwfst$Fsts,file = file.path(outdir,paste0(analysis_name,"_species_fst_full.txt")),sep=",")

## F statistics with STAMPP -  pop
pwfst2 <-stamppFst(pop.gl, nboots=1000, percent=95, nclusters=1)
pwfst2$Fsts
pwfst2$Pvalues
write.table(pwfst2$Fsts,file = file.path(outdir,paste0(analysis_name,"_pop_fst_full.txt")),sep=",")

#####pairwise species fst heatmap

#create species fst objects
species.matFst <- pwfst$Fsts 
species.pvalues <- pwfst$Pvalues

#make melted versions & prep data
melted_species <- melt(species.matFst, na.rm=TRUE)
melted_p_species = melt(species.pvalues, na.rm=TRUE)
melted_species2 <- melted_species %>% mutate_if(is.numeric, round, digits=3)
melted_species3 <- melted_species2[order(melted_species2$value),]
levels(melted_species3$Var1) <- c("Evans Crown","P. gilesii","P. phylicifolia")
levels(melted_species3$Var2) <- c("Evans Crown","P. gilesii","P. phylicifolia")

##ggplot
p = ggplot(data = melted_species3, aes(x=Var1, y=Var2, fill=value)) + 
  labs(x="",y="") + 
  geom_tile(color="white") + 
  scale_fill_distiller(palette = "YlOrRd", na.value = "white", trans = "reverse") + 
  theme_minimal() + 
  theme(axis.text.y=element_text(angle=90,size=12, face="bold.italic")) +
  theme(axis.text.x = element_text(size=12, face="bold.italic")) +
  geom_text(aes(Var1, Var2, label = format(round(value, digits=3), 
                nsmall = 3), color = ifelse(value < 0.5, "white", "black")), size = 8) +
  scale_color_manual(values = c("#FFFFCC", "black"), guide = "none") +
  labs(fill = "Fst") + theme(legend.title = element_text(size = 10, face = "bold.italic"))

p


#####pairwise pop fst heatmap

#create pop fst objects
pop.matFst=pwfst2$Fsts 
pop.pvalues = pwfst2$Pvalues

#make melted versions & prep data
melted_pop <- melt(pop.matFst, na.rm=TRUE)
melted_p_pop = melt(pop.pvalues, na.rm=TRUE)
melted_pop2 <- melted_pop %>% mutate_if(is.numeric, round, digits=3)
melted_pop3 <- melted_pop2[order(melted_pop2$value),]
levels(melted_pop3$Var1) <- c("Evans Crown","P. gilesii","P. phylicifolia\n(Adaminaby)","P. phylicifolia\n(Dangelong NR)")
levels(melted_pop3$Var2) <- c("Evans Crown","P. gilesii","P. phylicifolia\n(Adaminaby)","P. phylicifolia\n(Dangelong NR)")

#make melted versions & prep data
p = ggplot(data = melted_pop3, aes(x=Var1, y=Var2, fill=value)) + 
  labs(x="",y="") + 
  geom_tile(color="white") + 
  scale_fill_distiller(palette = "YlOrRd", na.value = "white", trans = "reverse") + 
  theme_minimal() + 
  theme(axis.text.y=element_text(angle=90,size=12, face="bold.italic")) +
  theme(axis.text.x = element_text(size=12, face="bold.italic")) +
  geom_text(aes(Var1, Var2, label = format(round(value, digits=3), nsmall = 3), 
                color = ifelse(value < 0.2, "white", "black")), size = 8) +
  scale_color_manual(values = c("#FFFFCC", "black"), guide = "none") +
  labs(fill = "Fst") + theme(legend.title = element_text(size = 10, face = "bold.italic"))
p


##################################################################
## SNMF Analysis
##################################################################


## Extract coordinates#####

##generate list of retained samples

keep = ECgl7@ind.names

## Read metadata file and sort according to sample order in dartfile
EC.metadata2=read.csv("DPro21-6054_coanalysis_metadata5.csv", stringsAsFactors = T) #all populations delineated
EC.metadata2=EC.metadata2[match(ECgl7@ind.names,EC.metadata2$id),]
EC.metadata2$fullLabel=paste(EC.metadata2$sp, EC.metadata2$id, sep="_")
EC.metadata3=EC.metadata2[match(ECgl7@ind.names,EC.metadata2$id),]
EC.metadata4=EC.metadata3[EC.metadata3$id %in% keep,]
EC.metadata4=EC.metadata4[match(ECgl7@ind.names,EC.metadata4$id),]
coords <- as.data.frame(EC.metadata4[c("id","lat","long")])
coords2 <- coords
ECgl8 <- ECgl7[order(ECgl7@ind.names),] ##sort gl by id
coords2 = coords[match(ECgl8@ind.names,coords$id),] %>% subset(select=c(lat,long))
coords3 = coords[match(ECgl9@ind.names,coords$id),] ## create a version with IDs for the snmf mapping function later
colnames(coords2)<-c("X", "Y")
ECgl8@other$xy<-coords2
ECgl9 <- ECgl8[order(ECgl8@other$xy$X, ECgl8@other$xy$Y),]

###run snmf analysis####

###create snmf directory on desktop - for some reason it only works if the directory is on the desktop
desktop = ""  ####replace with whatever your desktop path is
if(!dir.exists(paste0(desktop,"SNMF_", analysis_name, sep=''))){
  dir.create(paste0(desktop,"SNMF_", analysis_name, sep=''))
  snmf.dir=paste0(desktop,"SNMF_", analysis_name, sep='')
  setwd(snmf.dir) 
} else {dir.create(paste0(desktop,"SNMF_", analysis_name,"_2", sep=''))
  snmf.dir=paste0(desktop,"SNMF_", analysis_name,"_2", sep='')
  setwd(snmf.dir) 
}

#Step 1. Convert genlight object 'gl' to a structure file using package dartR.

gl2structure(ECgl9, indNames = ECgl9@ind.names, addcolumns = ECgl9@pop, ploidy = 2, 
             exportMarkerNames = TRUE, outfile = paste0(analysis_name,".str", sep=''), outpath = getwd(), verbose = 5)

#Step 2. In Notepad++ first add two new columns names 'ind' and 'pop', 
# separated by tabs to push the locus names by two columns.

#Step 3. Convert structure file to geno object using package LEA.

geno <- struct2geno(input.file= paste0(analysis_name,".str", sep=''), ploidy = 2, FORMAT=2, extra.row=1, extra.column=2) 

##perform snmf runs
# remove.snmfProject("***.str.snmfProject") ## if you need to clear your snmf runs and start from scratch
#project1 = load.snmfProject("***.str.snmfProject") ###use to load existing snmf project from file
project1 = NULL
project1 = snmf(paste0(analysis_name,".str.geno", sep=''), CPU = 3, K = 1:8 ,entropy = TRUE, repetitions = 50,project = "new")


####post analysis data inspection#####

## inspect cross entropy values
##plot cross entropy
dev.off()
plot(project1, col = "blue4", cex = 1.2, pch = 19)
project1

### find lowest average cross entropy values
ce.values = function(k, r, filename){
  ce.values = matrix(nrow=r, ncol=k)
  kvalues <- c(paste("K",1:k, sep=""))
  run.no = c(paste("RUN",1:r, sep="_"))
  colnames(ce.values) <- kvalues
  rownames(ce.values) <- run.no
  for (i in 1:k){
    ce.values[,i] = as.matrix(cross.entropy(project1, K=i))}
  mean = colMeans(ce.values)
  min = which.min(mean)
  ce.values = rbind(ce.values, mean)
  ce.values = list(cross.entropy = ce.values, best.k = min)
  if(missing(filename)){
    print(ce.values)
  } else {
    print(ce.values)
    write.csv(ce.values, filename)
  }
  
}

ce.values(8,50, paste0(analysis_name, '_CE.csv', sep = ''))

##################################################################
## SNMF bar plots
##################################################################

###create snmf plot function
snmf.plot = function(snmf, gl, K){
  dev.off() ## reset graphics engine
  par(mfrow=c(1,(length(K)+1)))
  # my.colors = brewer.pal(8, "Set2")
  my.colors = c("#66C2A5","#8DA0CB","#FC8D62","#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3") ##force the order of colours so species match
  for (i in K){
    best = which.min(cross.entropy(snmf, K = i))
    barchart(snmf, K = i, sort.by.Q = FALSE, lab=TRUE, run = best, border = 1, space = 0,col = my.colors, horiz=TRUE, main = paste("K=",i, sep=''))
  }
  # axis(4, at = 1:length(gl@ind.names),labels = paste(gl@pop, gl@ind.names, sep="_"),las=2,cex.axis = 0.9)
  axis(4, at = 1:length(gl@ind.names),labels = gl@ind.names,las=2,cex.axis = 0.9)
}  ###snmf object, genlight with population details, desired values of K

snmf.plot(project1, ECgl9, 2:8)

##################################################################
## SNMF mapping
##################################################################

##create allscores function

k.allscores.file = function(x){
  k=x
  for (i in 1:k) {
    best = which.min(cross.entropy(project1, K = i))
    path <- paste0(getwd(),"/",analysis_name,".str.snmf/K",i,"/run",best,"/", analysis_name,".str_r",best,".",i,".Q", sep="")
    df = read.csv(path, header = FALSE, sep = ' ')
    mapping_graph = cbind.data.frame(coords3,df)
    headerKallscores <- c("id","Lat","Lon",paste("LEA",i,"_K",1:i, sep=""))
    colnames(mapping_graph) <- headerKallscores
    path2 = paste(getwd(),"/", analysis_name,"_TOTAL_BASIC_PCA_allscores_K",i,".csv",sep="")
    write.table(mapping_graph, path2, sep=",",row.names = F)
  }
}

k.allscores.file(8)  ###change to however many values of K were calculated

#create leaflet map function

piemap.snmf <- function(x){
  pcapath <- paste(analysis_name,'_TOTAL_BASIC_PCA_allscores_K',x,'.csv', sep='')
  mappinggraph <- read.csv(pcapath)
  Speciescolors <- c("#66C2A5","#8DA0CB","#FC8D62","#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  MAP <- leaflet() %>% 
    addTiles() %>% 
    setView( lng = 151.66069, lat = -27.07511, zoom = 4.5) %>%
    addProviderTiles("Esri.WorldImagery")# 
  MAP %>% addMarkers (lng=mappinggraph$Lon, lat=mappinggraph$Lat, label = mappinggraph$id)
  
  PIEMAP <- MAP %>%
    addMinicharts (
      mappinggraph$Lon, mappinggraph$Lat,
      type = "pie",
      showLabels = FALSE,
      labelMinSize = 10,
      popup = popupArgs(noPopup = FALSE, showTitle = TRUE ),
      chartdata = mappinggraph[,4:ncol(mappinggraph)],
      colorPalette = Speciescolors, 
      opacity = 0.9,
      width = 25, 
      height = 25,
      legendPosition = "bottomright") 
  
  PIEMAP
}

piemap.snmf(3) ##change x to whatever value you want to inspect 

#########################################SECONDARY DATASET FOR PHYLOGENY###################
##################################################################
## Setup
##################################################################
## Specify main parameters
  
setwd("") ##set to the dart delivery folder
overwrite=FALSE
indir=getwd()
analysis_name="" ##change to whatever you want

## Calculate other parameters
dartname=list.files(indir,"Report_DPro21-6054_3_moreOrders_SNP_2.csv")
dartfile=file.path(indir,dartname)
ECphy_metadataname=list.files(indir,"DPro21-6054_coanalysis_scutmetadata18.csv")
ECphy_metadatafile=file.path(indir,ECphy_metadataname)
outdir=file.path("output",analysis_name)
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("temp")) dir.create("temp")
if(dir.exists(outdir)){
  if(overwrite){
    print("Overwriting existing directory")
  }else print("Directory already exists")
} else dir.create(outdir)

##################################################################
## Load data
##################################################################
## Read data file. Note that metadata are not read in here, as
##  the order must be the same, but the dartfile is not in a
##  sensible order.

ECphy_gl=gl.read.dart(filename = "Report_DPro21-6054_3_moreOrders_SNP_2.csv", ind.metafile = "DPro21-6054_coanalysis_metadata_051022.csv",)
ECphy_metadata
## Read metadata file and sort according to sample order in dartfile
ECphy_metadata=read.csv("DPro21-6054_coanalysis_metadata_051022.csv", stringsAsFactors = T) #all populations delineated
ECphy_metadata=ECphy_metadata[match(ECphy_gl@ind.names,ECphy_metadata$id),]

## Add sample metadata to ECgl
ECphy_gl@pop=factor(ECphy_metadata$pop)
levels(ECphy_gl@pop)
ECphy_gl <- ECphy_gl[order(ECphy_gl$pop),]
ECphy_gl$ind.names <- as.character(ECphy_gl$pop)
##################################################################
## Filter data - first pass
##################################################################
## Explore data quality
gl.report.bases(ECphy_gl)
gl.report.callrate(ECphy_gl)
gl.report.callrate(ECphy_gl,method = "ind")
gl.report.secondaries(ECphy_gl)
gl.report.monomorphs(ECphy_gl)

##select individuals for further analysis

ECphy_gl2 = gl.keep.ind(ECphy_gl, ind.list = c("EvansCrown_GMT853c","EvansCrown_RPO29a",
                                               "EvansCrown_RPO55","gilesiiTowac_HZ1","gilesiiTowac_JJB3602a",
                                               "gilesiiWalls_RWM7111", "gilesiiWalls_TCW 577",
                                               "phylicifoliaAda_RWM217100a1","phylicifoliaAda_RWM217100b1",            
                                               "phylicifoliaCobrunga_NGW8979","phylicifoliaCobrunga_NGW9013",           
                                               "phylicifoliaDeuaNP_NSW1057814","phylicifoliaDNR_RWM217099a1",            
                                               "phylicifoliaDNR_RWM217099b1","phylicifoliaGelantipy_ROM984","phylicifoliaKNP_RWM2170961a",            
                                               "phylicifoliaKNP_RWM2170962a","phylicifoliaNullica_GPP219","phylicifoliaNullica_MP9561",             
                                               "phylicifoliaTinderry_RPO61","phylicifoliaTinderry_RPO62",  
                                               "scutellarioidesTypeCastlereaghNR_TCW215","graniticaWarrumbungles_KT162",
                                               "graniticaWarrumbungles_NSW1057815","densa_NSW1040587",
                                               "densa_NSW203752","marifolia_NSW841013","marifolia_NSW1038335",
                                               "marifolia_NSW1057839"), recalc = TRUE, mono.rm = TRUE)

##filter secondary SNPs
ECphy_gl2 <-gl.filter.secondaries(ECphy_gl2)
##filter by read depth
ECphy_gl3 <- gl.filter.rdepth(ECphy_gl2)
##filter by reproducibility (repavg)
ECphy_gl4 <- gl.filter.reproducibility(ECphy_gl3,  threshold=0.98)
##filter by locus callrate
ECphy_gl5 <- gl.filter.callrate(ECphy_gl4, method="loc",threshold=0.80)
##filter by individual callrate
ECphy_gl6 <- gl.filter.callrate(ECphy_gl5, method="ind",threshold=0.60)

##check number of inds and loci
nInd(ECphy_gl6)
nLoc(ECphy_gl6)

##export snp matrix to be analysed in svdquartets
gl2svdquartets(ECphy_gl6, outpath=getwd(),".nex",method=2) ##name it what you want
save.image(".Rdata") ##name it what you want

##################################################################
## Plot phylogeny
##################################################################
#setwd() - you can change wd at this point if you want
tree=read.tree("_svd_boot.tre") ##put in your _svd_boot.tre file
plot(tree)
data = read.csv("tree_labels.csv", header=TRUE, 
                row.names=1, stringsAsFactors = TRUE) ##this is a csv with original taxon labels and revised labels for figures

##change labels
orig_labels = tree$tip.label
data=data[match(orig_labels,rownames(data)),]

##add label data to tree object
d <- data.frame(label=orig_labels,label2=data)
tree2 <- full_join(tree,d,by="label")

plotTree(tree,node.numbers=T)
dev.off()

##plot tree - this does the clade highlighting and rounding but I just exported
##the tree and text and made it look nicer in illustrator
P <- ggtree(tree2) + geom_tiplab(aes(label=label2),align=TRUE,linesize=.5,hjust=-.02)  + 
xlim(0,1500) + geom_highlight(node=43, fill="#FC8D62", alpha=.3,type="roundrect") +
geom_highlight(node=40, fill="#66C2A5", alpha=.3,type="roundrect") +
geom_highlight(node=38, fill="#8DA0CB", alpha=.3,type="roundrect") +
geom_highlight(node=29, fill="grey", alpha=.60,type="roundrect") +
geom_text(aes(x=branch,label=round(branch.length),fontface="bold", vjust=-.5),size=3)
P
ggsave(file=".pdf",device = "pdf")


