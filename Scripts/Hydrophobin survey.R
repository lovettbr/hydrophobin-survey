####Load packages####
library(seqinr)
library(effectR)
library(biomartr)
library(bio3d)
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(taxize)
library(gtable)
library(ggforestplot)
library(xlsx)
library(grid)
library(ggpubr)
library(grid)
library(gridExtra)
library(viridis)
library(broom)
library(cowplot)
library(patchwork)
library(agricolae)
library(UpSetR)
library(scales)
library(ggforce)

##In Terminal##
#Load files into folder Proteomes
#Names fasta files with 3 letter abbreviations
#Run rename_headers.sh in Proteomes

####Pull motifs####
#Motif finding function
pull.hyd=function(protein_file, out){
  proteins=seqinr::read.fasta(file=protein_file, forceDNAtolower = F)
  hyd_pattern_strict="CC[^C]+C[^C]+C[^C]+CC"
  hydrophobins=regex.search(proteins, motif = "custom", reg.pat = hyd_pattern_strict)
  seqinr::write.fasta(sequences=hydrophobins, names=names(hydrophobins), nbchar=80, file.out=out)
  print(paste(protein_file, length(hydrophobins), sep=": "))
}

dir="Proteomes/renamed/"
names=list.files(dir)
names=names[grep(".faa", names)]
names=names[grep("candidates", names, invert=T)]

#Pull protein counts
prot.counts=data.frame()
temp=data.frame()

for(i in names){
  genome=i
  proteins=seqinr::read.fasta(file=paste(dir, genome, sep=""), forceDNAtolower = F)
  temp=data.frame(Genome=genome, Count=length(proteins))
  prot.counts=bind_rows(prot.counts, temp)
}

#Pull 6C candidates
for(i in names){
  genome=i
  #pull.hyd(paste(dir, genome, sep=""), paste(dir, "/Hyd_candidates_", genome, sep=""))
}

##In Terminal##
#Run find_hydrophobins.sh (including candidates) arguments "renamed" "renamed_out"
#Move 6C hydrophobin candidates into renamed_out
#Run tablemaker.py in renamed_out > C.stat.tab (Based off biopython cookbook, https://gist.github.com/wdecoster/5782affc0753c9b89a05306fe942a021 and https://stackoverflow.com/questions/69696077/i-want-to-parse-sequences-and-sequence-ids-from-a-fasta-file-and-assign-them-to)
#Combined .faa files into single fall All_candidates.faa in renamed_out
#Run sequence_cleaner_name.py All_candidates.faa producing unique_All_candidates.faa (Based off biopython cookbook "Sequence_Cleaner")

####Stats on motifs motifs with annotations####
C.stat <- read.csv("Proteomes/renamed_out/C.Stat.tab", header=F, sep=" ")
colnames(C.stat)=c("Source", "Accession", "Cysteines", "Length", "Per.C", "Sequence")

C.stat2 = C.stat %>%
  separate(Source, as.character(1:4), sep="_") %>%
  mutate(Accession=str_replace_all(Accession, "#", "_")) %>%
  mutate(Accession=str_replace_all(Accession, "\\|", "_")) %>%
  select(-`1`, -`3`) %>%
  pivot_longer(cols=c(`2`, `4`), names_to="Annotation") %>%
  select(-Annotation) %>%
  filter(!is.na(value)) %>%
  mutate(presence=T) %>%
  pivot_wider(names_from=value, values_from=presence, values_fn = {unique}) %>%
  separate(Accession, into='Genome', sep="_", remove=F) %>%
  relocate(Genome) %>%
  rename(Candidates.6C=candidates, Pfam.Hyd1=hyd1.faa , Pfam.Hyd2=hyd2.faa)

C.stat2[is.na(C.stat2)]=F

####Incorporating exclusion criteria####
exclusion.dat <- read.csv("Input/Final Exclusion.csv")

C.stat.exclusion = left_join(C.stat2, exclusion.dat)

####Genome key and cladistic order####
key.dat <- read.csv("Input/Genome Candidates Hydrophobins Key.csv")

key.dat2 = key.dat %>%
  select(Abbreviation, Kingdom, Phylum, Species, Lifestyle, Host, Niche) %>%
  rename("Genome"=Abbreviation)

class.ids=get_ids(key.dat2$Species, db="gbif", row=1)
first.pass=classification(class.ids$gbif, db="gbif")

key.dat2$Species[is.na(names(first.pass))]="Apophysomyces"

class.ids2=get_ids(key.dat2$Species, db="gbif", row=1)
second.pass=classification(class.ids2$gbif, db="gbif")

clado=class2tree(second.pass)

is_tip <- clado$phylo$edge[,2] <= length(clado$phylo$tip.label)
ordered_tips <- clado$phylo$edge[is_tip, 2]

phylo=key.dat2$Genome[ordered_tips]

phylo=phylo[c(1:8, 10:12, 9, 13:length(phylo))]

##Online
#Submit candidates to SignalP 5.0 server (Eukarya, Long output)
#https://services.healthtech.dtu.dk/service.php?SignalP-5.0

####Incorporating secretion predictions and genome key####
sig.dat <- read.delim("Proteomes/sigp/output_protein_type.txt", header=FALSE, comment="#")

sig.dat2 = sig.dat %>%
  rename(Accession=V1, Prediction=V2, Likelihood=V3, Cleavage=V5) %>%
  mutate(Accession=str_replace(Accession, "\\.1_.*", "\\.1")) %>%
  select(-V4) %>%
  mutate(Prediction=recode_factor(as.factor(Prediction), `SP(Sec/SPI)`="Secreted", `OTHER`="Not"))

C.stat.combined = left_join(C.stat.exclusion, sig.dat2) %>%
  left_join(key.dat2)

C.stat.secreted = C.stat.combined %>%
  filter(Prediction=="Secreted" & Exclusion=="No") %>%
  separate(Cleavage, into=c("First", "Second", "Third", "Start")) %>%
  select(-First, -Second, -Third) %>%
  mutate(Trimmed=str_sub(Sequence, Start, nchar(Sequence))) %>%
  rename(`Cleavage`=`Start`)

C.stat.unsecreted = C.stat.combined %>%
  filter(Prediction=="Not" & Exclusion=="No") %>%
  mutate(Trimmed=Sequence)

C.stat.combined2=bind_rows(C.stat.secreted, C.stat.unsecreted)

C.stat.excluded= C.stat.combined2 %>%
  filter(Exclusion=="No")

#Write out trimmed sequences
#for(i in C.stat.secreted$Genome){
#seqinr::write.fasta(sequences=as.list(C.stat.secreted[C.stat.secreted$Genome==i,]$Trimmed), names=paste("Secreted", C.stat.secreted[C.stat.secreted$Genome==i,]$Accession, sep="_"), file.out=paste("Proteomes/sigp/", i, "_secreted_trimmed.fa", sep=""))
#}

#seqinr::write.fasta(sequences=as.list(C.stat.secreted$Trimmed), names=paste("Secreted", C.stat.secreted$Accession, sep="_"), file.out="Proteomes/sigp/All_secreted_trimmed.fa")

#Write out unsecreted sequences
#for(i in C.stat.unsecreted$Genome){
#  seqinr::write.fasta(sequences=as.list(C.stat.unsecreted[C.stat.unsecreted$Genome==i,]$Sequence), names=C.stat.unsecreted[C.stat.unsecreted$Genome==i,]$Accession, file.out=paste("Proteomes/sigp/", i, "_unsecreted_trimmed.fa", sep=""))
#}

#seqinr::write.fasta(sequences=as.list(C.stat.unsecreted$Sequence), names=paste("Unsecreted", C.stat.unsecreted$Accession, sep="_"), file.out="Proteomes/sigp/All_unsecreted_trimmed.fa")

#Sequences were combined with 29 indicators and aligned via MAFFT in Geneious

####Clustering hydrophobins####
dist <- read.csv("Proteomes/grouping/Excluded distances with all indicators.csv")

dist2=dist %>%
  pivot_longer(cols= 2:length(colnames(dist)), names_to="Match2", values_to="Identity") %>%
  filter(!is.na(Identity))

match_thresh=function(mat, thresh=50){
  colnames(mat)[1]="Match"
  mat2 = mat %>%
    pivot_longer(cols= 2:length(colnames(mat)), names_to="Match2", values_to="Identity") %>%
    filter(!is.na(Identity)) %>%
    filter(Identity>=thresh) %>%
    rowwise() %>%
    mutate(pair=toString(sort(c(Match, Match2)))) %>%
    distinct(pair, .keep_all=T) %>%
    select(Match, Match2) %>%
    mutate(Match=str_replace_all(Match, pattern="[.]", replacement="-"), Match2=str_replace_all(Match2, pattern="[.]", replacement="-"))
}

unionize=function(matches){
  
  first=list(empty="empty")
  for(i in 1:length(matches$Match)){
    append="yes"
    for(j in 1:length(first)){
      first[[j]]=unique(first[[j]])
      
      if(matches[i,]$Match %in% first[[j]]){
        first[[j]]=c(first[[j]], as.character(matches[i,]))
        append="no"
      }
      if(matches[i,]$Match2 %in% first[[j]]){
        first[[j]]=c(first[[j]], as.character(matches[i,]))
        append="no"
      }
    }
    if(append=="yes"){
      first=append(first, list(as.character(matches[i,])))
    }
  }
  
  first=first[2:length(first)]
  names(first)=paste("Group", 1:length(first))
  
  return(first)
}

collapse=function(first){
  second=list(empty="empty")
  for(x in 1:length(first)){
    united="no"
    for(y in 1:length(second)){
      if(length(intersect(second[[y]], first[[x]]))>0){
        second[[y]]=sort(union(second[[y]], first[[x]]))
        second[[y]]=unique(second[[y]])
        united="yes"
      }
    }
    if(united=="no"){
      second[[length(second)+1]]=sort(first[[x]])
      second[[length(second)]]=unique(second[[length(second)]])
    }
  }
  
  second=second[2:length(second)]
  names(second)=paste("Group", 1:length(second))
  return(second)
}

dist2=match_thresh(dist, thresh=25)

first.25=unionize(dist2)

out=collapse(first.25)

group.key=read.csv("Proteomes/grouping/group_key.csv")

out2=data.frame(Group=names(out), Values=unlist(as.character(out))) %>%
  separate(Values, into=as.character(1:1000), sep="[^A-z0-9_.-]") %>%
  select(Group, c(as.character(2:1000))) %>%
  pivot_longer(cols=-Group, names_to="col", values_to="Accession") %>%
  drop_na(Accession) %>%
  filter(Accession!=" " & Accession!="") %>%
  select(-col) %>%
  left_join(group.key) %>%
  select(-Group) %>%
  rename(Group="Old") %>%
  separate(Accession, into=c("Prediction", "Genome"), sep="_", remove=F) %>%
  separate(Accession, into=c("toss", "Accession"), sep="ted_") %>%
  select(-toss) %>%
  mutate(Prediction=str_replace(Prediction, "Unsecreted", "Not")) %>%
  distinct()

out3=out2 %>%
  distinct() %>%
  group_by(Group, Prediction) %>%
  summarize(n=length(Accession), Genomes=toString(unique(Genome)), Genome_count=length(unique(Genome)), Accessions=toString(Accession))

ind=read.csv("Input/Indicator_groups.csv")

indicators=ind %>%
  group_by(Group) %>%
  summarize(Indicators=toString(Indicator), Indicators_n=length(unique(Indicator)))

#No group: ANI_DewA, AFU_RodD, AFU_RodG, AFU_RodB, AFU_RodF, SLA_HYDX, WIC_HYDX     
#Self match: AFU_RodE (Secreted/Unsecreted; Group 39 - dropped), NCR_EAS

theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=30), axis.text.x = element_text(size=12), axis.text.y = element_text(size=15), title = element_text(size=35), legend.title = element_text(size=25), legend.text = element_text(size=20))

plt_sec_group=ggplot(subset(out3, n>2), aes(x=reorder(Group, -n), fill=Prediction, y=n, labels=str_wrap(Accessions)))+geom_col()+theme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Proteins")+xlab("")+scale_fill_viridis_d(option="G")+labs(fill="Secretion")

plt_sec_group

C.stat.excluded2=C.stat.excluded %>%
  mutate(Accession=str_replace_all(Accession, pattern="[.]", replacement="-")) %>%
  mutate(Pfam.Hyd1=recode_factor(as.factor(Pfam.Hyd1), `TRUE`='Hyd1', `FALSE`="NA"),
         Pfam.Hyd2=recode_factor(as.factor(Pfam.Hyd2), `TRUE`='Hyd2', `FALSE`="NA"),
         Candidates.6C=recode_factor(as.factor(Candidates.6C), `TRUE`='6C', `FALSE`="NA")) %>%
  mutate(Groups=str_c(Candidates.6C, Pfam.Hyd1, Pfam.Hyd2, sep="/"))

####Pfam inclusion####
#Proteins in batches of secreted or not-secreted were annotated using PfamA.hmm via hmmscan

secreted_pfam <- read.table("Proteomes/Pfam/secreted_pfam.out", quote="\"")
colnames(secreted_pfam)[c(1,3:5)]=c("Accession", "Pfam_name", "Pfam", "Eval")

secreted_pfam2 = secreted_pfam %>%
  select("Accession", "Pfam_name", "Pfam", "Eval") %>%
  mutate(Accession=str_remove(Accession, "Secreted_")) %>%
  group_by(Accession, Pfam_name) %>%
  tally() %>%
  pivot_wider(names_from=Pfam_name, values_from=n)

unsecreted_pfam <- read.table("Proteomes/Pfam/unsecreted_pfam.out", quote="\"")
colnames(unsecreted_pfam)[c(1,3:5)]=c("Accession", "Pfam_name", "Pfam", "Eval")

unsecreted_pfam2 = unsecreted_pfam %>%
  select("Accession", "Pfam_name", "Pfam", "Eval") %>%
  mutate(Accession=str_remove(Accession, "Unsecreted_")) %>%
  group_by(Accession, Pfam_name) %>%
  tally() %>%
  pivot_wider(names_from=Pfam_name, values_from=n)

all_pfam=rbind(secreted_pfam2, unsecreted_pfam2)
all_pfam[is.na(all_pfam)]=0

length_n=function(col){
  if(length(unique(col))==1 & is.na(unique(col))){0}else(length(unique(col[!is.na(col)])))
}

enlist=function(col){
  if(length(unique(col))==1 & is.na(unique(col))){NA}else(toString(unique(col[!is.na(col)])))
}

all_pfam2=all_pfam %>%
  mutate(Accession=str_replace_all(Accession, pattern="[.]", replacement="-")) %>%
  pivot_longer(cols=colnames(all_pfam)[2:length(all_pfam)], names_to="Pfam", values_to="count") %>%
  filter(count!=0) 

pfam_group=left_join(C.stat.excluded2, out2) %>%
  left_join(all_pfam2) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  left_join(indicators) %>%
  group_by(Group) %>%
  summarize(Pfams=enlist(Pfam), Pfams_n=length_n(Pfam))

pfam_acc=left_join(C.stat.excluded2, out2) %>%
  left_join(all_pfam2) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  left_join(indicators) %>%
  group_by(Accession) %>%
  summarize(Pfams=enlist(Pfam), Pfams_n=length_n(Pfam))

ex_group=left_join(C.stat.excluded2, out2) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  left_join(indicators)

ex_group2=ex_group %>%
  select(Genome, Length, Indicators, Indicators_n, Cysteines, Accession, Phylum, Lifestyle, Prediction, Host, Niche, Groups, Group, Sequence) %>%
  distinct() %>%
  group_by(Group) %>%
  summarize(Phyla_n=length_n(Phylum), Phyla=enlist(Phylum),
            Lifestyles_n=length_n(Lifestyle), Lifestyles=enlist(Lifestyle),
            Hosts_n=length_n(Host), Hosts=enlist(Host),
            Niches_n=length(unique(Niche)), Niches=toString(unique(Niche)),
            Groups_n=length_n(Groups), Groups=enlist(Groups),
            Genomes_n=length_n(Genome), Genomes=enlist(Genome),
            Accessions_n=length_n(Accession), Accessions=enlist(Accession),
            Predictions_n=length_n(Prediction), Predictions=enlist(Prediction),
            Indicators_n=unique(Indicators_n), Indicators=unique(Indicators),
            Indicators_n=ifelse(is.na(Indicators_n), 0, Indicators_n),
            avg_length=round(mean(Length), digits=0),
            avg_C=round(mean(Cysteines), digits=0)) %>%
  left_join(pfam_group) %>%
  filter(!is.na(Group))

#write.csv(ex_group2, "Tables/Grouped_Stats.csv", row.names = F)

ex_ungrouped=ex_group %>%
  select(Genome, Length, Indicators, Indicators_n, Cysteines, Accession, Phylum, Lifestyle, Prediction, Host, Niche, Groups, Group, Sequence) %>%
  group_by(Group, Accession) %>%
  summarize(Phyla=enlist(Phylum),
            Lifestyles=enlist(Lifestyle),
            Hosts=enlist(Host),
            Niches=toString(unique(Niche)),
            Groups=enlist(Groups),
            Genomes=enlist(Genome),
            Predictions=enlist(Prediction),
            avg_length=round(mean(Length), digits=0),
            avg_C=round(mean(Cysteines), digits=0)) %>%
  filter(is.na(Group)) %>%
  ungroup() %>%
  left_join(pfam_acc) %>%
  select(-Group)

#write.csv(ex_ungrouped, "Tables/Ungrouped_Stats.csv", row.names = F)

ex_group3 = ex_group %>%
  mutate(Prediction=ifelse(Prediction=="Not", "Unsecreted", "Secreted")) %>%
  mutate(ID=paste(Prediction, "_", Accession, sep="")) %>%
  left_join(pfam_acc) %>%
  distinct()

Hydrophobin_groups = ex_group %>%
  select(Genome, Length, Indicators, Indicators_n, Cysteines, Accession, Phylum, Lifestyle, Prediction, Host, Niche, Groups, Group, Sequence, Pfam.Hyd1, Pfam.Hyd2) %>%
  filter(!is.na(Indicators) | Pfam.Hyd1=="Hyd1" | Pfam.Hyd2=="Hyd2") %>%
  group_by(Group) %>%
  summarize(Phyla_n=length_n(Phylum), Phyla=enlist(Phylum),
            Lifestyles_n=length_n(Lifestyle), Lifestyles=enlist(Lifestyle),
            Hosts_n=length_n(Host), Hosts=enlist(Host),
            Niches_n=length(unique(Niche)), Niches=toString(unique(Niche)),
            Groups_n=length_n(Groups), Groups=enlist(Groups),
            Genomes_n=length_n(Genome), Genomes=enlist(Genome),
            Accessions_n=length_n(Accession), Accessions=enlist(Accession),
            Predictions_n=length_n(Prediction), Predictions=enlist(Prediction),
            Indicators_n=unique(Indicators_n), Indicators=unique(Indicators)) %>%
  ungroup()

Hyd_groups=unique(Hydrophobin_groups$Group)[!is.na(unique(Hydrophobin_groups$Group))]

Hydrophobins = ex_group %>%
  select(Genome, Length, Indicators, Indicators_n, Cysteines, Accession, Phylum, Lifestyle, Prediction, Host, Niche, Groups, Group, Sequence, Pfam.Hyd1, Pfam.Hyd2) %>%
  group_by(Accession, Phylum, Lifestyle, Host, Niche, Groups, Genome, Prediction, Length, Cysteines) %>%
  summarize(Cluster_n=length_n(Group), Cluster=enlist(Group),
            Indicators_n=unique(Indicators_n), Indicators=unique(Indicators),
            Hyd=length(grep("Hyd", Pfam.Hyd1, fixed=T))+length(grep("Hyd", Pfam.Hyd2, fixed=T))) %>%
  filter(Indicators_n>0 | Hyd>0 | Groups %in% Groups[grep("Hyd", Groups)] | Cluster %in% Hyd_groups) %>%
  ungroup() %>%
  left_join(pfam_acc)

Hydrophobin_accessions=Hydrophobins$Accession
Hydrophobin_groups=Hydrophobins$Cluster

Hydrophobin_groups2 = Hydrophobins %>%
  rename(Group="Cluster") %>%
  group_by(Group) %>%
  summarize(Phyla_n=length_n(Phylum), Phyla=enlist(Phylum),
            Lifestyles_n=length_n(Lifestyle), Lifestyles=enlist(Lifestyle),
            Hosts_n=length_n(Host), Hosts=enlist(Host),
            Niches_n=length(unique(Niche)), Niches=toString(unique(Niche)),
            Groups_n=length_n(Groups), Groups=enlist(Groups),
            Genomes_n=length_n(Genome), Genomes=enlist(Genome),
            Accessions_n=length_n(Accession), Accessions=enlist(Accession),
            Predictions_n=length_n(Prediction), Predictions=enlist(Prediction),
            Indicators_n=unique(Indicators_n), Indicators=unique(Indicators)) %>%
  ungroup()

#write.csv(Hydrophobin_groups2, "Tables/Hydrophobin_groups.csv", row.names=F)

Hyd_org=Hydrophobins %>%
  mutate(Grouped=ifelse(is.na(Cluster), "Ungrouped", "Grouped")) %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  group_by(Groups, Genome) %>%
  tally

method.plt=ggplot(Hyd_org, aes(Groups, Genome, fill=n))+geom_tile()+theme+scale_fill_viridis_c()+xlab("Search Method")+theme+theme(axis.title.y=element_blank(), strip.background=element_rect(fill="white"), strip.text=element_text(size=20), axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+labs(fill="Count")+geom_text(aes(label=n), color="white", size=7)
method.plt

seq.trim=ex_group %>%
  select(Group, Accession, Prediction, Trimmed) %>%
  rowwise() %>%
  mutate(Accession=if(Prediction=="Not"){paste("Unsecreted", Accession, sep="_")}else(paste("Trim", Prediction, Accession, sep="_"))) %>%
  filter(Group %in% Hyd_groups) %>%
  mutate(Seq=list(s2c(Trimmed)))

seq.un=ex_group %>%
  select(Group, Accession, Prediction, Trimmed) %>%
  rowwise() %>%
  filter((Accession %in% Hydrophobin_accessions) & is.na(Group)) %>%
  mutate(Accession=if(Prediction=="Not"){paste("Unsecreted", Accession, sep="_")}else(paste("Trim", Prediction, Accession, sep="_"))) %>%
  mutate(Seq=list(s2c(Trimmed)))

seq2=bind_rows(seq.trim, seq.un)
  
#for(i in seq.trim$Group){
#   seqinr::write.fasta(sequences=seq.trim[seq.trim$Group==i,]$Seq, names=seq.trim[seq.trim$Group==i,]$Accession, file.out=paste("Proteomes/grouping/Trimmed", i, "sequences.fa", sep=" "))
#}

#seqinr::write.fasta(sequences=seq2$Seq, names=seq2$Accession, file.out="All_strict_candidates.fa")

####8C Pattern Finding####
eight_C="C[^C]+CC[^C]+C[^C]+C[^C]+CC[^C]+C"

eight_dat = ex_group3 %>%
  rowwise() %>%
  mutate(C8=ifelse(grepl(eight_C, Trimmed), "Yes", "No"), C6=ifelse(Candidates.6C=="6C", "Yes", "No"), Strict=ifelse(Accession %in% Hydrophobin_accessions, "Yes", "No")) %>%
  arrange(desc(Strict), desc(`C6`), desc(`C8`))

ftable(eight_dat$C6~eight_dat$C8)

StrictC6notC8=eight_dat %>%
  filter(`C6`=="Yes", `C8`=="No", Strict=="Yes")

#write.csv(eight_dat, "Tables/eightC_dat.csv", row.names=F)

missing_stat= key.dat2 %>% mutate(Status=ifelse(Genome %in% unique(Hydrophobins$Genome), "Present", "Missing")) %>% group_by(Phylum, Status) %>% summarize(n=length(unique(Species)), list=toString(unique(Species))) %>%
  pivot_wider(id_cols=Phylum, names_from=Status, values_from=c(n, list)) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  mutate(n_Total=n_Missing+n_Present, Percent=n_Present/n_Total) %>%
  select(Phylum, n_Present, n_Missing, n_Total, Percent, list_Present, list_Missing)

#write.csv(missing_stat, "Tables/missing_stat.csv", row.names=F)

missing_stat_niche= key.dat2 %>%
  filter(Kingdom=="Fungi") %>% mutate(Status=ifelse(Genome %in% unique(Hydrophobins$Genome), "Present", "Missing")) %>% group_by(Niche, Status) %>% summarize(n=length(unique(Species)), list=toString(unique(Species))) %>%
  pivot_wider(id_cols=Niche, names_from=Status, values_from=c(n, list)) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  mutate(n_Total=n_Missing+n_Present, Percent=n_Present/n_Total) %>%
  select(Niche, n_Present, n_Missing, n_Total, Percent, list_Present, list_Missing)

missing_stat_style= key.dat2 %>%
  filter(Kingdom=="Fungi") %>% mutate(Status=ifelse(Genome %in% unique(Hydrophobins$Genome), "Present", "Missing")) %>% group_by(Lifestyle, Status) %>% summarize(n=length(unique(Species)), list=toString(unique(Species))) %>%
  pivot_wider(id_cols=Lifestyle, names_from=Status, values_from=c(n, list)) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  mutate(n_Total=n_Missing+n_Present, Percent=n_Present/n_Total) %>%
  select(Lifestyle, n_Present, n_Missing, n_Total, Percent, list_Present, list_Missing)

Hyd_group_count=Hydrophobins %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host) %>%
  mutate(tot=length(Accession)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host, Cluster) %>%
  summarize(count=length(Accession))

Hyd_group_per=Hydrophobins %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host) %>%
  mutate(tot=length(Accession)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host, Cluster) %>%
  summarize(per=length(Accession)/tot) %>%
  distinct()

Hyd_genome_tot=Hydrophobins %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host, Prediction) %>%
  summarize(count=length(Accession)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host) %>%
  mutate(tot=sum(count)) %>%
  rename("Secretion"=Prediction) %>%
  mutate(Secretion=as.character(Secretion), Secretion=replace(Secretion, Secretion=="Not", "Unsecreted"))

library(moments)
Hyd_genome_tot2=Hydrophobins %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  select(Genome, Cluster, Phylum, Prediction, Accession) %>%
  group_by(Genome, Cluster, Phylum, Prediction) %>%
  summarize(count=length(Accession)) %>%
  group_by(Genome, Phylum, Prediction) %>%
  summarize(tot=sum(count), clust_n=length(unique(Cluster)), skew=skewness(count), sd=sd(count), range=max(count)-min(count)) %>%
  rename("Secretion"=Prediction) %>%
  mutate(Secretion=as.character(Secretion), Secretion=replace(Secretion, Secretion=="Not", "Unsecreted"))

genome.mini.key=Hydrophobins %>%
  select(Genome, Lifestyle, Niche, Phylum, Host) %>%
  distinct() %>%
  left_join(key.dat2)

theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=25), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), title = element_text(size=35), legend.title = element_text(size=15), legend.text = element_text(size=10))

plt=ggplot(Hyd_genome_tot, aes(x=Genome, y=count, fill=Secretion))+theme+
  geom_col()+geom_text(data=(Hyd_genome_tot %>% select(Genome, tot) %>% distinct()), aes(x=Genome, y=tot+15, label = tot), fontface ="plain", size = 7, inherit.aes=F)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = c(0.9, 0.7))+ylab("Number")+scale_fill_viridis_d(option="G", direction = -1)+scale_y_continuous(limits=c(0, 155))
plt

plt2=ggplot(Hyd_group_per, aes(Genome, per, fill=Cluster))+theme+
  geom_col(aes(color=Phylum), size=2)+scale_fill_viridis_d()+
  theme(legend.position="bottom", legend.title=element_blank())+guides(fill=guide_legend(nrow=3,byrow=TRUE))+geom_text(data=Hyd_group_per[Hyd_group_per$per>0.05,], aes(label=str_replace(Cluster, "Group ", "")), col=rev(grey.colors(length(Hyd_group_per[Hyd_group_per$per>0.05,]$Cluster))), position = position_stack(vjust = 0.5), fontface ="plain", size = 7)+scale_y_continuous(labels = scales::percent)+ylab("Percent")+scale_color_viridis_d(option="turbo")+theme(axis.text.x=element_text(angle = 90, hjust=0.5), axis.title.x=element_blank())
plt2

#Grobbing from here: https://genviz.org/module-07-appendix/0007/01/01/advancedggplot2/

grob1=ggplotGrob(plt)
grob2=ggplotGrob(plt2)

maxWidth <- unit.pmax(grob1$widths, grob2$widths)

grob1$widths= maxWidth
grob2$widths= maxWidth

layout <- rbind(c(1),
                c(2))

grid.arrange(grob1, grob2, layout_matrix=layout, heights=c(0.3, 0.7))

HLPs = ex_group %>%
  select(Genome, Length, Indicators_n, Indicators, Cysteines, Accession, Phylum, Lifestyle, Prediction, Host, Niche, Groups, Group, Sequence, Pfam.Hyd1, Pfam.Hyd2) %>%
  group_by(Group) %>%
  mutate(Genom_n=length(unique(Genome))) %>%
  group_by(Accession, Phylum, Lifestyle, Host, Niche, Groups, Genome, Prediction, Length, Cysteines) %>%
  summarize(Cluster_n=length_n(Group), Cluster=enlist(Group),
            Indicators_n=unique(Indicators_n), Indicators=unique(Indicators),
            Hyd=length(grep("Hyd", Pfam.Hyd1, fixed=T))+length(grep("Hyd", Pfam.Hyd2, fixed=T))) %>%
  filter(!(Accession %in% Hydrophobins$Accession)) %>%
  ungroup()

missing_stat_HLPs= key.dat2 %>% mutate(Status=ifelse(Genome %in% unique(HLPs[!is.na(HLPs$Cluster),]$Genome), "Present", "Missing")) %>% group_by(Phylum, Status) %>% summarize(n=length(unique(Species)), list=toString(unique(Species))) %>%
  pivot_wider(id_cols=Phylum, names_from=Status, values_from=c(n, list)) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  mutate(n_Total=n_Missing+n_Present, Percent=n_Present/n_Total) %>%
  select(Phylum, n_Present, n_Missing, n_Total, Percent, list_Present, list_Missing)

HLPs_group_per=HLPs %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo[c(1:2, 50:51, 3:5, 52, 6:49)], levels(Genome))])) %>%
  filter(!is.na(Cluster)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host) %>%
  mutate(tot=length(Accession)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host, Cluster) %>%
  summarize(per=length(Accession)/tot) %>%
  distinct()

HLPs_genome_tot=HLPs %>%
  mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo[c(1:2, 50:51, 3:5, 52, 6:49)], levels(Genome))])) %>%
  filter(!is.na(Cluster)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host, Prediction) %>%
  summarize(count=length(Accession)) %>%
  group_by(Genome, Lifestyle, Niche, Phylum, Host) %>%
  mutate(tot=sum(count)) %>%
  rename("Secretion"=Prediction) %>%
  mutate(Secretion=as.character(Secretion), Secretion=replace(Secretion, Secretion=="Not", "Unsecreted"))

theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=25), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), title = element_text(size=35), legend.title = element_text(size=15), legend.text = element_text(size=10))

plt3=ggplot(HLPs_genome_tot, aes(x=Genome, y=count, fill=Secretion))+theme+
  geom_col()+geom_text(aes(x=Genome, y=tot+10, label = tot), fontface ="plain", size = 7)+scale_y_continuous(limits=c(0, 160))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = c(0.9, 0.7))+ylab("Number")+scale_fill_viridis_d(option="G")

plt3

plt4=ggplot(HLPs_group_per, aes(Genome, per, fill=Cluster))+theme+
  geom_col(aes(color=Phylum), size=2)+scale_fill_viridis_d(option="G")+
  theme(legend.position="bottom", legend.box="vertical", legend.title=element_blank())+guides(fill="none", color=guide_legend(order=1))+geom_text(aes(label=ifelse(per>0.15, str_replace(Cluster, "Group ", ""), NA)), size=6, position=position_stack(0.5), angle=90, color="grey50")+scale_y_continuous(labels = scales::percent)+scale_color_viridis_d(option="turbo")+theme(axis.text.x=element_text(angle = 90, vjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())
plt4

grob3=ggplotGrob(plt3)
grob4=ggplotGrob(plt4)

maxWidth <- unit.pmax(grob3$widths, grob4$widths)

grob3$widths= maxWidth
grob4$widths= maxWidth

layout <- rbind(c(1),
                c(2))

grid.arrange(grob3, grob4, layout_matrix=layout, heights=c(0.3, 0.7))

####Class determination####
many <- read.table("Proteomes/grouping/Hmms/out_strict_hyd.txt", quote="\"")

many2 = many %>%
  select(V1, V3, V6) %>%
  rename("Category"=V3, "score"=V6) %>%
  filter(Category!="AllStrict") %>%
  mutate(Accession=V1) %>%
  mutate(Accession=str_replace(Accession, "Trim_", "")) %>%
  separate(V1, into=c("Secretion", "Genome", "extra"), sep="_") %>%
  mutate(Accession=str_replace(Accession, ".*ecreted_", "")) %>%
  select(-extra) %>%
  group_by(Accession) %>%
  mutate(Max=max(score)) %>%
  filter(Max==score)

ex_class = Hydrophobins %>%
  select(Accession) %>%
  left_join(ex_group) %>%
  select(Accession, Group) %>%
  left_join(many2) %>%
  filter(!(Accession %in% ind$Indicator)) %>%
  mutate(Category=str_replace(Category, "up", "up ")) %>%
  select(Accession, Group, Category) %>%
  mutate(Match=Group==Category)

many3 = many %>%
  select(V1, V3, V5) %>%
  rename("Category"=V3, "score"=V5) %>%
  filter(Category!="AllStrict") %>%
  mutate(Accession=V1) %>%
  separate(V1, into=c("Secretion", "Genome", "extra"), sep="_") %>%
  select(-extra) %>%
  group_by(Category) %>%
  summarize(Accessions=list(Accession))

updat=many3$Accessions
names(updat)=many3$Category

upset(fromList(updat), order.by="freq", nsets=20, nintersects=29)

many4=data.frame(Names=names(updat))
for(i in 1:length(names(updat))){
  tmp=updat[i]
  array=NULL
  for(j in 1:length(names(updat))){
    array[j]=length(intersect(unlist(tmp), unlist(updat[j])))
  }
  many4=bind_cols(many4, array)
  colnames(many4)[i+1]=names(updat)[i]
}

many5=many4 %>%
  pivot_longer(-Names, names_to="Names2", values_to="Count") %>%
  group_by(Names2) %>%
  mutate(Total=max(Count)) %>%
  #filter(Names!=Names2) %>%
  rowwise() %>%
  mutate(Overlap=Count/Total)

ggplot(many5, aes(Names, Names2, fill=Overlap))+geom_tile()+theme+scale_fill_viridis_c(labels = scales::percent)+theme(axis.text.x=element_text(angle=90))

many.dist=as.matrix(many4[2:length(many4)])
names(many.dist)=many4$Names

many.hclust=hclust(d = dist(many.dist))

ord=many.hclust$order

strict_ind_group=ex_group2 %>%
  select(Group, Indicators_n, Phyla) %>%
  distinct() %>%
  mutate(Names2=str_replace(Group, " ", ""), Names=Names2, Indicators_n=as.character(Indicators_n)) %>%
  select(-Group) %>%
  rowwise() %>%
  mutate(Indicators_n=ifelse(Names=="Group82", "EAS", Indicators_n)) %>%
  mutate(Indicators_n=ifelse(Names=="Group1", "CU", Indicators_n)) %>%
  bind_rows(data.frame(Indicators_n=rep("Pfam", times=2), Names=c("Hydrophobin", "Hydrophobin_2"), Names2=c("Hydrophobin", "Hydrophobin_2"), Phyla=c("Ascomycota, Basidiomycota", "Ascomycota")))

many6 = many5 %>%
  left_join(strict_ind_group) %>%
  mutate(Names=as.factor(str_replace(Names, "up", "up ")), Names2=as.factor(str_replace(Names2, "up", "up "))) %>%
  mutate(Names=factor(Names, levels=levels(Names)[ord]),
         Names2=factor(Names2, levels=levels(Names2)[ord])) %>%
  rowwise() %>%
  mutate(Name_c=toString(sort(c(Names, Names2)))) %>%
  arrange(Names, Names2) %>%
  mutate(Overlap=ifelse(Overlap==0, NA, Overlap)) %>%
  mutate(Phyla=str_replace(Phyla, "Ascomycota, Basidiomycota", "Both")) %>%
  arrange(desc(is.na(Phyla)))

pair.plt=ggplot(many6, aes(Names2, Names, fill=Overlap, color=Phyla))+geom_tile(size=2, na.rm=T)+theme+scale_fill_viridis_c(labels = scales::percent, option="D")+theme(axis.text.x=element_text(angle=-45, vjust=0.5, hjust=0), legend.title=element_text(size=20), legend.text=element_text(size=20))+scale_y_discrete(position = "right")+scale_color_manual(values=viridis(length(unique(many6$Phyla)), option="C"), limits = c(unique(many6[!is.na(many6$Phyla),]$Phyla)), na.value="#000000")+
  geom_text(inherit.aes=F, aes(x=Names2, y=Names, label=Indicators_n), size=5)+xlab("")+ylab("")

ggplotly(pair.plt)

many7 = many6 %>%
  select(Names, Names2, Overlap) %>%
  mutate(Overlap=replace_na(Overlap, 0), Overlap=round(Overlap*100), Overlap=ifelse(Names==Names2, NA, Overlap)) %>%
  rename("X"=Names) %>%
  pivot_wider(names_from=Names2, values_from=Overlap)

dist3=match_thresh(many7, thresh=100)

last.100=unionize(dist3)

out4=collapse(last.100)

out5=data.frame(Class=names(out4), Values=unlist(as.character(out4))) %>%
  mutate(Values=str_replace_all(Values, "p ", "p_")) %>%
  separate(Values, into=as.character(1:1000), sep="[^A-z0-9.-]") %>%
  select(Class, c(as.character(2:1000))) %>%
  pivot_longer(cols=-Class, names_to="col", values_to="Accession") %>%
  drop_na(Accession) %>%
  filter(Accession!=" " & Accession!="") %>%
  rename("Group"=Accession) %>%
  select(-col) %>%
  mutate(Group=str_replace_all(Group, "p_", "p "),
         Class=str_replace_all(Class, "Group", "Class "))

many8= many7 %>%
  select(X) %>%
  rename("Group"=X) %>%
  merge(out5, all.x=T) %>%
  mutate(Class=str_replace(Class, "3", "4"),
         Class=str_replace(Class, "2", "3"),
         Class=str_replace(Class, "4", "2"),
         Class=replace_na(Class, "Class 4"),
         Class=str_replace(Class, "  ", " "))

Hyd_group_per2=Hydrophobins %>%
  dplyr::mutate(Genome=as.factor(Genome), Genome=factor(Genome, levels=levels(Genome)[match(phylo, levels(Genome))])) %>%
  dplyr::group_by(Genome, Lifestyle, Niche, Phylum, Host,) %>%
  dplyr::mutate(tot=length(Accession)) %>%
  dplyr::group_by(Genome, Lifestyle, Niche, Phylum, Host, Cluster) %>%
  mutate(Cluster=ifelse(is.na(Cluster), ifelse(grepl("Hyd1", Groups), "I", "II"), Cluster)) %>%
  dplyr::summarize(per=length(Accession)/tot) %>%
  distinct() %>%
  dplyr::rename(Group="Cluster") %>%
  left_join(many8) %>%
  dplyr::mutate(CG=str_c(Class, Group, sep="; "))  %>%
  dplyr::arrange(Class, Group)

Hyd_with_class = eight_dat %>%
  left_join(many8) %>%
  dplyr::mutate(CG=str_c(Class, Group, sep=";"))  %>%
  dplyr::arrange(Class, Group) %>%
  dplyr::mutate(CG=str_replace(CG, ";", "\n")) %>%
  dplyr::mutate(CG=str_replace(CG, "Class 4", "No Class")) %>%
  dplyr::mutate(CG=replace_na(CG, "No Class\nNo Group")) %>%
  dplyr::mutate(Class=str_replace(Class, "Class 4", "NA")) %>%
  filter(Strict=="Yes")

#write.csv(Hyd_with_class, "Tables/Strict_hyd_with_class.csv", row.names=F)

plt2_class=ggplot(Hyd_group_per2, aes(Genome, per, fill=Class))+theme+
  geom_col(aes(color=Phylum), size=2)+scale_fill_viridis_d()+
  theme(legend.position="bottom", legend.title=element_blank())+guides(fill=guide_legend(nrow=3,byrow=TRUE))+geom_text(data=Hyd_group_per2[Hyd_group_per2$per>0.05,], aes(label=str_replace(Class, "Class ", "")), col=rev(grey.colors(length(Hyd_group_per2[Hyd_group_per2$per>0.05,]$Class))), position = position_stack(vjust = 0.5), fontface ="plain", size = 7)+scale_y_continuous(labels = scales::percent)+ylab("Percent")+scale_color_viridis_d(option="turbo")+theme(axis.text.x=element_text(angle = 90, hjust=0.5), axis.title.x=element_blank())

plt2_class

find_hue2=function(var){
  len=length(unique(var))
  options=3*(2^(ceiling(log(plyr::round_any(len, 3, ceiling)/3, base=2))))
  ratio=options/3
  hues=as.numeric()
  for(rnd in 0:(ratio-1)){
    start=rnd*ratio+rnd
    for(i in 0:2){
      x=start+ratio*i
      while(x>options){
        x=x-options
      }
      hues=c(hues, x)
    }
  }
  return(hues[1:len]/options)
}

nest_colors_shuffle2=function(Base, Fade, Offset=0){
  dat=data.frame(First=Base, Second=Fade)
  dat2=dat %>%
    distinct() %>%
    left_join(data.frame(First=unique(Base), h=find_hue2(unique(Base)))) %>%
    dplyr::group_by(First) %>%
    dplyr::mutate(v=rev(as.numeric(as.factor(Second)))/length(unique(Second)),
                  hsv=hsv(h+Offset-floor(h+Offset), (1-v)*.8+0.2, v*.8+0.2))
  return(dat2)  
}

pal3=nest_colors_shuffle2(Hyd_group_per2[!is.na(Hyd_group_per2$Class),]$Class, Hyd_group_per2[!is.na(Hyd_group_per2$Class),]$Group, Offset=0.2) %>%
  dplyr::mutate(CG=ifelse(!grepl(";", Second), str_c(First, Second, sep="; "), Second)) %>% ungroup()

Hyd_group_per3 = Hyd_group_per2 %>%
  dplyr::mutate(CG=str_replace(CG, ";", "\n")) %>%
  dplyr::mutate(CG=str_replace(CG, "Class 4", "No Class")) %>%
  dplyr::mutate(CG=replace_na(CG, "No Class\nNo Group"))

plt2_CG=ggplot(Hyd_group_per3, aes(Genome, per, fill=CG))+theme+
  geom_col(aes(color=Phylum), size=2)+scale_fill_manual(values=c(pal3$hsv, "grey"))+
  theme(legend.position="bottom", legend.title=element_blank())+guides(fill=guide_legend(nrow=3,byrow=TRUE), color=guide_legend(nrow=2, byrow=T))+geom_text(data=Hyd_group_per3[Hyd_group_per3$per>0.05,], aes(label=str_replace(Group, "Group ", "")), col="black", position = position_stack(vjust = 0.5), fontface ="plain", size = 7)+scale_y_continuous(labels = scales::percent)+ylab("Percent")+scale_color_viridis_d(option="turbo")+theme(axis.text.x=element_text(angle = 90, hjust=0.5, size=20), axis.text.y=element_text(size=20), axis.title.x=element_blank())
plt2_CG

SSTsub=Hyd_group_per3 %>%
  filter(Genome=="SST" & (Group!="Group 10" | is.na(Group)))

pal4=c(pal3$hsv, "grey")
names(pal4)=unique(Hyd_group_per3$CG)

plt_SSTsub=ggplot(SSTsub, aes(Genome, per, fill=CG))+theme+
  geom_col(aes(color=Phylum), size=2)+scale_fill_manual(values=pal4)+
  theme(legend.position="bottom", legend.title=element_blank())+guides(fill=guide_legend(nrow=3,byrow=TRUE))+geom_text(data=SSTsub, aes(label=str_replace(Group, "Group ", "")), col="black", position = position_stack(vjust = 0.5), fontface ="plain", size = 7)+ scale_y_continuous(labels = scales::label_percent(accuracy = 1L))+scale_color_manual(values=viridis(n=4, option="turbo")[2])+theme(axis.text.x=element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=20), legend.position="none")+coord_flip()

plt_SSTsub

grob1=ggplotGrob(plt)
grob2=ggplotGrob(plt2_CG)
grob3=ggplotGrob(plt_SSTsub)

maxWidth <- unit.pmax(grob1$widths, grob2$widths)

grob1$widths= maxWidth
grob2$widths= maxWidth
grob3$widths= maxWidth

layout <- rbind(c(1),
                c(2),
                c(3))

grid.arrange(grob1, grob2, grob3, layout_matrix=layout, heights=c(0.2, 0.7, 0.1))

####PCA####
ind.seq=read.csv("Input/Indicator Sequences.tsv", sep="\t")

indicators2=ind %>%
  mutate(ID=paste("Secreted","_", Indicator, sep=""), Group2=Group) %>%
  bind_rows(data.frame(Indicator="AFU_RodE_Trim", Group="Group 39", ID="Unsecreted_AFU_RodE", Group2="Group 39")) %>%
  left_join(ind.seq) %>%
  select(Indicator, ID, Group2) %>%
  rename(Accession="ID")

dist.pca=dist %>%
  mutate(Accession=X) %>%
  column_to_rownames('X') %>%
  mutate(Accession=str_trim(Accession, "right")) %>%
  replace(is.na(.), 100) %>%
  left_join(indicators2) %>%
  mutate(Accession=str_replace(Accession, ".*ecreted_", "")) %>%
  mutate(Accession=str_replace_all(Accession, "\\.", "-")) %>%
  left_join(ex_group3) %>%
  left_join(indicators2) %>%
  mutate(Strict=ifelse(Accession %in% Hydrophobin_accessions | Accession %in% str_replace(indicators2$Accession, ".*ecreted_", ""),
                       "Strict", "Not")) %>%
  mutate(Group=ifelse(!is.na(Group2), Group2, Group), Prediction=ifelse(!is.na(Group2), "Indicator", Prediction)) %>%
  mutate(Phylum=replace_na(Phylum, "Indicator")) %>%
  filter(Strict=="Strict")

pca.dist=dist.pca %>%
  select(1:1263) %>%
  prcomp()

#Adjusting guides and indicator information
pal_named=c(pal3$hsv, "white", "red")
names(pal_named)=str_replace(c(pal3$CG, "No Class; No Group", "Indicator"), "Class 4", "No Class")

dist.pca.plt=pca.dist %>%
  augment(dist.pca) %>%
  left_join(Hyd_group_per2) %>%
  mutate(CG=ifelse(!is.na(Indicator), "Indicator", CG)) %>%
  mutate(CG=ifelse(is.na(CG), "No Class; No Group", CG)) %>%
  mutate(CG=str_replace(CG, "Class 4", "No Class")) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color=CG, label=Accession, text=Group, shape=Phylum))+
  geom_point(size = 3) + ggrepel::geom_label_repel(aes(label=Indicator), fill= alpha(c("white"),0.5), segment.color="white", color="black", segment.curvature=0.3, box.padding = 0.5)+ background_grid()+guides(label="none")+theme_dark()+scale_color_manual(values=pal_named, na.value = "white")+theme(legend.text=element_text(color="white", size=12), legend.background=element_rect(fill="grey50", color="white"), legend.position = c(0.62, 0.62))+guides(shape=guide_legend("Phyla", nrow=1, order=2, color="white", title.position="left"), color=guide_legend("Class; Group", ncol=1, order=1))+theme(text=element_text(size=20), legend.box = "horizontal",  legend.box.just="bottom", legend.title=element_text(color="white"))+xlab("PC1")+ylab("PC2")

dist.pca.plt2=pca.dist %>%
  augment(dist.pca) %>%
  left_join(Hyd_group_per2) %>%
  mutate(CG=ifelse(!is.na(Indicator), "Indicator", CG)) %>%
  filter(.fittedPC1<(-200) & .fittedPC2<0) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color=CG, text=Group, shape=Phylum))+
  geom_point(size = 3) + ggrepel::geom_label_repel(aes(label=Indicator), fill= alpha(c("white"),0.5), segment.color="white", color="black", segment.curvature=0.3, box.padding = 0.5)+ background_grid()+guides(label="none")+theme_dark()+scale_color_manual(values=pal_named, na.value = "white")+theme(legend.text=element_text(color="white", size=15), legend.background=element_rect(fill="grey50", color="white"))+guides(color="none", shape="none")+theme(text=element_text(size=20))+xlab("PC1")+ylab("PC2")

dist.pca.plt+inset_element(dist.pca.plt2, 0.39, 0.335, 1, 1)

#Niche PCA
niche.pca.plt=pca.dist %>%
  augment(dist.pca) %>%
  mutate(Host=replace_na(Host, "None")) %>%
  mutate(Niche=replace_na(Niche, "Indicator")) %>%
  left_join(Hyd_group_per2) %>%
  mutate(CG=ifelse(!is.na(Indicator), "Indicator", CG)) %>%
  mutate(CG=ifelse(is.na(CG), "No Class; No Group", CG)) %>%
  mutate(CG=str_replace(CG, "Class 4", "No Class")) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color=Niche, label=Accession, text=Group, shape=Host))+
  geom_point(size = 3) + ggrepel::geom_label_repel(aes(label=Indicator), fill= alpha(c("white"),0.5), segment.color="white", color="black", segment.curvature=0.3, box.padding = 0.5)+ background_grid()+guides(color=guide_legend(ncol=1, order=1), label="none")+theme_dark()+scale_color_manual(values = c("red", "blue", "green", "yellow"))+theme(legend.text=element_text(color="white", size=12), legend.background=element_rect(fill="grey50", color="white"), legend.position = c(0.51, 0.38))+guides(shape=guide_legend("Host", nrow=1, order=2, color="white", title.position="left"))+theme(text=element_text(size=20), legend.box = "horizontal",  legend.box.just="bottom", legend.title=element_text(color="white"))+xlab("PC1")+ylab("PC2")

ggplotly(niche.pca.plt)

niche.pca.plt2=pca.dist %>%
  augment(dist.pca) %>%
  mutate(Host=replace_na(Host, "None")) %>%
  mutate(Niche=replace_na(Niche, "Indicator")) %>%
  left_join(Hyd_group_per2) %>%
  mutate(CG=ifelse(!is.na(Indicator), "Indicator", CG)) %>%
  filter(.fittedPC1<(-200) & .fittedPC2<0) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color=Niche, text=Group, shape=Host))+
  geom_point(size = 3) + ggrepel::geom_label_repel(aes(label=Indicator), fill= alpha(c("white"),0.5), segment.color="white", color="black", segment.curvature=0.3, box.padding = 0.5)+ background_grid()+guides(label="none")+theme_dark()+scale_color_manual(values = c("red", "blue", "green", "yellow"))+theme(legend.text=element_text(color="white", size=15), legend.background=element_rect(fill="grey50", color="white"))+guides(color="none", shape="none")+theme(text=element_text(size=20))+xlab("PC1")+ylab("PC2")

niche.pca.plt+inset_element(niche.pca.plt2, 0.39, 0.335, 1, 1)

#Centered, nested proportional area charts
#Data from Mulder and Wessels 1986
sc_dat <- read.csv("Input/SC hydrophobin expression data.csv") %>%
  pivot_longer(cols=c("Expression", "Total"), names_to="Category", values_to="Values") %>%
  mutate(Category=as.factor(Category)) %>%
  mutate(Category=factor(Category, levels=c("Total", "Expression")))

ggplot() +
  geom_circle(aes(x0 = 1, y0 = 1, r = sqrt(Values), fill = Category), data = sc_dat)+theme_void()+facet_grid(cols=vars(Hydrophobin), rows=vars(Type), switch="y")+coord_fixed()+geom_label(x=sqrt(40), y=1, aes(label=paste(Values, "%", sep="")), size=10, data=sc_dat %>% filter(Category=="Expression"))+theme(text=element_text(size=40), legend.title=element_blank(), legend.position="bottom", strip.text.y=element_text(margin=margin(r=10)))+scale_fill_viridis_d()
