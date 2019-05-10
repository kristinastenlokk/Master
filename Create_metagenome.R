library("microseq")

Genome_names <- c("Brevundimonas_abyssalis","Stenotrophomonas_acidaminiphila","Brachybacterium_alimentarium","Acinetobacter_apis","Methylobacterium_aquaticum","Microbacteriaceae_bacterium.fna","Micrococcaceae_bacterium_C1-50","Sphingomonas_adhaesiva","Sphingobacterium_cellulitidis","Massilia_alkalitolerans","Bacillaceae_bacterium_B16-10","Janthinobacterium_agaricidamnosum","Human")
matrisene <- as.character(c("Brevundimonas_abyssalis_F1","Stenotrophomonas_acidaminiphila_F1","Brachybacterium_alimentarium_F1","Acinetobacter_apis_F1","Methylobacterium_aquaticum_F1","Microbacteriaceae_bacterium.fna_F1","Micrococcaceae_bacterium_C1-50_F1","Sphingomonas_adhaesiva_F1","Sphingobacterium_cellulitidis_F1","Massilia_alkalitolerans_F1","Bacillaceae_bacterium_B16-10_F1","Janthinobacterium_agaricidamnosum_F1","Human_F1"))

abundance <- as.vector(c(523459,201840,66421,60534,27047,25575,15087,10304,9016,8648,8464,6808,36799))
names(abundance) <- Genome_names
Genomes.list <- as.list(0)


#Laster inn alle fastq_f1-filene
for(i in 1:length(Genome_names)){
  Genomes.list[[i]] <- readFastq(paste0("Silico_air/",Genome_names[i],"1.fq",sep= ""))
}

#Setter sammen etter abundance:

Newgenome <- "" %>% data.frame

#Velger ut tilfeldige sekvenser (etter antall "phylotyper" (Abundance)), og limer dette sammen til en dataframe
#Må bruke set.seed for at reads-parene skal kunne slås sammen

for(i in 1:length(Genome_names)){
  set.seed(1234)
  idx <- base::sample(nrow((Genomes.list[[i]])), size = abundance[i], replace = T)
  y <- Genomes.list[[i]][idx,]
  Newgenome <- bind_rows(Newgenome,y)
}


#Trikser litt for å få i "rett" fastq-format. Ny header-linje som starter med +
SeqHeader <- as.data.frame(gsub("@", "+", Newgenome$Header))
Final <- bind_cols(SeqHeader,Newgenome)
Final <- Final[-1,c(3,4,1,5)]
colnames(Final) <- c("Header","Sequences","QualityHeader","Quality")



#Samme med fastq_f2-filene:

Genome_names <- c("Brevundimonas_abyssalis","Stenotrophomonas_acidaminiphila","Brachybacterium_alimentarium","Acinetobacter_apis","Methylobacterium_aquaticum","Microbacteriaceae_bacterium.fna","Micrococcaceae_bacterium_C1-50","Sphingomonas_adhaesiva","Sphingobacterium_cellulitidis","Massilia_alkalitolerans","Bacillaceae_bacterium_B16-10","Janthinobacterium_agaricidamnosum","Human")
matrisene <- as.character(c("Brevundimonas_abyssalis_F2","Stenotrophomonas_acidaminiphila_F2","Brachybacterium_alimentarium_F2","Acinetobacter_apis_F2","Methylobacterium_aquaticum_F2","Microbacteriaceae_bacterium.fna_F2","Micrococcaceae_bacterium_C1-50_F2","Sphingomonas_adhaesiva_F2","Sphingobacterium_cellulitidis_F2","Massilia_alkalitolerans_F2","Bacillaceae_bacterium_B16-10_F2","Janthinobacterium_agaricidamnosum_F2","Human_F2"))



#Abundance i samme rekkefølge som Genome_names
tot.nr.of.reads <- 1000000
parts <- c(2845,1097,361,329,147,139,82,56,49,47,46,372,276)
sum.parts <- (sum(parts))
perc <- (parts/sum.parts)

names(abundance) <- Genome_names
Genomes.list <- as.list(0)

#Laster inn alle fastq_f2-filene
for(i in 1:length(Genome_names)){
  Genomes.list[[i]] <- readFastq(paste0("Silico_air/",Genome_names[i],"2.fq",sep= ""))
}

#Setter sammen etter abundance:

Newgenome <- "" %>% data.frame

#Velger ut tilfeldige sekvenser (etter antall "phylotyper" (Abundance)), og limer dette sammen til en dataframe
#Må bruke set.seed for at reads-parene skal kunne slås sammen

for(i in 1:length(Genome_names)){
  set.seed(1234)
  idx <- base::sample(nrow((Genomes.list[[i]])), size = abundance[i], replace = T)
  y <- Genomes.list[[i]][idx,]
  Newgenome <- bind_rows(Newgenome,y)
}


#Trikser litt for å få i "rett" fastq-format. Ny header-linje som starter med +
SeqHeader <- as.data.frame(gsub("@", "+", Newgenome$Header))
Final2 <- bind_cols(SeqHeader,Newgenome)
Final2 <- Final2[-1,c(3,4,1,5)]
colnames(Final2) <- c("Header","Sequences","QualityHeader","Quality")


fjern1 <- str_detect(Final$Sequences, "[^ACGT]")
fjern2 <- str_detect(Final2$Sequences, "[^ACGT]")

m <- which(fjern1)
n <- which(fjern2)
p <- unique(c(m,n))

Final <- Final[-p,]
Final2 <- Final2[-p,]

colnames(Final) <- c("Header","Sequence","QualityHeader","Quality")   
colnames(Final2) <- c("Header","Sequence","QualityHeader","Quality")   

writeFastq(Final, "Pos_silico_dataset_R1.fastq")
writeFastq(Final2, "Pos_silico_dataset_R2.fastq")

