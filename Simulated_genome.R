library(dplyr)
library(stringr)
library(RLinuxModules)
moduleInit(modulesHome="/local/genome/Modules/3.2.10") 
module("load art")
module("load samtools")
source("ART/artFun.R")
library("microseq")


#ref.file <- "./sim_dir/Acinetobacter_apis.fna"
#sampleID <- "Acinetobacter_apis"
out.dir <- "Silico_air/"

#Ref-filene er referansegenomer lastet ned fra NCBI sin RefSeq:


Genome_names <- c("Brevundimonas_abyssalis","Stenotrophomonas_acidaminiphila","Brachybacterium_alimentarium","Acinetobacter_apis","Methylobacterium_aquaticum","Microbacteriaceae_bacterium.fna","Micrococcaceae_bacterium_C1-50","Sphingomonas_adhaesiva","Sphingobacterium_cellulitidis","Massilia_alkalitolerans","Bacillaceae_bacterium_B16-10","Janthinobacterium_agaricidamnosum")
ref.file <- c("./sim_dir/Brevundimonas_abyssalis.fna","./sim_dir/Stenotrophomonas_acidaminiphila.fna","./sim_dir/Brachybacterium_alimentarium.fna","./sim_dir/Acinetobacter_apis.fna","./sim_dir/Methylobacterium_aquaticum.fna","./sim_dir/Microbacteriaceae_bacterium.fna","./sim_dir/Micrococcaceae_bacterium_C1-50.fna","./sim_dir/Sphingomonas_adhaesiva.fna","./sim_dir/Sphingobacterium_cellulitidis.fna","./sim_dir/Massilia_alkalitolerans.fna","./sim_dir/Bacillaceae_bacterium_B16-10.fna","./sim_dir/Janthinobacterium_agaricidamnosum.fna")

#Den mest abundant arten er simulert med flere reads for å ha nok å trekke fra

ff <- artmachine(ref.file[1], out.dir = out.dir, sampleID = Genome_names[1],
                 sequencing.depth = 10, N.readpairs = 1000000,
                 sequencing.technology = "MSv3", read.length = 250,
                 error.free = F)

#Resten er simulet med halvparten av readsene

for(i in 2:length(Genome_names)){
  ff <- artmachine(ref.file[i], out.dir = out.dir, sampleID = Genome_names[i],
                   sequencing.depth = 10, N.readpairs = 500000,
                   sequencing.technology = "MSv3", read.length = 250,
                   error.free = F)
}



