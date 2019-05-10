library("microseq")
library("gtools")

source("shufflefunction.R")


names <- as.character(c("Brevundimonas_abyssalis",
                            "Stenotrophomonas_acidaminiphila",
                            "Brachybacterium_alimentarium",
                            "Acinetobacter_apis",
                            "Methylobacterium_aquaticum",
                            "Microbacteriaceae_bacterium",
                            "Micrococcaceae_bacterium_C1-50",
                            "Sphingomonas_adhaesiva",
                            "Sphingobacterium_cellulitidis",
                            "Massilia_alkalitolerans",
                            "Bacillaceae_bacterium_B16-10",
                            "Janthinobacterium_agaricidamnosum",
                            "Human"))
input <- as.vector(0)
output <- as.vector(0)

for(i in 1:length(names)){
  input[i] <- paste0("sim_dir/",names[i],".fna", collapse = "")
  output[i] <- paste0("Final_shuffled_genomes_25_50/",names[i],"_shuffled_25.fasta", collapse = "")
}

for(i in 1:length(names)){
  shuffle.function(input[i],output[i],25)
}


input2 <- as.vector(0)
output2 <- as.vector(0)

for(i in 1:length(names)){
  input2[i] <- paste0("sim_dir/",names[i],".fna", collapse = "")
  output2[i] <- paste0("Final_shuffled_genomes_25_50/",names[i],"_shuffled_50.fasta", collapse = "")
}

for(i in 1:length(names)){
  shuffle.function(input2[i],output2[i],50)
}
