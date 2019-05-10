

shuffle.function <- function(input.file, output.file, w){
  library("microseq")
  library("gtools")
  library("microseq")
  
  d <- readFasta(input.file)
  all <- paste(d$Sequence, collapse = "")
  vec <- 1:(as.integer(nchar(all)/w)+1)
  thevector <- as.vector(0)
  
  for(i in vec){
    thevector[i] <- substring(all, (i-1)*w+1, (i)*w)
  }
  
  permuted <- permute(thevector)
  
  glued <- paste0(permuted, collapse = "")
  new.header <- paste0(">",d$Header[1],"_shuffled", collapse = "")
  new.data <- data.frame(Header = new.header, Secuene = glued)
  write.table(new.data, output.file, sep = "\n", col.names = F, quote = F, row.names = F)
  
}
