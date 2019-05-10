


artmachine <- function(ref.file, out.dir = ".", sampleID = "",
                       sequencing.depth = 10, N.readpairs = NULL,
                       sequencing.technology = "MSv3", read.length = 250,
                       error.free = TRUE,
                       paired.end = TRUE, insert.size = 750, insert.size.std = 100){
  
  paired <- if_else(paired.end, "--paired", "")
  err.free <- if_else(error.free, "--errfree", "")
  if(is.null(N.readpairs)){
    fold.coverage <- sequencing.depth
  } else {
    microseq::readFasta(ref.file) %>% 
      mutate(Length = str_length(Sequence)) -> rf
    fold.coverage <- (N.readpairs * 2 * read.length)/sum(rf$Length)
  }
  if(str_length(sampleID) > 0){
    sid <- paste("--id", sampleID)
    prefix <- sampleID
  } else{
    sid <- ""
    prefix <- "artsim"
  }
  
  cmd <- paste("art_illumina",
               "--rndSeed", sample(1:1000000, 1),
               "--quiet",
               "--noALN",
               paired,
               err.free,
               sid,
               "--seqSys", sequencing.technology,
               "--mflen", insert.size,
               "--sdev", insert.size.std,
               "--fcov", fold.coverage,
               "--len", read.length,
               "--in",  ref.file,
               "--out", file.path(out.dir, prefix))
  system(cmd)
  
  if(error.free){
    cmd <- paste("samtools view -bS",
                 file.path(out.dir, paste0(prefix, "_errFree.sam")),
                 ">",
                 file.path(out.dir, paste0(prefix, "_errFree.bam")))
    system(cmd)
    cmd <- paste("samtools fastq",
                 "-1", file.path(out.dir, paste0(prefix, "1_EF.fq")),
                 "-2", file.path(out.dir, paste0(prefix, "2_EF.fq")),
                 file.path(out.dir, paste0(prefix, "_errFree.bam")))
    system(cmd)
    ok <- file.remove(file.path(out.dir, paste0(prefix, "_errFree.bam")))
    ok <- file.remove(file.path(out.dir, paste0(prefix, "_errFree.sam")))
    ok <- file.remove(file.path(out.dir, paste0(prefix, ".sam")))
  }
  return(file.path(out.dir, prefix))
}



