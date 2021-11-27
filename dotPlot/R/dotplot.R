


##############################################
# input="sample.json"
# input="samples"



#' Plot a dot matrix.  
#' This function plots a dot matrix as a raster image.  
#'  
#' @param  filepath: either a directory contains all alignments (DNA/AA), or a json file
#' @param  isDir: `TRUE` or `FALSE`. It is true by default that our filepath is a directory. otherwise, it is a json format.
#' @return A raster plot
#' @export
plot_dot_matrix = function(input,isDir=T){
  
  #generate alignment from a folder 
  aln_from_dir <- function(inD){
    Files = list.files(inD,full.names=T)
    seqD  = data.frame()
    for(i in 1:length(Files)){
      seqs = Biostrings::readBStringSet(Files[i],use.names=FALSE)
      seqV = as.vector(seqs)
      seqD[i,'Seq1'] = seqV[1]
      seqD[i,'Seq2'] = seqV[2]
    }
    aln = seqD %>% count(Seq1,Seq2)
    return(aln)
  }
  
  #generate alignment from a json file
  aln_from_json <- function(inF){
    dat <- jsonlite::fromJSON(inF)
    dat <- dat %>% unpack(aln)
    aln <- dat %>% count(Seq1,Seq2)
    return(aln)
  }
  
  dot_matrix <- function(a, b, n=1) {
    ab  <- stringr::str_split(c(a,b),"",simplify=TRUE)
    ab  <- t(ab != "-")
    len <- nrow(ab) - colSums(!ab)
    m   <- matrix(0, nrow = 2*len[1]+1, ncol = 2*len[2]+1)
    ai  <- 1
    bi  <- 1
    for(i in seq_len(nrow(ab))) {
      m[ai,bi] <- n
      if(ab[i,1] == TRUE && ab[i,2] == TRUE) {#match
        m[ai+1,bi+1] <- n
        ai <- ai+2
        bi <- bi+2
      } else if(ab[i,1] == TRUE) {#deletion
        m[ai+1,bi] <- n
        ai <- ai+2
      } else if(ab[i,2] == TRUE) {#insertion
        m[ai,bi+1] <- n
        bi <- bi+2
      }
    }
    m[ai,bi] <- n
    m
  }
  
  
  if(isDir){
    aln = aln_from_dir(input)
  }else{
    aln = aln_from_json(input)
    }
  
  mat <- 0
  for(i in seq_len(nrow(aln))) {
    mat <- mat + dot_matrix(aln$Seq1[i], aln$Seq2[i], aln$n[i])
  }
  plot(as.raster(1-mat/sum(aln$n)), interpolate=FALSE)
}




###################
# body <- list()
# innerBody <- list()
# innerBody$a <- "A"
# innerBody$b <- "B"
# body$c <- innerBody
# jsonlite::toJSON(body, pretty=TRUE, auto_unbox=TRUE)







