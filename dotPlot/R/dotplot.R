#' This functions takes all pairwise alignments and creates a dot matrix
#' representing the alignment's path.
#' @param filepath: either a directory contains all alignments (DNA/AA), or a json file contains all alignments of a single raw DNA/AA sequence.
#' @return a single dot matrix.
#' @export
generate_dot_matrix = function(input){
  #read input
  aln = read_input(input)

  #initialize matrix
  mat = 0
  for(i in seq_len(nrow(aln))) {
    mat <- mat + dot_matrix(aln$Seq1[i], aln$Seq2[i], aln$n[i])
  }
  return(mat)
}


#' This function plots a dot matrix as a raster image.
#'
#' @param  filepath:    either a directory contains all alignments (DNA/AA), or a json file contains all alignments of a single raw DNA/AA sequence.
#' @param  use_ggplot: `TRUE` or `FALSE`. It is false by default. Otherwise, setup use_ggplot=TRUE which returns a raster ggplot.
#' @return A raster plot.
#' @export
plot_dot_matrix = function(input,use_ggplot=FALSE){

  aln = read_input(input)
  mat = generate_dot_matrix(input)
  create_plot(mat, aln, use_ggplot)
}

#' Produce a data frame given a dot matrix. The dot matrix can be gnerated by generate_dot_matrix(). 
#'
#' @param mat: dot matrix.
#' @return a data frame (Within the dataframe, The integeter values in either columns indiciate where the real match/gap occurs, the non-integer values indicates the locus of potential insertion space.)
#' @export
matrix_to_dataframe <- function(mat) {
  nrows = nrow(mat)
  ncols = ncol(mat)
  len = length(mat)

  #coordinates for seq1 & seq2
  #seq1: increments of 0.5, each repeated by number of columns
  #seq2: increments of 0.5, repeated overall by number of rows
  seq1 = seq(from=0.5, by=0.5, length.out=nrows)
  seq1 = rep(seq1, each=ncols)
  seq2 = seq(from=0.5, by=0.5, length.out=ncols)
  seq2 = rep(seq2, times=nrows)

  #magnitude
  magnitude = c(mat)

  #data frame
  data = data.frame(Seq1=seq1, Seq2=seq2, Magnitude=magnitude)
  data
}


################################################################################
#Helper functions
#convert a single alignment to a dot matix
dot_matrix <- function(a,b,n=1) {
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

#plot dot matrix
create_plot = function(mat,aln,use_ggplot) {
  if(!use_ggplot) {
    plot(as.raster(1-mat/sum(aln$n)), interpolate = FALSE)
  }else{
    dm  = dim(mat)
    am  = 1:dm[2]
    bm  = dm[1]:1
    abm = expand.grid(am,bm)
    abm = cbind(abm,c(t(mat)))
    dmat = data.frame(abm)
    colnames(dmat)=c('Seq1','Seq2','Magnitude')

    ggplot(dmat)+ geom_raster(mapping=aes(x=Seq1, y=Seq2, fill=Magnitude)) +
      scale_fill_gradientn(colours= rev(terrain.colors(10)), name='Magnitude')
   }

}


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>reading funcs
#read the DNA/AA alignments from a file or directory
read_input = function(input_list) {
  #define data frame
  data = data.frame(Seq1 = character(), Seq2 = character())

  #if only one argument and is a dir, look for all fasta files inside dir
  if(dir.exists(input_list)) {
    fasta_list = read_dir(input_list)
    if(length(fasta_list) > 0) {
      input_list = fasta_list
    } else {
      stop("No fasta files found in current directory.")
    }
  }

  #read json files
  for(input in input_list) {
    if(tools::file_ext(input) == "fasta" || tools::file_ext(input) == "fa") {
      # read fasta, concatenate to data frame
      data[nrow(data)+1, ] = read_fasta(input)
    } else if(tools::file_ext(input) == "json") {
      # read json, concatenate to data frame
      data = rbind(data, read_json(input))
    } else {
      stop("Input format not supported.")
    }
  }

  #count unique alignments and return data frame
  data = data %>% dplyr::count(Seq1, Seq2)
  data
}

#find all fasta files in `dir_path` and return a list of files
read_dir = function(dir_path) {
  files = list.files(path = dir_path, pattern = "\\.fasta$", full.names = TRUE)
  files = append(files, list.files(path = dir_path, pattern = "\\.fa$", full.names = TRUE))
  files
}

#read json file (coati-sample output format) and return as data frame
read_json = function(json_path) {
    data = jsonlite::fromJSON(json_path)
    data = data %>% tidyr::unpack(aln)
    data
}

#read fasta file and return as vector
read_fasta = function(fasta_path) {
  data = data.frame(Seq1 = character(), Seq2 = character(), weight = double(), log_weight = double())
  fasta = seqinr::read.fasta(fasta_path, set.attributes = FALSE)
  seq1 =  paste(stringr::str_to_upper(fasta[[1]]), collapse = "")
  seq2 =  paste(stringr::str_to_upper(fasta[[2]]), collapse= "")
  # follow data frame format: Seq1, Seq2, weight, log_weight
  c(seq1, seq2, NA, NA)
}

