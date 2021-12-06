#' Read input alignments.
#' This function reads DNA/AA alignments from a file or directory.
#'
#' @param file_path: path to directory containing fasta alignments or path to json file.
#' @return data frame with unique alignments and count.
#' @export
read_input = function(input) {
    # define data frame
    data = data.frame(Seq1 = character(), Seq2 = character())

    # if argument is a dir, look for all fasta files inside
    if(dir.exists(input)) {
        fasta_list = read_dir(input)
        if(length(fasta_list) < 1) {
            stop("No fasta files found in current directory.")
        }

        # read fasta files
        for(input in fasta_list) {
            # read fasta, concatenate to data frame
            data[nrow(data)+1, ] = read_fasta(input)
        }
    } else if(tools::file_ext(input) == "json") {
      # read json file
      data = read_json(input)
    } else {
      stop("Input format not supported.")
    }

    # check that all alignments without gaps are identical
    # count unique alignments and return data frame
    data = data %>% dplyr::count(Seq1, Seq2)
    data
}

#' Convert alignment to dot matrix.
#' This functions takes a pairwise alignment and creates a dot matrix
#' representing the alignment's path.
#'
#' @param a: string with DNA/AA sequence.
#' @param b: string with DNA/AA sequence.
#' @param n: number of ocurrences of alignment in list.
#' @return dot matrix.
#' @export
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

#' Plot a dot matrix.
#' This function plots a dot matrix as a raster image.
#'
#' @param  filepath: either a directory contains all alignments (DNA/AA), or a json file
#' @param  isDir: `TRUE` or `FALSE`. It is true by default that our filepath is a directory. otherwise, it is a json format.
#' @return A raster plot
#' @export
plot_dot_matrix = function(input, use_ggplot = FALSE){

    #read input
    aln = read_input(input)

    #initialize matrix
    mat <- 0

    #create dot matrix with all alignment paths
    for(i in seq_len(nrow(aln))) {
      mat <- mat + dot_matrix(aln$Seq1[i], aln$Seq2[i], aln$n[i])
    }

    #plot
    create_plot(mat, aln, use_ggplot)
}

#' Plot dot matrix.
#'
#' @param mat: dot matrix.
#' @param use_ggplot: logical. If TRUE, plot using ggplot2 package. Otherwise use base-R function.
#' @return A raster plot
#' @export
create_plot = function(mat, aln, use_ggplot) {
    #if(!use_ggplot) {
    plot(as.raster(1-mat/sum(aln$n)), interpolate = FALSE)
    #} # else {
      # use ggplot2 package
    # }

}

################################################################################
# Helper functions

# find all fasta files in `dir_path` and return a list of files
read_dir = function(dir_path) {
    files = list.files(path = dir_path, pattern = "\\.fasta$", full.names = TRUE)
    files = append(files, list.files(path = dir_path, pattern = "\\.fa$", full.names = TRUE))
    files
}

# read json file (coati-sample output format) and return as data frame
read_json = function(json_path) {
    data = jsonlite::fromJSON(json_path)
    data = data %>% tidyr::unpack(aln)
    data[c("Seq1", "Seq2")]
}

# read fasta file and return as vector
read_fasta = function(fasta_path) {
    data = data.frame(Seq1 = character(), Seq2 = character(), weight = double(), log_weight = double())
    fasta = seqinr::read.fasta(fasta_path, set.attributes = FALSE, as.string = TRUE, forceDNAtolower = FALSE)
    seq1 =  stringr::str_to_upper(fasta[[1]])
    seq2 =  stringr::str_to_upper(fasta[[2]])
    # follow data frame format: Seq1, Seq2, weight, log_weight
    c(seq1, seq2)
}

# check that, after removing gaps, all alignments are identical
check_seqs = function(data) {
    seq1 = gsub(pattern = "-", replacement = "", x = data$Seq1[1])
    seq2 = gsub(pattern = "-", replacement = "", x = data$Seq2[1])
    for(i in seq(from = 2, to = nrow(data))) {
        stopifnot(seq1 == gsub("-", "", data$Seq1[i]), seq2 == gsub("-", "", data$Seq2[i]))
    }
}
