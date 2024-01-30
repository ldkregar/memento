.datatable.aware=TRUE
#' Compute probabilites of possible assignments
#'
#' @param priors_per_tree list obtained from generate priors
#' @param m_reads_matrix matrix of methylated reads: ensure column and row names
#'  are present
#' @param u_reads_matrix matrix of unmethylated reads: ensure column and row
#' names are present
#' @param matrix_list matrix of phased reads: ensure column and row names are
#' present for each matrix and names of matrices are MT_M, MT_U, PT_M, PT_U,
#' U_U and U_M
#' @param ro numeric: error rate (only required if gamma and delta not prvided)
#' @param gamma numeric: proportion of unmethylated reads if inferred state is
#' methylated
#' @param delta numeric: proportion of methylated reads if inferred state is
#' unmethylated
#' @param verbose logical: whether to print messages
#' @param n_top_br integer: number of top branches to return
#' @param hash_name character: name of current combination (only for mle step)
#' @param current_comb only for mle step
#' @param site_specific logical: if using gammas and deltas, are they site or
#' colony specific?
#' @param use_ff_fully logical: whether to use ff fully or just for computing
#' probabilities
#'
#' @return data.table of probabilities and solutions
#' @export
#'
#' @examples
#' \dontrun{
#'   test <- compute_probabilities(priors_per_tree = result,
#'   m_reads_matrix = m_mtx, u_reads_matrix = u_mtx, ro = 0.95)
#' }
compute_probabilities <- function(priors_per_tree,
                                  m_reads_matrix=NULL,
                                  u_reads_matrix = NULL,
                                  matrix_list = NULL,
                                  ro=NA,
                                  gamma = NA,
                                  delta = NA,
                                  verbose = TRUE,
                                  n_top_br = 2,
                                  hash_name = NA,
                                  current_comb = NA,
                                  site_specific = TRUE,
                                  use_ff_fully = FALSE){
  chr_pos = branch = only_branch = starting_state = rank = n_assignments = Chrom = Pos = NULL
  if(length(priors_per_tree) == 1){
    only_one <- TRUE
    priors_per_tree <- priors_per_tree[[1]]
  }else{
    only_one <- FALSE
  }
  branch_lengths <- priors_per_tree$priors
  one_tree <- priors_per_tree$tree_for_priors
  genotype_vectors <- priors_per_tree$genotype_vectors
  branches <- priors_per_tree$branches
  n_permuted_b <- unique(lengths(priors_per_tree$affected_branches))

  # new gnt column names
  current_col_names <- lapply(genotype_vectors, function(g){
    return(colnames(g))
  })
  new_col_names <- vector("list", length(current_col_names))
  for(t in 1:length(current_col_names)){
    new_col_names[[t]] <- paste0(names(current_col_names)[t], ".", current_col_names[[t]])
  }
  vector_col_names <- unlist(new_col_names)
  genotype_matrix <- do.call("cbind", genotype_vectors)
  colnames(genotype_matrix) <- vector_col_names

  if(verbose){
    cat(sprintf("Now working on tree with %i tips: %s\n",
                length(one_tree$tip.label),
                paste(one_tree$tip.label, collapse = ", ")))
  }

  # prepare a matrix for each state, where tips is state S have 1 and 0 elsewhere
  prep_gnt_tf <- lapply(c(1:the$number_of_states), function(si){
    yy <- genotype_matrix == si
    mode(yy) <- "numeric"
    return(yy)
  })
  names(prep_gnt_tf) <- as.character(c(1:the$number_of_states))

  if(the$number_of_states < 4){
    nsites=dim(m_reads_matrix)[1]
    all_sites <- rownames(m_reads_matrix)
    methcount = m_reads_matrix
    wtcount = u_reads_matrix
    #if(mode(methcount) != "numeric"){
    mode(methcount) <- "numeric"
    #}
    #if(mode(wtcount) != "numeric"){
    mode(wtcount) <- "numeric"
    #}
  }else{
    nsites=dim(matrix_list[[1]])[1]
    all_sites=rownames(matrix_list[[1]])
    for(idx in 1:length(matrix_list)){
      mode(matrix_list[[idx]]) <- "numeric"
    }
  }

  only_prior_values_1 <- unlist(branch_lengths, recursive = FALSE)
  only_prior_values_2 <- sapply(only_prior_values_1, "[[", 1)

  if(all(is.na(gamma)  & is.na(delta))){
    stopifnot("Must provide gamma and delta or ro values." = !is.na(ro))
    a=log(ro)
    b=log(1-ro)
    c=log(0.5)
    d=log(0.5)
    e=log(1-ro)
    f=log(ro)
  }else{
    a=log(1-gamma)
    b=log(gamma)
    c=log((1+delta-gamma)/2)
    d=log((1+gamma-delta)/2)
    e=log(delta)
    f=log(1-delta)
  }

  prior1=matrix(rep(only_prior_values_2, nsites),
                nrow=nsites,
                ncol=length(only_prior_values_2), byrow = TRUE)

  if(the$number_of_states == 3){
    if(verbose){
      cat(sprintf("Dim meth: %i, %i\nDim unmeth: %i, %i\nDim genotype: %i, %i\n",
                  dim(methcount)[1], dim(methcount)[2],
                  dim(wtcount)[1], dim(wtcount)[2],
                  dim(prep_gnt_tf[["3"]])[1], dim(prep_gnt_tf[["3"]])[2]))
    }
    if(length(a) == 1 & length(e) == 1 & length(c) == 1){
      if(the$use_ff){
        # p = (ff::as.ff(a*methcount + b*wtcount))[] %*% (ff::as.ff(prep_gnt_tf[["3"]]))[] +
        #   (ff::as.ff(c*(methcount+wtcount)))[] %*% (ff::as.ff(prep_gnt_tf[["2"]]))[] +
        #   (ff::as.ff((b*methcount+a*wtcount)))[] %*% (ff::as.ff(prep_gnt_tf[["1"]]))[]
        #profvis({
          p=methods::as((a*methcount + b*wtcount), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["3"]], "HDF5Matrix") +
            methods::as((c*methcount + d*wtcount), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["2"]], "HDF5Matrix") +
            methods::as((e*methcount + f*wtcount), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["1"]], "HDF5Matrix")
          #})
      }else if(the$use_Rfast){
        p=Rfast::mat.mult((a*methcount + b*wtcount), prep_gnt_tf[["3"]]) +
          Rfast::mat.mult((c*methcount + d*wtcount), prep_gnt_tf[["2"]]) +
          Rfast::mat.mult((e*methcount + f*wtcount), prep_gnt_tf[["1"]])
      }else{
        p=eigenMapMatMult2((a*methcount + b*wtcount), prep_gnt_tf[["3"]]) +
          eigenMapMatMult2((c*methcount + d*wtcount), prep_gnt_tf[["2"]]) +
          eigenMapMatMult2((e*methcount + f*wtcount), prep_gnt_tf[["1"]])
      }
    }else{
      if(site_specific){
        if(length(a) > nrow(methcount)){
          a <- a[match(rownames(methcount), names(a))]
          b <- b[match(rownames(methcount), names(b))]
          c <- c[match(rownames(methcount), names(c))]
          d <- d[match(rownames(methcount), names(d))]
          e <- e[match(rownames(methcount), names(e))]
          f <- f[match(rownames(methcount), names(f))]
        }
        if(all(names(a) == names(e)) & all(names(a) == rownames(methcount)) & all(names(a) == rownames(wtcount))){
          # transpose matrix so you can do eachrow next
          methcount <- t(methcount)
          wtcount <- t(wtcount)
        }else{
          stop("Ensure names of gamma and delta match rownames of methcount and wtcount")
        }
      }else{
        if(any(is.null(names(a)) |
               names(a) != colnames(methcount) |
               names(a) != colnames(wtcount))){
          stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
        }else if(any(is.null(names(b)) |
                     names(b) != colnames(methcount) |
                     names(b) != colnames(wtcount))){
          stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
        }else{
          if(verbose){
            cat("Names of gamma & delta match column names of meth/wt count matrix.\n")
          }
        }
      }

      tmp1 <- Rfast::eachrow(methcount, a, "*")
      tmp2 <- Rfast::eachrow(wtcount, b, "*")
      tmp3 <- Rfast::eachrow(methcount, c, "*")
      tmp4 <- Rfast::eachrow(wtcount, d, "*")
      tmp5 <- Rfast::eachrow(methcount, e, "*")
      tmp6 <- Rfast::eachrow(wtcount, f, "*")

      if(site_specific){
        tmp1 <- t(tmp1)
        tmp2 <- t(tmp2)
        tmp3 <- t(tmp3)
        tmp4 <- t(tmp4)
        tmp5 <- t(tmp5)
        tmp6 <- t(tmp6)
      }

      if(the$use_ff){
        # new_p <- (ff::as.ff(tmp1+tmp2))[] %*% (ff::as.ff(prep_gnt_tf[["3"]]))[] +
        #   (ff::as.ff(tmp3))[] %*% (ff::as.ff(prep_gnt_tf[["2"]]))[] +
        #   (ff::as.ff(tmp4 + tmp5))[] %*% (ff::as.ff(prep_gnt_tf[["1"]]))[]

        new_p <- methods::as((tmp1+tmp2), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["3"]], "HDF5Matrix") +
          methods::as((tmp3+tmp4), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["2"]], "HDF5Matrix") +
          methods::as((tmp5 + tmp6), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["1"]], "HDF5Matrix")
      }else if(the$use_Rfast){
        new_p <- Rfast::mat.mult((tmp1+tmp2), prep_gnt_tf[["3"]]) +
          Rfast::mat.mult((tmp3+tmp4), prep_gnt_tf[["2"]]) +
          Rfast::mat.mult((tmp5 + tmp6), prep_gnt_tf[["1"]])
      }else{
        new_p <- eigenMapMatMult2((tmp1+tmp2), prep_gnt_tf[["3"]]) +
          eigenMapMatMult2((tmp3+tmp4), prep_gnt_tf[["2"]]) +
          eigenMapMatMult2((tmp5 + tmp6), prep_gnt_tf[["1"]])
      }
      p <- new_p
      }
    }else if(the$number_of_states == 4){
      if(the$use_ff){
        stop("low imprint code for 4 states not yet implemented, set the$use_ff to FALSE")
      }else if(the$use_Rfast){
        if(length(a) == 1 & length(b) == 1){
          state1 <- Rfast::mat.mult(matrix_list[["MT_M"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["PT_M"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["U_M"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["MT_U"]], prep_gnt_tf[["2"]]) +
            Rfast::mat.mult(matrix_list[["PT_M"]], prep_gnt_tf[["2"]]) +
            Rfast::mat.mult(matrix_list[["PT_U"]], prep_gnt_tf[["3"]]) +
            Rfast::mat.mult(matrix_list[["MT_M"]], prep_gnt_tf[["3"]]) +
            Rfast::mat.mult(matrix_list[["MT_U"]], prep_gnt_tf[["1"]]) +
            Rfast::mat.mult(matrix_list[["PT_U"]], prep_gnt_tf[["1"]]) +
            Rfast::mat.mult(matrix_list[["U_U"]], prep_gnt_tf[["1"]])

          state2 <- Rfast::mat.mult(matrix_list[["MT_U"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["PT_U"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["U_U"]], prep_gnt_tf[["4"]]) +
            Rfast::mat.mult(matrix_list[["MT_M"]], prep_gnt_tf[["2"]]) +
            Rfast::mat.mult(matrix_list[["PT_U"]], prep_gnt_tf[["2"]]) +
            Rfast::mat.mult(matrix_list[["PT_M"]], prep_gnt_tf[["3"]]) +
            Rfast::mat.mult(matrix_list[["MT_U"]], prep_gnt_tf[["3"]]) +
            Rfast::mat.mult(matrix_list[["PT_M"]], prep_gnt_tf[["1"]]) +
            Rfast::mat.mult(matrix_list[["MT_M"]], prep_gnt_tf[["1"]]) +
            Rfast::mat.mult(matrix_list[["U_M"]], prep_gnt_tf[["1"]])

          state3 <- Rfast::mat.mult(matrix_list[["U_M"]], prep_gnt_tf[["2"]]) +
            Rfast::mat.mult(matrix_list[["U_U"]], prep_gnt_tf[["2"]]) +

            Rfast::mat.mult(matrix_list[["U_M"]], prep_gnt_tf[["3"]]) +
            Rfast::mat.mult(matrix_list[["U_U"]], prep_gnt_tf[["3"]])

          p <- state1 * a + state2 * b + state3 * c
        }else{
          if(site_specific){
            if(length(a) > nrow(matrix_list[[1]])){
              a <- a[match(rownames(matrix_list[[1]]), names(a))]
              b <- b[match(rownames(matrix_list[[1]]), names(b))]
            }
            if(all(names(a) == names(b)) & all(names(a) == rownames(matrix_list[[1]]))){
              # transpose matrix so you can do eachrow next
              matrix_list <- lapply(matrix_list, function(onem){
                return(t(onem))
              })
            }else{
              stop("Ensure names of gamma and delta match rownames of matrix_list matrices")
            }
          }else{
            if(any(is.null(names(a)) |
                   names(a) != colnames(matrix_list[[1]]) |
                   names(a) != colnames(matrix_list[[1]]))){
              stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
            }else if(any(is.null(names(b)) |
                         names(b) != colnames(matrix_list[[1]]) |
                         names(b) != colnames(matrix_list[[1]]))){
              stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
            }else{
              if(verbose){
                cat("Names of gamma & delta match column names of meth/wt count matrix.\n")
              }
            }
          }

          c_vec <- rep(c, ncol(matrix_list[[1]]))
          # sweep(methcount, 2, a, "*")
          if(site_specific){
            state1 <- Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_U"]], a, "*")), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_M"]], a, "*")), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_U"]], a, "*")), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_M"]], a, "*")), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_U"]], a, "*")), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_U"]], a, "*")), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_U"]], a, "*")), prep_gnt_tf[["1"]])

            state2 <- Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_M"]], b, "*")), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_U"]], b, "*")), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_M"]], b, "*")), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_U"]], b, "*")), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["PT_M"]], b, "*")), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["MT_M"]], b, "*")), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_M"]], b, "*")), prep_gnt_tf[["1"]])

            state3 <- Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*")), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*")), prep_gnt_tf[["2"]]) +

              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*")), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(t(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*")), prep_gnt_tf[["3"]])
          }else{
            state1 <- Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_U"]], a, "*"), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_M"]], a, "*"), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_U"]], a, "*"), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_M"]], a, "*"), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_U"]], a, "*"), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_U"]], a, "*"), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_U"]], a, "*"), prep_gnt_tf[["1"]])

            state2 <- Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_M"]], b, "*"), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_U"]], b, "*"), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_M"]], b, "*"), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_U"]], b, "*"), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["PT_M"]], b, "*"), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["MT_M"]], b, "*"), prep_gnt_tf[["1"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_M"]], b, "*"), prep_gnt_tf[["1"]])

            state3 <- Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*"), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*"), prep_gnt_tf[["2"]]) +

              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*"), prep_gnt_tf[["3"]]) +
              Rfast::mat.mult(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*"), prep_gnt_tf[["3"]])
          }

          p <- state1 + state2 + state3
        }
      }else{
        if(length(a) == 1 & length(b) == 1){
          state1 <- eigenMapMatMult2(matrix_list[["MT_M"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["PT_M"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["U_M"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["MT_U"]], prep_gnt_tf[["2"]]) +
            eigenMapMatMult2(matrix_list[["PT_M"]], prep_gnt_tf[["2"]]) +
            eigenMapMatMult2(matrix_list[["PT_U"]], prep_gnt_tf[["3"]]) +
            eigenMapMatMult2(matrix_list[["MT_M"]], prep_gnt_tf[["3"]]) +
            eigenMapMatMult2(matrix_list[["MT_U"]], prep_gnt_tf[["1"]]) +
            eigenMapMatMult2(matrix_list[["PT_U"]], prep_gnt_tf[["1"]]) +
            eigenMapMatMult2(matrix_list[["U_U"]], prep_gnt_tf[["1"]])

          state2 <- eigenMapMatMult2(matrix_list[["MT_U"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["PT_U"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["U_U"]], prep_gnt_tf[["4"]]) +
            eigenMapMatMult2(matrix_list[["MT_M"]], prep_gnt_tf[["2"]]) +
            eigenMapMatMult2(matrix_list[["PT_U"]], prep_gnt_tf[["2"]]) +
            eigenMapMatMult2(matrix_list[["PT_M"]], prep_gnt_tf[["3"]]) +
            eigenMapMatMult2(matrix_list[["MT_U"]], prep_gnt_tf[["3"]]) +
            eigenMapMatMult2(matrix_list[["PT_M"]], prep_gnt_tf[["1"]]) +
            eigenMapMatMult2(matrix_list[["MT_M"]], prep_gnt_tf[["1"]]) +
            eigenMapMatMult2(matrix_list[["U_M"]], prep_gnt_tf[["1"]])

          state3 <- eigenMapMatMult2(matrix_list[["U_M"]], prep_gnt_tf[["2"]]) +
            eigenMapMatMult2(matrix_list[["U_U"]], prep_gnt_tf[["2"]]) +

            eigenMapMatMult2(matrix_list[["U_M"]], prep_gnt_tf[["3"]]) +
            eigenMapMatMult2(matrix_list[["U_U"]], prep_gnt_tf[["3"]])

          p <- state1 * a + state2 * b + state3 * c
        }else{
          if(site_specific){
            if(length(a) > nrow(matrix_list[[1]])){
              a <- a[match(rownames(matrix_list[[1]]), names(a))]
              b <- b[match(rownames(matrix_list[[1]]), names(b))]
            }
            if(all(names(a) == names(b)) & all(names(a) == rownames(matrix_list[[1]]))){
              # transpose matrix so you can do eachrow next
              matrix_list <- lapply(matrix_list, function(onem){
                return(t(onem))
              })
            }else{
              stop("Ensure names of gamma and delta match rownames of matrix_list matrices")
            }
          }else{
            if(any(is.null(names(a)) |
                   names(a) != colnames(matrix_list[[1]]) |
                   names(a) != colnames(matrix_list[[1]]))){
              stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
            }else if(any(is.null(names(b)) |
                         names(b) != colnames(matrix_list[[1]]) |
                         names(b) != colnames(matrix_list[[1]]))){
              stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
            }else{
              if(verbose){
                cat("Names of gamma & delta match column names of meth/wt count matrix.\n")
              }
            }
          }

          c_vec <- rep(c, ncol(matrix_list[[1]]))
          # sweep(methcount, 2, a, "*")
          if(site_specific){
            state1 <- eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_M"]], a, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_U"]], a, "*")), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_M"]], a, "*")), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_U"]], a, "*")), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_M"]], a, "*")), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_U"]], a, "*")), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_U"]], a, "*")), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_U"]], a, "*")), prep_gnt_tf[["1"]])

            state2 <- eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_U"]], b, "*")), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_M"]], b, "*")), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_U"]], b, "*")), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_M"]], b, "*")), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_U"]], b, "*")), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["PT_M"]], b, "*")), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["MT_M"]], b, "*")), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_M"]], b, "*")), prep_gnt_tf[["1"]])

            state3 <- eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*")), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*")), prep_gnt_tf[["2"]]) +

              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*")), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(t(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*")), prep_gnt_tf[["3"]])
          }else{
            state1 <- eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_M"]], a, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_U"]], a, "*"), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_M"]], a, "*"), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_U"]], a, "*"), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_M"]], a, "*"), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_U"]], a, "*"), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_U"]], a, "*"), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_U"]], a, "*"), prep_gnt_tf[["1"]])

            state2 <- eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_U"]], b, "*"), prep_gnt_tf[["4"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_M"]], b, "*"), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_U"]], b, "*"), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_M"]], b, "*"), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_U"]], b, "*"), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["PT_M"]], b, "*"), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["MT_M"]], b, "*"), prep_gnt_tf[["1"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_M"]], b, "*"), prep_gnt_tf[["1"]])

            state3 <- eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*"), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*"), prep_gnt_tf[["2"]]) +

              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_M"]], c_vec, "*"), prep_gnt_tf[["3"]]) +
              eigenMapMatMult2(Rfast::eachrow(matrix_list[["U_U"]], c_vec, "*"), prep_gnt_tf[["3"]])
          }

          p <- state1 + state2 + state3
        }
      }
      }else{
        if(the$use_ff){
          #stop("low imprint code for 2 states not yet implemented, set the$use_ff to FALSE")
          if(length(a) == 1 & length(b) == 1){
            p=methods::as((a*methcount + b*wtcount), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["2"]], "HDF5Matrix") +
              methods::as((b*methcount + a*wtcount), "HDF5Matrix") %*% methods::as(prep_gnt_tf[["1"]], "HDF5Matrix")
          }else{
            if(site_specific){
              if(length(a) > nrow(methcount)){
                a <- a[match(rownames(methcount), names(a))]
                b <- b[match(rownames(methcount), names(b))]
              }
              if(all(names(a) == names(b)) & all(names(a) == rownames(methcount)) & all(names(a) == rownames(wtcount))){
                # transpose matrix so you can do eachrow next
                methcount <- t(methcount)
                wtcount <- t(wtcount)
              }else{
                stop("Ensure names of gamma and delta match rownames of methcount and wtcount")
              }
            }else{
              if(any(is.null(names(a)) |
                     names(a) != colnames(methcount) |
                     names(a) != colnames(wtcount))){
                stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
              }else if(any(is.null(names(b)) |
                           names(b) != colnames(methcount) |
                           names(b) != colnames(wtcount))){
                stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
              }else{
                if(verbose){
                  message("Yay, names of gamma & delta match column names of meth/wt count matrix.\n")
                }
              }
            }
            if(site_specific){
              p_new = methods::as(t((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*"))), "HDF5Matrix") %*%
                methods::as(prep_gnt_tf[["2"]], "HDF5Matrix") +
                methods::as(t((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*"))), "HDF5Matrix") %*%
                methods::as(prep_gnt_tf[["1"]], "HDF5Matrix")
            }else{
              p_new = methods::as((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*")), "HDF5Matrix") %*%
                methods::as(prep_gnt_tf[["2"]], "HDF5Matrix") +
                methods::as((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*")), "HDF5Matrix") %*%
                methods::as(prep_gnt_tf[["1"]], "HDF5Matrix")
            }
            p <- p_new
          }
        }else if(the$use_Rfast){

          if(length(a) == 1 & length(b) == 1){
            p = Rfast::mat.mult((a*methcount + b*wtcount), prep_gnt_tf[["2"]]) +
              Rfast::mat.mult((b*methcount + a*wtcount), prep_gnt_tf[["1"]])
          }else{
            if(site_specific){
              if(length(a) > nrow(methcount)){
                a <- a[match(rownames(methcount), names(a))]
                b <- b[match(rownames(methcount), names(b))]
              }
              if(all(names(a) == names(b)) & all(names(a) == rownames(methcount)) & all(names(a) == rownames(wtcount))){
                # transpose matrix so you can do eachrow next
                methcount <- t(methcount)
                wtcount <- t(wtcount)
              }else{
                stop("Ensure names of gamma and delta match rownames of methcount and wtcount")
              }
            }else{
              if(any(is.null(names(a)) |
                     names(a) != colnames(methcount) |
                     names(a) != colnames(wtcount))){
                stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
              }else if(any(is.null(names(b)) |
                           names(b) != colnames(methcount) |
                           names(b) != colnames(wtcount))){
                stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
              }else{
                if(verbose){
                  message("Yay, names of gamma & delta match column names of meth/wt count matrix.\n")
                }
              }
            }
            if(site_specific){
              p_new <- Rfast::mat.mult(t((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*"))), prep_gnt_tf[["2"]]) +
                Rfast::mat.mult(t((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*"))), prep_gnt_tf[["1"]])
            }else{
              p_new <- Rfast::mat.mult((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*")), prep_gnt_tf[["2"]]) +
                Rfast::mat.mult((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*")), prep_gnt_tf[["1"]])
            }
            p <- p_new
          }

        }else{
          if(length(a) == 1 & length(b) == 1){
            p = eigenMapMatMult2((a*methcount + b*wtcount), prep_gnt_tf[["2"]]) +
              eigenMapMatMult2((b*methcount + a*wtcount), prep_gnt_tf[["1"]])
          }else{
            if(site_specific){
              if(length(a) > nrow(methcount)){
                a <- a[match(rownames(methcount), names(a))]
                b <- b[match(rownames(methcount), names(b))]
              }
              if(all(names(a) == names(b)) & all(names(a) == rownames(methcount)) & all(names(a) == rownames(wtcount))){
                # transpose matrix so you can do eachrow next
                methcount <- t(methcount)
                wtcount <- t(wtcount)
              }else{
                stop("Ensure names of gamma and delta match rownames of methcount and wtcount")
              }
            }else{
              if(any(is.null(names(a)) |
                     names(a) != colnames(methcount) |
                     names(a) != colnames(wtcount))){
                stop("Names of gamma vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of meth/wt count matrix")
              }else if(any(is.null(names(b)) |
                           names(b) != colnames(methcount) |
                           names(b) != colnames(wtcount))){
                stop("Names of delta vector are missing or do not match names of meth or unmeth count matrices. Ensure gamma and delta values are ordered according to column names of the meth/wt count matrix")
              }else{
                if(verbose){
                  message("Yay, names of gamma & delta match column names of meth/wt count matrix.\n")
                }
              }
            }
            if(site_specific){
              p_new <- eigenMapMatMult2(t((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*"))), prep_gnt_tf[["2"]]) +
                eigenMapMatMult2(t((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*"))), prep_gnt_tf[["1"]])
            }else{
              p_new <- eigenMapMatMult2((Rfast::eachrow(methcount, a, "*") + Rfast::eachrow(wtcount, b, "*")), prep_gnt_tf[["2"]]) +
                eigenMapMatMult2((Rfast::eachrow(methcount, b, "*")+Rfast::eachrow(wtcount, a, "*")), prep_gnt_tf[["1"]])
            }
            p <- p_new
            }
          }
      }
  if(is.null(colnames(p))){
    colnames(p) <- colnames(prep_gnt_tf[["1"]])
  }

  if(is.null(rownames(p))){
    rownames(p) <- all_sites
  }
  if(the$use_ff){
    raw_val <- p + methods::as(prior1, "HDF5Matrix")
  }else{
    raw_val <- p + prior1
  }

  if(!is.na(hash_name)){
    the$global_hash_list[[hash_name]][["raw_vals"]][[as.character(current_comb)]] <- list(raw_val)
  }

  if(use_ff_fully & the$use_ff){
    raw_rowsum <- DelayedArray::rowSums(raw_val)

    max_raw <- apply(raw_val, 1, max)

    max_subtracted <- sweep(as(raw_val, "Matrix"), 1, max_raw)

    max_prob <- exp(max_subtracted) / apply(exp(max_subtracted), 1, sum)

    how_many <- min(n_top_br, ncol(max_prob))
    comb_name <- colnames(max_prob)

    max_prob_both <- DelayedArray::blockApply(x = max_prob,
                                   FUN = function(row1){
                                     y <- order(row1, decreasing = TRUE)[1:how_many]
                                     indi <- comb_name[y]
                                     valu <- row1[y]
                                     return(list("vals" = valu,
                                                 "inds" = indi))
                                   },
                                   grid = S4Arrays::RegularArrayGrid(refdim = dim(max_prob),
                                                           spacings = c(1L, ncol(max_prob))))
    max_prob_values <- lapply(max_prob_both, "[[", 1)
    max_prob_indices <- lapply(max_prob_both, "[[", 2)
  }else{
    if(the$use_ff){
      raw_val <- as.matrix(raw_val)
    }
    raw_rowsum <- rowSums(raw_val)

    max_raw <- apply(raw_val, 1, max)

    max_subtracted <- sweep(raw_val, 1, max_raw)
    max_prob <- exp(max_subtracted) / apply(exp(max_subtracted), 1, sum)

    how_many <- min(n_top_br, ncol(max_prob))
    comb_name <- colnames(max_prob)

    max_prob_both <- apply(max_prob, 1, function(x){
      sort(x, decreasing = TRUE)[1:how_many]
    }, simplify = FALSE)

    max_prob_values <- max_prob_both

    max_prob_indices <- lapply(max_prob_both, function(oneset){
      return(names(oneset))
    })
  }

  dt_prob <- data.table::data.table(do.call("rbind", max_prob_values))
  colnames(dt_prob) <- as.character(c(1:n_top_br))
  dt_prob[, chr_pos := all_sites]

  dt_solution <- data.table::data.table(do.call("rbind", max_prob_indices))
  colnames(dt_solution) <- as.character(c(1:n_top_br))
  dt_solution[, chr_pos := all_sites]

  dt_prob_molten <- data.table::melt(dt_prob, id.vars = c("chr_pos"),
                                     value.name = "probability",
                                     variable.name = "rank")

  dt_sol_molten <- data.table::melt(dt_solution, id.vars = c("chr_pos"),
                                       value.name = "branch",
                                       variable.name = "rank")

  dt_sol_molten[, c("only_branch", "starting_state", "branch_state") := tstrsplit(branch, "\\.|:")]
  dt_sol_molten[, starting_state := as.numeric(starting_state)]
  dt_new <- data.table::merge.data.table(dt_prob_molten, dt_sol_molten, by = c("chr_pos", "rank"))
  dt_new[, rank := as.numeric(rank)]

  if(!("Chrom" %in% colnames(dt_new))){
    dt_new[, c("Chrom", "Pos") := tstrsplit(chr_pos, "_")]
    dt_new[, Pos := as.numeric(Pos)]
  }

  if(only_one){
    list_rate <- lapply(dt_new$branch,
                      get_rate,
                      priors_w_idx = branch_lengths,
                      full_tree = one_tree)
    names(list_rate) <- dt_new$chr_pos

    rates_df <- do.call("rbind", list_rate)
    rates_df2 <- data.table::data.table("chr_pos" = rownames(rates_df),
                                        "rank" = dt_new$rank,
                                        rates_df)
    data.table::setkey(dt_new, chr_pos, rank)
    dt_new2 <- data.table::merge.data.table(dt_new, rates_df2, by = c("chr_pos", "rank"))
    dt_new <- dt_new2
  }
  dt_new[, n_assignments := sapply(only_branch, function(x){
    y <- strsplit(x, "_")[[1]]
    if(all(y == "n")){
      return(0)
    }else{
      return(length(y))
    }
  })]
  return(dt_new)
}

#' Get MLE solution for branches at reads
#'
#' @param branches character of candidate branches
#' @param reads named list with items mm and um that holds methylated and
#' unmethlyated reads per site for all tips
#' @param forced_starting_states starting state (if NULL - will try all
#' starting states)
#' @param return_edge_changes logical
#' @param priors_per_branch list of priors per branch
#' @param full_tree ape phylo object of full tree
#' @param ro conversion rate
#' @param gamma error rate: prop. of unmethylated reads when state is meth
#' @param delta error rate: prop. of methylated reads when state is unmeth
#' @param verbose_ logical: if verbose
#' @param n_top_branches integer number of top branches to report
#' @param delete_intermediate_matrices logical: if intermediate matrices should
#' be deleted
#' @param max_loops_ max number of branches to drop
#' @param site_specific whether gamma and delta are site specific (alternative
#' is colony specific)
#' @param loop_through_all logical - whether to drop all candidate branches until n_to_permute = 0
#'
#' @return data.table of most likely solution
#' @export
#'
#' @examples
#' \dontrun{
#'     tmp_df <- get_mle(branches = assignment_space[[counter]],
#'     reads = assignment_reads[[counter]],
#'     priors_per_branch = priors_per_branch, full_tree = full_tree,
#'     ro = ro, gamma = gamma, delta = delta, verbose_ = verbose,
#'     n_top_branches = n_top_br,
#'     delete_intermediate_matrices = delete_matrices_from_hash,
#'     max_loops_ = max_loops,
#'     site_specific = site_specific,
#'     forced_starting_states = fss[[counter]])
#' }
get_mle <- function(branches,
                    reads,
                    forced_starting_states = NULL,
                    return_edge_changes = TRUE,
                    priors_per_branch,
                    full_tree,
                    ro = NA,
                    gamma = NA,
                    delta = NA,
                    verbose_ = FALSE,
                    n_top_branches = 1,
                    delete_intermediate_matrices = TRUE,
                    max_loops_ = NULL,
                    site_specific = TRUE,
                    loop_through_all = FALSE){
  # branches <- assignment_space[[4]]
  # sites <- assignment_sites[[4]]
  chr_pos = branch = starting_state  = probability = only_branch = n_assignments = branch_state = NULL
  if(branches != "n"){
    branches_split <- as.numeric(strsplit(branches, "_|,")[[1]])
    branches_split <- sort(branches_split)
  }else{
    branches_split <- branches
  }

  # while probability is highest for smallest set of branches
  # try combination with one less branch to see if fit improves
  p_last <- TRUE
  if(branches == "n"){
    n_to_permute <- 0
  }else{
    n_to_permute <- length(branches_split)
  }
  while_counter <- 0
  if(is.null(max_loops_)){
    while_counter_limit <- Inf
  }else{
    while_counter_limit <- max_loops_
    if(max_loops_ < 1){
      stop("Max loops must be at least one")
    }
  }
  n_branches <- length(branches_split)
  # while(p_last & while_counter < while_counter_limit){
  while(p_last){
    while_counter <- while_counter + 1
    if(verbose_) cat(sprintf("n branches to permute: %i\t", n_to_permute))

    if(branches != "n" & n_to_permute > 0){
      # generate all possible permutations of up to n_branches
      br_permut <- generate_permutations(n = n_branches,
                                         m = n_to_permute,
                                         only_one = TRUE)

      colnames(br_permut) <- as.character(branches_split)

      affected_branches_full <- apply(br_permut, 1, function(x){
        if(all(x == 0)){
          return(c())
        }else{
          return(as.numeric(names(which(x == 1))))
        }
      }, simplify = FALSE)

      affected_branches_strings <- sapply(affected_branches_full, function(x){
        if(length(x) > 0){
          return(paste0(x, collapse = "_"))
        }else{
          return("n")
        }
      })

      names(affected_branches_full) <- affected_branches_strings

    }else{
      affected_branches_full <- list("n" = c())
    }
    # compute priors
    if(!is.null(forced_starting_states)){
      ss_of_comb <- unique(forced_starting_states)
      if(n_to_permute == 0){
        ss_of_comb <- the$all_possible_states
      }
    }else{
      ss_of_comb <- NULL
    }
    br_priors_w_idx <- tryCatch(lapply(affected_branches_full,
                                       generate_priors_hemimeth,
                                       tree_ = full_tree,
                                       priors_p_branch = priors_per_branch,
                                       return_idx = TRUE,
                                       force_ss = ss_of_comb),
                                warning = function(w) {
                                  cat(paste0("warning!! : "))
                                })

    full_tree_aslist <- lapply(1:length(br_priors_w_idx), function(x) return(full_tree))

    br_full_affected_tips <- mapply(map_branch_changes_to_affected_tips_hemi,
                                    pr = br_priors_w_idx,
                                    tree_ = full_tree_aslist,
                                    SIMPLIFY = FALSE)

    br_genotype_vectors <- lapply(br_full_affected_tips,
                                  generate_genotypes,
                                  original_tree = full_tree,
                                  root_tree = NULL)

    tree_for_priors <- full_tree

    full_result <- list("priors" = br_priors_w_idx,
                        "affected_tips" = br_full_affected_tips,
                        "affected_branches" = affected_branches_full,
                        "tree_for_priors" = tree_for_priors,
                        "genotype_vectors" = br_genotype_vectors,
                        "branches" = branches)

    # generate probabilities for the tips
    if(the$number_of_states < 4){
      br_probs_new <- compute_probabilities(priors_per_tree = full_result,
                                            m_reads_matrix=reads$mm,
                                            u_reads_matrix = reads$um,
                                            ro = ro,
                                            verbose = FALSE,
                                            gamma = gamma,
                                            delta = delta,
                                            hash_name = branches,
                                            current_comb = n_to_permute,
                                            site_specific = site_specific)
    }else{
      br_probs_new <- compute_probabilities(priors_per_tree = full_result,
                                            matrix_list = reads,
                                            ro = ro,
                                            verbose = FALSE,
                                            gamma = gamma,
                                            delta = delta,
                                            hash_name = branches,
                                            current_comb = n_to_permute,
                                            site_specific = site_specific)
    }

    states_rates <- lapply(br_probs_new$branch,
                           get_rate,
                           priors_w_idx = br_priors_w_idx,
                           full_tree = full_tree)
    names(states_rates) <- br_probs_new$chr_pos
    the$global_hash_list[[branches]][["states"]][[as.character(n_to_permute)]] <- states_rates

    if(branches != "n"){ # combination has assignments
      # combine raw values, find max
      # elements <- global_hash_list[[branches]][["raw_vals"]]
      # eligible_el <- elements[elements != "" & elements != "states"]
      eligible_el <- the$global_hash_list[[branches]][["raw_vals"]]
      if(length(eligible_el) == 1 & while_counter < while_counter_limit){ # only one done, must try a smaller set

      }else{
        all_raw <- eligible_el #global_hash_list[[branches]][["raw_vals"]][eligible_el]
        all_raw_matrices <- lapply(all_raw, "[[", 1)
        if(the$use_ff){
          all_raw_bound <- as(do.call(DelayedArray::cbind, all_raw_matrices), "Matrix")
        }else{
          all_raw_bound <- do.call("cbind", all_raw_matrices)
        }
        max_per_row <- apply(all_raw_bound, 1, max)
        max_per_row_ind <- apply(all_raw_bound, 1, which.max)

        tmp_per_part <- exp(sweep(all_raw_bound, 1, max_per_row)) / apply(exp(sweep(all_raw_bound, 1, max_per_row)), 1, sum)

        names_pp <- colnames(tmp_per_part)
        this_many <- min(n_top_branches, ncol(tmp_per_part))
        # sort
        ranked_max_per_row_prob <- apply(tmp_per_part, 1, function(x){
          sort(x, decreasing = TRUE)[1:this_many]
        }, simplify = FALSE)

        ranked_max_per_row_ind_prob <- lapply(ranked_max_per_row_prob, function(item){
          return(names(item))
        })
        # ranked_max_per_row_ind_prob <- apply(tmp_per_part, 1, function(x){
        #   y <- order(x, decreasing = TRUE)[1:this_many]
        #   z <- names_pp[y]
        #   return(z)
        # }, simplify = FALSE)

        max_per_row_prob <- apply(tmp_per_part, 1, max)
        max_per_row_ind_prob <- apply(tmp_per_part, 1, which.max)
        #}
        # are all max in smallest set?
        indices_per_row <- colnames(all_raw_bound)[max_per_row_ind_prob]
        n_as_per_win_idx <- sapply(indices_per_row, function(full_name){
          ff <- strsplit(full_name, "\\.")[[1]][1]
          gg <- strsplit(ff, "_")[[1]]
          return(length(gg))
        })

        # if you're not in smallest set and you are happy with local optimal solution (you dont explore all -> exit while loop)
        if(all(n_as_per_win_idx > n_to_permute) & loop_through_all == FALSE){
          p_last <- FALSE
        }

        # loop_through_all
        # if(all(n_as_per_win_idx > n_to_permute) | while_counter == while_counter_limit | n_to_permute == 0){ #if all best set are not in the smallest set, exit while loop
        #   p_last <- FALSE
          # also update br_max_lik matrix
          # dt_prob <- data.table::data.table(do.call("rbind", ranked_max_per_row_prob))
          # colnames(dt_prob) <- as.character(c(1:n_top_branches))
          # dt_prob[, chr_pos := names(ranked_max_per_row_prob)]
          #
          # dt_solution <- data.table::data.table(do.call("rbind", ranked_max_per_row_ind_prob))
          # colnames(dt_solution) <- as.character(c(1:n_top_branches))
          # dt_solution[, chr_pos := names(ranked_max_per_row_ind_prob)]
          #
          # dt_prob_molten <- data.table::melt(dt_prob, id.vars = c("chr_pos"),
          #                                    value.name = "probability",
          #                                    variable.name = "rank")
          #
          # dt_sol_molten <- data.table::melt(dt_solution, id.vars = c("chr_pos"),
          #                                   value.name = "branch",
          #                                   variable.name = "rank")
          #
          # dt_sol_molten[, c("only_branch", "starting_state", "branch_state") := tstrsplit(branch, "\\.|:")]
          # dt_sol_molten[, starting_state := as.numeric(starting_state)]
          # dt_new <- data.table::merge.data.table(dt_prob_molten, dt_sol_molten, by = c("chr_pos", "rank"))
          # dt_new[, rank := as.numeric(rank)]
          # new_br_max_lik <- dt_new
        # }

      }

      # if(while_counter == while_counter_limit | n_to_permute == 0){ #if all best set are not in the smallest set, exit while loop
      #   p_last <- FALSE
      # }
      # if max is smallest set, compute for smaller perm
      n_to_permute <- n_to_permute - 1
      # if(n_to_permute < 0){
      if(n_to_permute < 0 | while_counter >= while_counter_limit){
        p_last <- FALSE
      }

    }else{
      # assigned to root, can exit while
      p_last <- FALSE
      # new_br_max_lik <- br_probs_new[, list(chr_pos, rank, probability, branch, only_branch, starting_state, branch_state)]
    }
  }

  if(branches != "n"){
    dt_prob <- data.table::data.table(do.call("rbind", ranked_max_per_row_prob))
    colnames(dt_prob) <- as.character(c(1:n_top_branches))
    dt_prob[, chr_pos := names(ranked_max_per_row_prob)]

    dt_solution <- data.table::data.table(do.call("rbind", ranked_max_per_row_ind_prob))
    colnames(dt_solution) <- as.character(c(1:n_top_branches))
    dt_solution[, chr_pos := names(ranked_max_per_row_ind_prob)]

    dt_prob_molten <- data.table::melt(dt_prob, id.vars = c("chr_pos"),
                                       value.name = "probability",
                                       variable.name = "rank")

    dt_sol_molten <- data.table::melt(dt_solution, id.vars = c("chr_pos"),
                                      value.name = "branch",
                                      variable.name = "rank")

    dt_sol_molten[, c("only_branch", "starting_state", "branch_state") := tstrsplit(branch, "\\.|:")]
    dt_sol_molten[, starting_state := as.numeric(starting_state)]
    dt_new <- data.table::merge.data.table(dt_prob_molten, dt_sol_molten, by = c("chr_pos", "rank"))
    dt_new[, rank := as.numeric(rank)]
    new_br_max_lik <- dt_new
  }else{
    new_br_max_lik <- br_probs_new[, list(chr_pos, rank, probability, branch, only_branch, starting_state, branch_state)]
  }

  new_br_max_lik[,n_assignments := sapply(only_branch, function(x){
    y <- strsplit(x, "_")[[1]]
    if(all(y == "n")){
      return(0)
    }else{
      return(length(y))
    }
  })]
  rates_for_win <- the$global_hash_list[[branches]][["states"]]
  rates_for_win2 <- vector("list", nrow(new_br_max_lik))
  for(q in 1:nrow(new_br_max_lik)){
    rates_for_win2[[q]] <- rates_for_win[[as.character(new_br_max_lik[q,n_assignments])]][[new_br_max_lik[q,chr_pos]]]
  }
  names(rates_for_win2) <- new_br_max_lik$chr_pos
  # remove from hash matrix
  the$global_hash_list[[branches]][["num_permuts"]] <- sapply(the$global_hash_list[[branches]][["raw_vals"]], function(x){
    ncol(x[[1]])
  })
  if(delete_intermediate_matrices){
    the$global_hash_list[[branches]][["raw_vals"]] <- NULL
  }
  rates_df <- do.call("rbind", rates_for_win2)
  rates_df2 <- data.table::data.table("chr_pos" = rownames(rates_df),
                          "rank" = new_br_max_lik$rank,
                          rates_df)
  data.table::setkey(new_br_max_lik, chr_pos, rank)
  lik_rates <- data.table::merge.data.table(new_br_max_lik, rates_df2, by = c("chr_pos", "rank"))
  return(lik_rates)
}


#' Find the most likely set of changed branches
#'
#' @param dt_pre data.table of assignments
#' @param full_tree "phylo" object of full tree
#' @param ro numeric: error rate (only required if gamma and delta not prvided)
#' @param gamma numeric value for proportion of unmethylated reads if inferred state is methylated
#' @param delta numeric value for proportion of methylated reads if inferred state is unmethylated
#' @param epi_fwd numeric value for forward rate (gain of methylation)
#' @param epi_bck numeric value for backward rate (loss of methylation)
#' @param Q_mtx Q matrix
#' @param fraction_estimate numeric value for how many sites to estimate MLE for
#' @param n_top_br integer: number of top branches to return
#' @param dt_anno data.table of site annotations, where chr_pos is key
#' @param file_to_write where to write solutions
#' @param m_mtx matrix of methylated reads
#' @param u_mtx matrix of unmethylated reads
#' @param mtx_list list of phased methylated and unmethylated reads
#' @param verbose logical whether to print messages
#' @param delete_matrices_from_hash logical whether to remove raw probability matrices from hash table
#' @param priors_per_branch custom priors per branch
#' @param max_loops integer maximum number of branches to drop when permuting candidate branches
#' @param site_specific logical: if providing gammas & deltas, are they colony or site specific
#' @param force_starting_state logical: should only certain start states be explored?
#' @param loop_through_all logical - whether to drop all candidate branches until n_to_permute = 0
#'
#' @return data.table of final assignments
#' @export
#'
#' @examples
#' \dontrun{
#' mle_branch_dt <- find_most_likely_set(dt_pre =  converted_dt,
#' full_tree = full_tree, ro = 0.95, n_top_br = 2, m_mtx = meth_mtx,
#' u_mtx = unmeth_mtx, fraction_estimate = 1, verbose = FALSE)
#' }
find_most_likely_set <- function(dt_pre,
                                 full_tree,
                                 ro = NA,
                                 gamma = NA,
                                 delta = NA,
                                 epi_fwd = NULL,
                                 epi_bck = NULL,
                                 Q_mtx = NULL,
                                 fraction_estimate = 1,
                                 n_top_br = 1,
                                 dt_anno = NULL,
                                 file_to_write = NULL,
                                 m_mtx = NULL,
                                 u_mtx = NULL,
                                 mtx_list = NULL,
                                 verbose = TRUE,
                                 delete_matrices_from_hash = TRUE,
                                 priors_per_branch = NULL,
                                 max_loops = NULL,
                                 site_specific = TRUE,
                                 force_starting_state = FALSE,
                                 loop_through_all = FALSE){
  chr_pos = combined_branch = Chrom = Pos = NULL
  if(all(is.na(gamma)  & is.na(delta))){
    stopifnot("Must provide gamma and delta or ro values." = !is.na(ro))
  }

  assignment_space <- unique(dt_pre[["combined_branch"]])
  assignment_space <- assignment_space[!is.na(assignment_space) & assignment_space != ""]

  assignment_sites <- lapply(assignment_space, function(b){
    positions <- dt_pre[combined_branch == b][["chr_pos"]]
    return(positions)
  })

  the$global_hash_list <- hash::hash(assignment_space, vector("list", length = length(assignment_space)))
  if(is.null(priors_per_branch)){
    totsnvs=sum(full_tree$edge.length)
    if(the$tree_branch_is_snv){
      tottime=totsnvs/the$snv_rate_per_year
    }else{
      tottime = totsnvs
    }
    weighted_branch_lenghts <- tottime* full_tree$edge.length/sum(full_tree$edge.length)

    if(is.null(epi_fwd) & is.null(epi_bck)){
      epi_mut_rate_fwd <- the$epi_fwd*the$human.snv.rate.per.base #human.epi.rate.per.cpg
      epi_mut_rate_bck <- the$epi_bck*the$human.snv.rate.per.base #human.epi.rate.per.cpg
    }else{
      cat("Using custom forward and backward transition rates.\n")
      epi_mut_rate_fwd <- epi_fwd*the$human.snv.rate.per.base #human.epi.rate.per.cpg
      epi_mut_rate_bck <- epi_bck*the$human.snv.rate.per.base #human.epi.rate.per.cpg
    }

    if(is.null(Q_mtx)){
      Q_mtx <- the$Q
    }else{
      cat("Using custom Q matrix")
    }

    if(the$useQ){
      p_list <- lapply(weighted_branch_lenghts,
                       getP,
                       Q = Q_mtx)
    }else{
      if(the$number_of_states == 2){
        p_list <- lapply(weighted_branch_lenghts,
                         getP2,
                         epi_mut_rate_bck = epi_mut_rate_bck,
                         epi_mut_rate_fwd = epi_mut_rate_fwd)
      }else if(the$number_of_states == 3){
        p_list <- lapply(weighted_branch_lenghts,
                         getP3,
                         epi_mut_rate_bck = epi_mut_rate_bck,
                         epi_mut_rate_fwd = epi_mut_rate_fwd)
      }else if(the$number_of_states == 4){
        p_list <- lapply(weighted_branch_lenghts,
                         getP4,
                         epi_mut_rate_bck = epi_mut_rate_bck,
                         epi_mut_rate_fwd = epi_mut_rate_fwd)
      }
    }

    p_vectors <- lapply(p_list, function(x){
      y <- as.vector(t(x))
      return(y)
    })
    priors_per_branch <- log(do.call("rbind", p_vectors))
  }

  if(fraction_estimate < 1){
    warning(sprintf("You are only computing MLE on the first %i%% of sites, sorted from easiest.", round(fraction_estimate*100, digits = 0)))
  }

  assignment_space_split <- lapply(assignment_space, function(x){
    strsplit(x, ",|_")[[1]]
  })
  reorder <- order(lengths(assignment_space_split))
  assignment_space <- assignment_space[reorder]
  assignment_sites <- assignment_sites[reorder]

  assignment_reads <- lapply(assignment_sites, function(site){
    if(the$number_of_states < 4){
      sites_mreads_vector <- convert_to_matrix(matrix = m_mtx, chr_pos = site)
      sites_ureads_vector <- convert_to_matrix(matrix = u_mtx, chr_pos = site)
      return(list("mm" = sites_mreads_vector,
                  "um" = sites_ureads_vector))
    }else{
      matrix_list_sub <- lapply(mtx_list, convert_to_matrix, chr_pos = site)
      return(matrix_list_sub)
    }
  })

  do_now <- ceiling(fraction_estimate*length(assignment_space))

  tmp_list <- vector("list", length = length(assignment_space))
  if("rank" %in% colnames(dt_pre)){
    dt_pre$rank <- NULL
  }
  for(counter in 1:do_now){
    if(verbose) cat(sprintf("%i of %i: computing MLE for branch %s\n", counter, do_now, assignment_space[[counter]]))
    if(force_starting_state){
      fss <- lapply(assignment_sites, function(set_of_sites){
        ss_root <- unique(dt_pre[match(set_of_sites, chr_pos),][["ss_root_tree"]])
        return(ss_root)
      })
      n_idx <- match("n", assignment_space)
      if(!is.na(n_idx)){
        fss[[n_idx]] <- the$all_possible_states
      }
    }else{
      fss <- NULL
    }
    tmp_df <- get_mle(branches = assignment_space[[counter]],
                      reads = assignment_reads[[counter]],
                      priors_per_branch = priors_per_branch,
                      full_tree = full_tree,
                      ro = ro,
                      gamma = gamma,
                      delta = delta,
                      verbose_ = verbose,
                      n_top_branches = n_top_br,
                      delete_intermediate_matrices = delete_matrices_from_hash,
                      max_loops_ = max_loops,
                      site_specific = site_specific,
                      forced_starting_states = fss[[counter]],
                      loop_through_all = loop_through_all)
    data.table::setnames(tmp_df, "n_assignments", "n_assignments_mle")
    dt_post <- data.table::merge.data.table(dt_pre, tmp_df, all.y = TRUE, by = "chr_pos")
    if(!is.null(dt_anno)){
      dt_post <- data.table::merge.data.table(dt_post, dt_anno, all.x = TRUE, by = "chr_pos")
    }
    if(!is.null(file_to_write)){
      if(file.exists(file_to_write)){
        data.table::fwrite(dt_post, file_to_write, append = TRUE)
      }else{
        data.table::fwrite(dt_post, file_to_write)
      }
    }

    tmp_list[[counter]] <- dt_post
  }

  branch_mle <- data.table::data.table(plyr::rbind.fill(tmp_list))

  data.table::setkey(branch_mle, chr_pos, rank)

  updated_branch_dt <- branch_mle

  if(!("Chrom" %in% colnames(updated_branch_dt))){
    updated_branch_dt[, list("Chrom", "Pos") := data.table::tstrsplit(chr_pos, "_")]
    updated_branch_dt[, Pos := as.numeric(Pos)]
  }

  return(updated_branch_dt)
}



