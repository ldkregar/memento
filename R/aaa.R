the <- new.env(parent = emptyenv())
the$number_of_states <- 3
the$human.snv.rate.per.base=20/3e9
the$epi_fwd = 1e5
the$epi_bck = 1e5
the$all_possible_states <- c(1:the$number_of_states)
the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                   the$number_of_states,
                                   the$number_of_states,
                                   byrow = TRUE)
colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

rebase <- function(b, base = 10){
  base <- as.integer(base)

  znj <- lapply(b, function(x){
    n   <- ceiling(log(x, base))
    vec <- numeric()
    val <- x

    while(n >= 0){
      rem <- val %/% base^n
      val <- val - rem * base^n
      vec <- c(vec, rem)
      n <- n - 1
    }

    while(vec[1] == 0 & length(vec) > 1) vec <- vec[-1]
    return(vec)
  })
  return(unlist(znj))
}

add_padding <- function(v, new_length){
  while(length(v) < new_length){
    v <- c(0, v)
  }
  return(v)
}

get_change_combinations_4s <- function(n_changes_){
  possible_perm <- 1:((get("number_of_states", envir = the) -1)**n_changes_)-1
  change_combination <- sapply(possible_perm, function(x){ rebase(x, (get("number_of_states", envir = the) -1))} )
  last_item <- change_combination[[length(change_combination)]]
  change_combination2 <- lapply(change_combination, add_padding, new_length = length(last_item))
  change_combination <- do.call("rbind", change_combination2)
  names_of_changes <- apply(change_combination, 1, function(x){
    y <- as.character(x)
    y[y == "0"] <- "A"
    y[y == "1"] <- "B"
    y[y == "2"] <- "C"
    y[y == "3"] <- "D"
    return(paste0(y, collapse = ""))
  })
  change_combination <- change_combination + 1 # this way you get indices
  return(list("change_combinations" = change_combination,
              "names_of_changes" = names_of_changes))
}

the$change_combs_list <- lapply(c(0:10), function(nc){
  cc <- get_change_combinations_4s(nc)
  return(cc)
})
names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))

getRateDf <- function(numStates){
  do.call("rbind",
          lapply(0:(numStates-1),
                 function(s) data.frame(s.parent=s,s.child=setdiff(0:(numStates-1),s))
          ))
}

the$df_rate <- getRateDf(the$number_of_states)
the$df_rate$rate <- rep(0.001, nrow(the$df_rate))

getQMatrix <- function(df){
  N=max(df$s.parent)+1
  out=matrix(NA,ncol=N,nrow=N)
  out[N*df$s.child+(df$s.parent+1)]=df$rate
  diag(out)=-rowSums(out,na.rm = TRUE)
  out
}

the$Q <- getQMatrix(the$df_rate)
the$useQ <- FALSE

#' Report my default settings
#'
#' @export
memento_settings <- function() {
  return(list("number_of_states" = the$number_of_states,
              "human.snv.rate.per.base" = the$human.snv.rate.per.base,
              "epi_fwd" = the$epi_fwd,
              "epi_bck" = the$epi_bck,
              "all_possible_states" = the$all_possible_states,
              "state_transition_mtx" = the$state_transition_mtx,
              "state_transition_vec_out" = the$state_transition_vec_out,
              "length_change_comb_list" = length(the$change_combs_list),
              "df_rate" = the$df_rate,
              "Q" = the$Q,
              "useQ" = the$useQ))
}

#' Set memento settings
#'
#' @param new_number_of_states integer value for number of states: 3 or 4 (requires phased data)
#' @param human.snv.rate.per.base human snv rate per base (usually 20/3e9)
#' @param new_epi_fwd numeric value for forward rate (gain of methylation): scales the human SNV rate, e.g. 1e5 times greater than SNV rate per base
#' @param new_epi_bck numeric value for backward rate (loss of methylation): scales the human SNV rate
#' @param new_rates data.frame of new rates
#' @param new_Q new Q matrix
#' @param new_useQ logical whether to use Q matrix
#'
#' @export
set_memento_settings <- function(new_number_of_states = NULL,
                                 human.snv.rate.per.base = NULL,
                                 new_epi_fwd = NULL,
                                 new_epi_bck = NULL,
                                 new_rates = NULL,
                                 new_Q = NULL,
                                 new_useQ = NULL){
  old_num_states <- the$number_of_states
  if(!is.null(new_number_of_states)){
    the$number_of_states <- new_number_of_states
    the$all_possible_states <- c(1:the$number_of_states)
    the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                       the$number_of_states,
                                       the$number_of_states,
                                       byrow = TRUE)
    colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

    the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

    the$change_combs_list <- lapply(c(0:10), function(nc){
      cc <- get_change_combinations_4s(nc)
      return(cc)
    })
    names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))
  }

  if(!is.null(human.snv.rate.per.base)){
    the$human.snv.rate.per.base <- human.snv.rate.per.base
  }
  if(!is.null(new_epi_fwd)){
    the$epi_fwd <- new_epi_fwd
  }
  if(!is.null(new_epi_bck)){
    the$epi_bck <- new_epi_bck
  }
  if(!is.null(new_rates)){
    the$df_rate <- new_rates
    the$Q <- getQMatrix(the$df_rate)
  }
  if(!is.null(new_useQ)){
    the$useQ <- new_useQ
  }
  if(!is.null(new_Q)){
    the$df_rate <- NULL
    the$Q <- new_Q
    cat("New Q matrix, but df_rate is NULL")
  }
  invisible(old_num_states)
}

#' Set memento settings to default
#'
#' @export
#'
#' @examples
#' set_memento_settings_to_default()
#' memento_settings()
set_memento_settings_to_default <- function(){
  the$number_of_states <- 3
  the$human.snv.rate.per.base=20/3e9
  the$epi_fwd = 1e5
  the$epi_bck = 1e5
  the$all_possible_states <- c(1:the$number_of_states)
  the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                     the$number_of_states,
                                     the$number_of_states,
                                     byrow = TRUE)
  colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

  the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

  the$change_combs_list <- lapply(c(0:10), function(nc){
    cc <- get_change_combinations_4s(nc)
    return(cc)
  })
  names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))

  the$df_rate <- getRateDf(the$number_of_states)
  the$df_rate$rate <- rep(0.001, nrow(the$df_rate))

  the$Q <- getQMatrix(the$df_rate)
  the$useQ <- FALSE
}

