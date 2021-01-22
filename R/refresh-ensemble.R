#' Refresh a model ensemble by replacing low performing members.
#' 
#' Identify the lowest performing members of an ensemble and replace them with 
#' random draws from the original MCMC.
#' 
#' @param modlist List of filter model fits by county.
#' @param mcmcdir Name of the directory containing the MCMC inputs
#' @param frac Fraction of ensemble members to refresh (default = 0.05)
#' @param replacement_ensemble_override List of replacement ensembles. This is a
#' debugging parameter and will be removed.
#' @export
refresh_ensemble <- function(modlist, mcmcdir, frac=0.05, 
                             replacement_ensemble_override = NULL)
{
  if(is.null(replacement_ensemble_override)) 
  {  
    ## First, get the replacement ensemble members.
    nens <- length(modlist[[1]]$ids)
    nrepl <- as.integer(ceiling(nens*frac))
    rndselect <- withr::with_seed(NULL, runif(nrepl))
    datemax <- max(modlist[[1]]$obsdata$date)
    
    newensembles <- lapply(modlist,
                           function(countymod) {
                             locality <- countymod$obsdata$locality[1]
                             message('Running ', locality)
                             filter_fit_locality(locality, mcmcdir, nrepl, rndselect,
                                                 datemax=datemax)
                           })
  }
  else {
    nrepl <- length(replacement_ensemble_override[[1]]$ids)
    newensembles <- lapply(replacement_ensemble_override,
                           function(neloc) {
                             bogusids <- unique(neloc$history$id)
                             ids <- seq_along(bogusids)
                             row.names(neloc$finalstate) <- paste0('result.', ids)
                             neloc$history$id <- match(neloc$history$id, bogusids)
                             stopifnot(all(neloc$history$id %in% ids))
                             neloc
                           })
  }
  ## Determine which ensemble members are going to be replaced
  eeval <- ensemble_eval(modlist)
  replid <- eeval$id[seq(1, nrepl)]    # These are ordered by logl
  
  ## For each locality, splice the new ensemble members in in place of the old.
  ## First, assign the new ensemble members the id numbers of the ensemble members
  ## being replaced.
  idtbl <- data.frame(id = seq(1,nrepl), newid = replid)
  newensembles <- lapply(newensembles,
                         function(neloc) {
                           row.names(neloc$finalstate) <- paste0('result.', replid)
                           hist <- dplyr::left_join(neloc$history, idtbl, by='id')
                           hist$id <- hist$newid
                           hist$newid <- NULL
                           neloc$history <- hist
                           neloc
                         })
  
  ## Now splice these new ensemble members into the old ensemble.  The finalstate
  ## and history objects are what need to change.
  lapply(seq_along(modlist), 
         function(i) {
           countymod <- modlist[[i]]
           newens <- newensembles[[i]]
           countymod$finalstate[replid, ] <- newens$finalstate
           newhist <- countymod$history[!countymod$history$id %in% replid, ]
           newhist <- dplyr::bind_rows(newhist, newens$history)
           countymod$history <- dplyr::arrange(newhist, time, id)
           countymod
         })
}
