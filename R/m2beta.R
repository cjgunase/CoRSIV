#' Convert m values to beta values
#'
#' This function takes input as dataframe with columns of probe data and rows are samples
#' @param m m is a data frame with columns are probe data and rows are samples
#' @keywords mvalues betavalues
#' @export
#' @examples
#' m2beta(m)
#' @return returns a dataframe
m2beta <- function(m){return (2^m/(2^m+1))}
