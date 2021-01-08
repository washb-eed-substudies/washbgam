

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
