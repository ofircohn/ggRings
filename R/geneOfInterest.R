#' Example gene-level dataset for circular ring plotting
#'
#' A toy dataset with 10 genetic variants and associated signals.
#'
#' @format A data frame with 10 rows and 4 columns:
#' \describe{
#'   \item{rs}{SNP ID (e.g., "rs10127495")}
#'   \item{ATAC}{ATAC-seq signal (0–1)}
#'   \item{R2}{Linkage disequilibrium (R-squared)}
#'   \item{H3K27ac}{H3K27ac ChIP-seq signal (0–1)}
#' }
#'
#' @examples
#' data(geneOfInterest)
#' head(geneOfInterest)
"geneOfInterest"
