#' Is the \code{rowRanges} empty?
#'
#' Check the \code{\link[SummarizedExperiment]{rowRanges}} of a \link[SummarizedExperiment]{RangedSummarizedExperiment} is empty, 
#' i.e., a \link[GenomicRanges]{GRangesList} with no ranges. 
#' 
#' @param x A \link[SummarizedExperiment]{RangedSummarizedExperiment} object or the contents of its \code{rowRanges}.
#'
#' @return A logical scalar indicating whether \code{x} has empty \code{rowRanges}.
#'
#' @details
#' Metadata in \code{\link[S4Vectors]{mcols}} is ignored for the purpose of this discussion, 
#' as this can be moved to the \code{\link[SummarizedExperiment]{rowData}(x)} of the base SummarizedExperiment class without loss.
#' In other words, non-empty \code{mcols} will not be used to determine that the \code{rowRanges} is not empty.
#' However, non-empty fields in the \code{\link[S4Vectors]{metadata}} or in the inner \code{mcols} of the \link[GenomicRanges]{GRanges} will trigger a non-emptiness decision.
#' 
#' @examples
#' emptyRowRanges(SummarizedExperiment())
#' emptyRowRanges(SummarizedExperiment(rowRanges=GRanges()))
#' emptyRowRanges(SummarizedExperiment(rowRanges=GRangesList()))
#' @export
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges PartitioningByEnd
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols<- mcols
#' @importFrom BiocGenerics relist
emptyRowRanges <- function(x) {
    if (is(x, "SummarizedExperiment")) {
        x <- rowRanges(x)
    }
    if (is(x, "GRangesList")) {
        # Creating an empty GRL and comparing it. This re-uses the same
        # code as the SE->RSE coerce method, more or less.
        partitioning <- PartitioningByEnd(integer(length(x)), names = names(x))
        rowRanges <- relist(GRanges(), partitioning)
        mcols(x) <- mcols(rowRanges, use.names = FALSE)
        identical(x, rowRanges)
    } else {
        FALSE
    }
}

