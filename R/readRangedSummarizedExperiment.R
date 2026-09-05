#' Read a RangedSummarizedExperiment from disk
#'
#' Read a \link[SummarizedExperiment]{RangedSummarizedExperiment} from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link[alabaster.base]{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link[alabaster.base]{saveObject}} method for RangedSummarizedExperiment objects.
#' @param metadata Named list of metadata for this object, see \code{\link[alabaster.base]{readObjectFile}} for details.
#' @param ... Further arguments passed to \code{\link{readSummarizedExperiment}} and internal \code{\link[alabaster.base]{altReadObject}} calls.
#' 
#' @return A RangedSummarizedExperiment object.
#'
#' @author Aaron Lun
#' @seealso
#' \code{"\link{saveObject,RangedSummarizedExperiment-method}"}, to save the RangedSummarizedExperiment to disk.
#'
#' @examples
#' # Mocking up an experiment:
#' mat <- matrix(rpois(10000, 10), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#'
#' gr <- GRanges("chrA", IRanges(1:1000, width=10))
#' se <- SummarizedExperiment(list(counts=mat), rowRanges=gr)
#' se$stuff <- LETTERS[1:10]
#' rowData(se)$blah <- runif(1000)
#' metadata(se)$whee <- "YAY"
#' 
#' tmp <- tempfile()
#' saveObject(se, tmp)
#' readObject(tmp)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<- rowRanges<-
#' @import alabaster.base
readRangedSummarizedExperiment <- function(path, metadata, ...) {
    # We don't try to respect application overrides when loading the base
    # instance. Application developers should just pretend that we copied the
    # code from readSummarizedExperiment, rather than trying to inject in
    # custom code at this point, which gets too complicated - see the
    # associated commentary for saveObject,SummarizedExperiment-method.
    se <- readSummarizedExperiment(path, metadata=metadata, ...)

    rrdir <- file.path(path, "row_ranges")
    if (file.exists(rrdir)) {
        # Avoid overriding the old rowData with the rowRanges's mcols.
        rr <- altReadObject(rrdir, ...)
        old.rd <- rowData(se)
        old.names <- rownames(se)
        rowRanges(se) <- rr
        rowData(se) <- old.rd
        rownames(se) <- old.names
    } else {
        se <- as(se, "RangedSummarizedExperiment")
    }

    se
}
