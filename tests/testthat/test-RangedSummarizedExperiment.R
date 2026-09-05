# library(testthat); library(alabaster.se); source("test-RangedSummarizedExperiment.R")

set.seed(100)

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SummarizedExperiment(list(counts=mat, cpm=mat/10), rowRanges=GRanges("chrA", IRanges(1:200, width=1)))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

test_that("saveObject works as expected for RSE objects", {
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(colData(se), colData(out2))
    expect_identical(rowRanges(se), rowRanges(out2))
    expect_true(file.exists(file.path(tmp, "row_ranges", "OBJECT")))
    expect_false(file.exists(file.path(tmp, "row_ranges", "range_annotations")))
})

test_that("saveObject preserves RSE rownames", {
    copy <- se
    rownames(copy) <- sprintf("GENE_%i", seq_len(nrow(copy)))

    tmp <- tempfile()
    saveObject(copy, tmp)
    out2 <- readObject(tmp)

    expect_identical(rownames(out2), rownames(copy))
    expect_identical(rowRanges(out2), rowRanges(copy))
    expect_identical(rowData(out2), rowData(copy))
})

test_that("saveObject auto-skips on empty rowRanges", {
    # GRL is empty.
    copy <- GRangesList(rep(list(GRanges()), nrow(se)))
    names(copy) <- rownames(se)
    rowRanges(se) <- copy
    expect_true(emptyRowRanges(se))

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "RangedSummarizedExperiment")
    expect_true(emptyRowRanges(out2))
    expect_false(file.exists(file.path(tmp, "row_ranges")))

    # Non-empty rowData but GRL is still empty.
    mcols(copy)$FOO <- 2
    rowRanges(se) <- copy

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "RangedSummarizedExperiment")
    expect_true(emptyRowRanges(out2))
    expect_identical(unique(rowData(out2)$FOO), 2)
    expect_false(file.exists(file.path(tmp, "row_ranges")))
})

test_that("saveObject handles GRLs", {
    rowRanges(se) <- splitAsList(rowRanges(se), seq_len(nrow(se)))
    rownames(se) <- paste0("FEATURE_", seq_len(nrow(se)))

    tmp <- tempfile()
    saveObject(se, tmp)
    out <- readObject(tmp)
    expect_identical(rowRanges(se), rowRanges(out))
})
