# library(testthat); library(alabaster.se); source("test-SummarizedExperiment.R")

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SummarizedExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

rownames(se) <- sprintf("GENE_%i", seq_len(nrow(se)))

test_that("saveObject works as expected for SE objects", {
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rownames(se), rownames(out2))
    expect_equal(rowData(se)[,1], rowData(out2)$whee)
    expect_equal(sum(assay(se)), sum(assay(out2)))
    expect_equal(sum(assay(se, 2)), sum(assay(out2, 2)))
})

test_that("saveObject works as expected with no row or column names", {
    dimnames(se) <- NULL

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("saveObject works as expected with no assays", {
    tmp <- tempfile()
    se <- SummarizedExperiment()
    saveObject(se, tmp)
    expect_identical(readObject(tmp), se)

    # Plus some dimensions, at least.
    tmp <- tempfile()
    se <- SummarizedExperiment(colData=DataFrame(row.names=LETTERS), rowData=DataFrame(row.names=1:10))
    saveObject(se, tmp)
    expect_identical(readObject(tmp), se)
})

test_that("saveObject works as expected with no row or column data, but still names", {
    tmp <- tempfile()
    dir.create(tmp)

    colData(se) <- colData(se)[,0]
    rowData(se) <- rowData(se)[,0]

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("saveObject works as expected with no row or column data at all", {
    tmp <- tempfile()
    dir.create(tmp)

    colData(se) <- colData(se)[,0]
    rowData(se) <- rowData(se)[,0]
    dimnames(se) <- list(NULL, NULL)

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
    expect_false(file.exists(file.path(tmp, "row_data")))
    expect_false(file.exists(file.path(tmp, "column_data")))
})

test_that("saveObject fails when the assay names are NULL or non-unique", {
    tmp <- tempfile()
    ass <- assays(se)
    names(ass) <- NULL
    assays(se) <- ass

    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "should be named")

    # Duplicated
    tmp <- tempfile()
    dir.create(tmp)
    assayNames(se) <- rep("FOO", length(assays(se)))

    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "unique")

    # Empty.
    tmp <- tempfile()
    dir.create(tmp)
    assayNames(se) <- c("", head(LETTERS, length(assayNames(se)) - 1))

    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "empty")
})

test_that("saveObject handles data frames in the assays", {
    tmp <- tempfile()
    dir.create(tmp)

    assay(se) <- as.data.frame(assay(se))
    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "should not contain data frames")

    # But now can handle DataFrame.
    se2 <- se
    assay(se2) <- as(assay(se2), "DataFrame")
    tmp <- tempfile()
    saveObject(se2, tmp, SummarizedExperiment.allow.dataframe.assay=TRUE)
    out2 <- readObject(tmp)
    expect_identical(assay(se2), assay(out2))
})

test_that("saveObject saves other metadata when necessary", {
    metadata(se)$YAY <- 1

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(metadata(out2), list(YAY=1))
})
