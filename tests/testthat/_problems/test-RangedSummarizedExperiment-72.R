# Extracted from test-RangedSummarizedExperiment.R:72

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "alabaster.se", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
set.seed(100)
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))
se <- SummarizedExperiment(list(counts=mat, cpm=mat/10), rowRanges=GRanges("chrA", IRanges(1:200, width=1)))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))
rowRanges(se) <- splitAsList(rowRanges(se), seq_len(nrow(se)))
rownames(se) <- paste0("FEATURE_", seq_len(nrow(se)))

# test -------------------------------------------------------------------------
tmp <- tempfile()
saveObject(se, tmp)
out2 <- readObject(tmp)
expect_identical(rowRanges(se), rowRanges(out))
