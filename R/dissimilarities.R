pr_Euclidean <- function(x, y) sqrt(crossprod(x - y))
pr_Euclidean_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Euclidean"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_euclidean_dist",
                names = c("Euclidean","L2"),
                PREFUN = "pr_Euclidean_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sqrt(sum_i (x_i - y_i)^2))",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Euclidean Distance (C implementation with compensation for excluded components)")

pr_Mahalanobis <- function(x, y, cov) sqrt(mahalanobis(x, y, cov))
pr_Mahalanobis_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (length(p) < 1) p <- list(cov(x, y))
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "pr_Mahalanobis",
                names = "Mahalanobis",
                PREFUN = "pr_Mahalanobis_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt((x - y) Sigma^(-1) (x - y))",
                reference = "Mahalanobis P.C. (1936), On the generalised distance in statistics, Proceedings of the National Institute of Science of India 12, pp. 49-55",
                description = "The Mahalanobis Distance. The Variance-Covariance-Matrix is estimated from the input data if unspecified.")

pr_Bhjattacharyya <- function(x, y) sqrt(crossprod(sqrt(x) - sqrt(y)))
pr_DB$set_entry(FUN = "pr_Bhjattacharyya",
                names = "Bhjattacharyya",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt(sum_i (sqrt(x_i) - sqrt(y_i))^2))",
                reference = "Bhattacharyya A. (1943). On a measure of divergence between two statistical populations defined by probability distributions, Bull. Calcutta Math. Soc., vol. 35, pp. 99--109",
                description = "The Bhjattacharyya Distance")

pr_Manhattan <- function(x, y) sum(abs(x - y))
pr_Manhattan_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Manhattan"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_manhattan_dist",
                names = c("Manhattan", "City-Block", "L1", "Taxi"),
                PREFUN = "pr_Manhattan_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Manhattan/City-Block/Taxi/L1-Distance (C implementation with compensation for excluded components)")

pr_supremum <- function(x, y) max(abs(x - y))
pr_supremum_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_supremum"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_maximum_dist",
                names = c("supremum", "max", "Tschebyscheff", "Chebyshev"),
                distance = TRUE,
                PREFUN = "pr_supremum_prefun",
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "max_i |x_i - y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Maximum/Supremum/Chebyshev Distance (C implementation)")

pr_Minkowski <- function(x, y, p = 2) (sum(abs(x - y) ^ p)) ^ (1/p)
pr_Minkowski_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (length(p) < 1)
        stop("Argument 'p' mandatory!")
    p <- p[[1]]
    if (p < 1)
        stop("p must not be smaller than 1.")
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Minkowski"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_minkowski_dist",
                names = c("Minkowski","Lp"),
                PREFUN = "pr_Minkowski_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "(sum_i (x_i - y_i)^p)^(1/p)",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Minkowski Distance (C implementation with compensation for excluded components)")

pr_Canberra <- function(x, y) {tmp <- abs(x - y) / abs(x + y); sum(tmp[!is.nan(tmp)])}
pr_Canberra_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Canberra"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_canberra_dist",
                names = "Canberra",
                PREFUN = "pr_Canberra_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / |x_i + y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Canberra Distance (C implementation with compensation for excluded components)")

pr_WaveHedges <- function(x, y) sum(1 - min(x, y) / max(x, y))
pr_DB$set_entry(FUN = "pr_WaveHedges",
                names = c("Wave", "Hedges"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i (1 - min(x_i, y_i) / max(x_i, y_i))",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Wave/Hedges Distance")

pr_Divergence <- function(x, y) {tmp <- (x - y)^2 / (x + y)^2; sum(tmp[!is.nan(tmp)])}
pr_DB$set_entry(FUN = "pr_Divergence",
                names = "Divergence",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i (x_i - y_i)^2 / (x_i + y_i)^2",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Divergence Distance")

pr_KullbackLeibler <- function(x, y) sum(x * log((x / sum(x) / (y / sum(y))) / sum(x)))
pr_DB$set_entry(FUN = "pr_KullbackLeibler",
                names = c("Kullback", "Leibler"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i [x_i * log((x_i / sum_j x_j) / (y_i / sum_j y_j)) / sum_j x_j)]",
                reference = "Kullback S., and Leibler, R.A. (1951). On information and sufficiency. The Annals of Mathematical Statistics, vol. 22, pp. 79--86",
                description = "The Kullback-Leibler-distance.")

pr_BrayCurtis <- function(x, y) sum(abs(x - y)) / sum(x + y)
pr_DB$set_entry(FUN = "pr_BrayCurtis",
                names = c("Bray","Curtis"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / sum_i (x_i + y_i)",
                reference = "Bray J.R., Curtis J.T. (1957). An ordination of the upland forest of the southern Winsconsin. Ecological Monographies, 27, pp. 325--349",
                description = "The Bray/Curtis Distance")

pr_Soergel <- function(x, y) sum(abs(x - y)) / sum(max(x, y))
pr_DB$set_entry(FUN = "pr_Soergel",
                names = "Soergel",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / sum_i max{x_i, y_i}",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Soergel Distance")

pr_Levenshtein_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!require("cba")) stop("Need package cba!")
    if (pairwise)
        stop("Pairwise distances not implemented by sdist()!")
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "sdists",
                names = "Levenshtein",
                PREFUN = "pr_Levenshtein_prefun",
                convert = "pr_dist2simil",
                distance = TRUE,
                loop = FALSE,
                abcd = FALSE,
                C_FUN = FALSE,
                formula = "Number of insertions, edits, and deletions between to strings",
                reference = "Levenshtein V.I. (1966). Binary codes capable of correcting deletions, insertions, and reversals. Soviet Physics Doklady 10, pp. 707--710",
                description = "Wrapper for sdists() in the cba-package (C implementation).")





