\name{getWindows}
\alias{getWindows}
\title{Computes sliding LD-based windows}
\description{
    For each SNP present in a genotpe matrix, obtains nearby SNPs in LD based
    on genotype correlations.
}
\usage{
getWindows(X, center_snp = 1:100, rSq = 0.1, maxGaps = 2, ids = 1:nrow(X))
}
\arguments{
    \item{X}{
        A matrix-like object, typically \code{@geno} of a \code{BGData} object.
    }
    \item{center_snp}{
        Columns of \code{X} to compute windows. By default, does not attempt to
        obtain LD windows surrounding each element of \code{X}.
    }
    \item{rSq}{
        R-squared threshold for determining genotype correlations.
    }
    \item{maxGaps}{
        The number of consequtive SNPs that may fall below \code{rSq} until LD
        is considered lost.
    }
    \item{ids}{
        Indices of individuals used to compute LD, if not all individuals in X
        are desired.
    }
}
\value{
    A list the length of the number of columns of \code{X} containing sliding
    LD-based windows. The elements of the list are named for each "center" SNP
    used in LD determination. Each element is a vector of SNP names--SNPs in LD
    with the "center" SNP.
}
