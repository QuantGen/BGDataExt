\name{getSegments}
\alias{getSegments}
\title{Find non-overlapping segments based on a summary statistic}
\description{
    Given a summary statistic \code{x} and a threshold, this function computes
    the number of non-overlapping segments, each defined as a discovery (i.e.,
    \code{x[i] <= threshold)} +/- a gap, in the same units as \code{bp} (often
    base-pair position).
}
\usage{
getSegments(x, chr, bp, names, threshold, gap, trim = FALSE, verbose = FALSE)
}
\arguments{
    \item{x}{
        A statistic (e.g., BFDR or p-values).
    }
    \item{chr}{
        A vector containing the chromosome for each value of \code{x}.
    }
    \item{bp}{
        A vector containing the base-pair positions for each value of \code{x}.
    }
    \item{names}{
        A vector containing the variant names for each value of \code{x}.
    }
    \item{threshold}{
        The threshold to determine 'significance' (e.g., \code{1e-5} for
        p-values).
    }
    \item{gap}{
        1/2 of the length of the desired segments.
    }
    \item{trim}{
        Whether to collapse segments that were artifically inflated by
        \code{gap}. Defaults to \code{FALSE}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
}
\value{
    A data frame containing the following information:
    \item{chr}{
        Chromosome
    }
    \item{start}{
        Index where segment starts within \code{x}.
    }
    \item{end}{
        Index where segment ends within \code{x}.
    }
    \item{length}{
        Length of segment.
    }
    \item{bpStart}{
        Base-pair position where segment starts.
    }
    \item{bpEnd}{
        Base-pair position where segment ends.
    }
    \item{bpLength}{
        Length of segment in base-pair positions.
    }
    \item{minName}{
        Name of variant with the smallest value of \code{x} within segment.
    }
    \item{minValue}{
        Smallest value of \code{x} within segment.
    }
}
\examples{
library(BGData)

# Load example data
bg <- BGData:::loadExample()

# Perform GWAS
pValues <- GWAS(
    formula = FT10 ~ 1,
    data = bg,
    method = "rayOLS"
)

# Determine segments within +/- 1MB from a significant variant
segments <- getSegments(
    x = pValues[, 4],
    chr = map(bg)$chromosome,
    bp = map(bg)$base_pair_position,
    names = map(bg)$snp_id,
    threshold = 1e-5,
    gap = 1e6,
    trim = FALSE,
    verbose = FALSE
)
}
