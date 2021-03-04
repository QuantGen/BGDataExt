\name{trimSegments}
\alias{trimSegments}
\title{Trim segments}
\description{
    Remove variants that are not internal to the segment containing
    GWAS-significant variants.
}
\usage{
trimSegments(segments, statistic, chr, bp, threshold)
}
\arguments{
    \item{segments}{
        Segments determined by \code{getSegments}.
    }
    \item{statistic}{
        The statistic used for \code{segments}.
    }
    \item{chr}{
        The chromosomes used for \code{segments}.
    }
    \item{bp}{
        The base-pair positions used for \code{segments}.
    }
    \item{threshold}{
        The threshold to determine 'significance' (e.g., \code{1e-5} for
        p-values) used for \code{segments}.
    }
}
\value{
    Trimmed segments, see \link{\code{getSegments}}.
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
    statistic = pValues[, 4],
    chr = map(bg)$chromosome,
    bp = map(bg)$base_pair_position,
    threshold = 1e-5,
    gap = 1e6,
    verbose = FALSE
)

# Trims segments
trimmedSegments <- trimSegments(
    segments = segments,
    statistic = pValues[, 4],
    chr = map(bg)$chromosome,
    bp = map(bg)$base_pair_position,
    threshold = 1e-5
)
}
