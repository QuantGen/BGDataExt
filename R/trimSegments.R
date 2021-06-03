trimSegments <- function(segments, statistic, threshold) {
    runStart <- segments[["start"]]
    runEnd <- segments[["end"]]
    for (curSeg in 1:nrow(segments)) {
        segFilter <- seq(runStart[curSeg], runEnd[curSeg])
        statisticSeq <- statistic[segFilter]
        # Determine which variants in the segment passed the threshold
        significantVariants <- which(statisticSeq <= threshold)
        # Set start of run to first significant variant and end of run to last
        # significant variant
        runStart[curSeg] <- segFilter[significantVariants[1]]
        runEnd[curSeg] <- segFilter[significantVariants[length(significantVariants)]]
    }
    segments[["start"]] <- runStart
    segments[["end"]] <- runEnd
    segments[["length"]] <- runEnd - runStart + 1
    return(segments)
}
