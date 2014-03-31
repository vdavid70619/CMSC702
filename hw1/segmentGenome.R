# implemets genome segmentation
# arguments:
#   y: state assignments, vector of 1 and 0
#   L: segment length used to estimate states
#   minSegmentLength: minimum genomic length of segments returned
#
# returns:
#   GRanges object for contiguous segments where y=1 with length greater than minSegmentLength
segmentGenome <- function(y, L, minSegmentLength = 200, seqnames="chr16") {
  seg_cand = Rle(y)
  seg_start = start(seg_cand)[seg_cand@values==1 & seg_cand@lengths>minSegmentLength/L]
  seg_end = end(seg_cand)[seg_cand@values==1 & seg_cand@lengths>minSegmentLength/L]
  return(GRanges(seqnames=seqnames, ranges=IRanges(start=(seg_start-1)*L+1,end=seg_end*L)))
}
