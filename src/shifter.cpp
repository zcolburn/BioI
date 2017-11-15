#include <Rcpp.h>
using namespace Rcpp;

//' Identify the local shift between two temporally separated images.
//'
//' This function determines the local shift between two temporally separated
//' images. It does so by minimizing the squared difference between a local
//' template centered on the supplied peg coordinates and adjacent matrices of
//' the same dimension.
//'
//' @title Identify the local shift between two temporally separated images.
//'
//' @param t0 A numeric matrix of intensity values for the first image.
//' @param t1 A numeric matrix of intensity values for the second image.
//' @param xSeq A numeric vector of the x coordinates defining template centers.
//' @param ySeq A numeric vector of the y coordinates defining template centers.
//' @param searchBoxRadius An integer defining the search box radius.
//' @param matchBoxRadius An integer defining the radius of templates.
//' @param shiftSkipper An integer defining the number of pixels away the next
//' template to image comparison will be. This value must be 1 or greater.
//' Larger values will dramatically reduce execution time, but can reduce data
//' quality.
//' @author Zach Colburn
//' @examples
//' # Function call
//' shifter(t0, t1, xSeq, ySeq, searchBoxRadius, matchBoxRadius, shiftSkipper)
//'
// [[Rcpp::export]]
List shifter(
    NumericMatrix t0,
    NumericMatrix t1,
    NumericVector xSeq,
    NumericVector ySeq,
    int searchBoxRadius,
    int matchBoxRadius,
    int shiftSkipper)
{
  NumericMatrix shiftX(t0.ncol(),t1.nrow());
  NumericMatrix shiftY(t0.ncol(),t1.nrow());
  xSeq = xSeq - 1;
  ySeq = ySeq - 1;
  int boxSize = pow(2 * searchBoxRadius + 1, 2);

  for (int xPegIndex = 0; xPegIndex < xSeq.size(); xPegIndex++) {
    int xPeg = xSeq[xPegIndex];
    for (int yPegIndex = 0; yPegIndex < ySeq.size(); yPegIndex++) {
      int yPeg = ySeq[yPegIndex];
      // Below here: x and y peg positions are denoted by xPeg and yPeg

      // Define template matrix
      NumericMatrix tMat = t0(Range(yPeg - searchBoxRadius,yPeg + searchBoxRadius),Range(xPeg - searchBoxRadius,xPeg + searchBoxRadius));

      // Define optimization parameters
      int bestX = 0;
      int bestY = 0;

      for (int relativeX = -matchBoxRadius+searchBoxRadius; relativeX <= matchBoxRadius-searchBoxRadius; relativeX = relativeX + shiftSkipper) {
        int shiftedXCenter = relativeX + xPeg;
        if ((shiftedXCenter - searchBoxRadius) < 0) {
          continue;
        }
        if ((shiftedXCenter + searchBoxRadius) > shiftX.ncol()) {
          continue;
        }
        for (int relativeY = -matchBoxRadius+searchBoxRadius; relativeY <= matchBoxRadius-searchBoxRadius; relativeY = relativeY + shiftSkipper) {
          int shiftedYCenter = relativeY + yPeg;
          if ((shiftedYCenter - searchBoxRadius) < 0) {
            continue;
          }
          if ((shiftedYCenter + searchBoxRadius) > shiftY.nrow()) {
            continue;
          }
          // Below here: x and y overlay centers are denoted by shiftedXCenter and shiftedYCenter

          // Define overlay matrix
          NumericMatrix oMat = t1(Range(shiftedYCenter - searchBoxRadius,shiftedYCenter + searchBoxRadius),Range(shiftedXCenter - searchBoxRadius,shiftedXCenter + searchBoxRadius));

          // Find differences
          float sum = 0.0;
          for (int i = 0; i < boxSize; i++) {
            float difference = oMat[i] - tMat[i];
            if (difference < 0.0) {
              difference = difference * -1;
            }
            sum = sum + difference;
          }

          // Above here: x and y overlay centers are denoted by shiftedXCenter and shiftedYCenter
        }
      }
      // Assign optimized parameters to output
      shiftX(yPeg,xPeg) = bestX;
      shiftY(yPeg,xPeg) = bestY;

      // Above here: x and y peg positions are denote by xPeg and yPeg
    }
  }

  List output;
  output["x"] = shiftX;
  output["y"] = shiftY;
  return output;
}
