#include <Rcpp.h>
using namespace Rcpp;

//' Group continuous, thresholded objects in rasterized objects.
//'
//' This function takes a matrix corresponding to a thresholded image and
//' returns a matrix of the same size, where all adjacent, thesholded pixels
//' are the same integer corresponding to that object's cluster number.
//'
//' @title Assign all neighboring pixels the same group number.
//'
//' @param img A thresholded matrix (where non-object pixels are assigned a
//' value of 0).
//' @param pixRange An integer number of pixels. Execution is faster when this
//' value is small. However, the value must be larger than the diameter of the
//' largest continuous object in the image.
//' @author Zach Colburn
//' @examples
//' # Function call
//' \dontrun{identify_thresholded_objects(img, pixRange)}
//' @export
// [[Rcpp::export]]
NumericMatrix identify_thresholded_objects(NumericMatrix img, int pixRange) {
  // Get the image dimensions
  int numRows = img.nrow();
  int numCols = img.ncol();

  // Assign all non-zero values a new number
  int gN = 1;
  for (int y = 0; y < numRows; y++){
    for (int x = 0; x < numCols; x++){
      if (img(y,x) != 0) {
        img(y,x) = gN;
        gN++;
      }
    }
  }

  // Iterate through the pixels in the image
  NumericVector xG = NumericVector::create(-1,-1,0,1);
  NumericVector yG = NumericVector::create(0,-1,-1,-1);
  for(int y=1; y<numRows; y++){
    for(int x=1; x<(numCols-1); x++){
      // If this pixel belongs to a thresholded object...
      if(img(y,x) != 0){
        bool doneF=false;
        int fixTo=0;
        // Check its neighbors to see if they belong to an object as well.
        for(int i=0; i<4; i++){
          int nX=xG[i]+x;
          int nY=yG[i]+y;
          int nG=img(nY,nX);
          // If the neighbor belongs to an object...
          if(nG > 0){
            // Check if a matched neighbor has already been converted...
            if(!doneF){
              // If no such conversion has been performed, then convert the
              // present pixel rather than its neighbors.
              doneF=true;
              fixTo=nG;
              img(y,x)=fixTo;
            }else{
              // If a conversion has been performed then convert other neighbors
              // to the the same value the previous conversion made.

              // Determine the bounds of the conversion region
              int ly = -pixRange+y;
              if(ly<0){ly=0;}
              int lx = -pixRange+x;
              if(lx<1){lx=1;}
              int hx = pixRange+x;
              if(hx>(numCols-1)){hx=numCols-1;}

              // Perform conversion of groups in the conversion region
              // identified above.
              for(int fY=ly;fY<y;fY++){
                for(int fX=lx;fX<hx;fX++){
                  if(img(fY,fX) == nG){
                    img(fY,fX) = fixTo;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Return a matrix in which all thresholded pixels are replaced with their
  // group number.
  return img;
}
