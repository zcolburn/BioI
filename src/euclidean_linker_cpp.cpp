#include <Rcpp.h>
using namespace Rcpp;

#include<cmath>


//' @title Return the group number for each localization.
//'
//' @description
//' Group PALM/iPALM localizations based on their physical separation distance
//'
//' PALM/iPALM data results in a list of spatial coordinates for fluorophore
//' localizations. This function groups nearby localizations if they are within
//' the provided critical distance from each other.
//'
//' @param input A numeric matrix where each row is a localization and each
//' column is a spatial axis.
//' @param critDist The critical distance for which localizations nearer than
//' this distance are deemed part of the same group.
//'
//' @author Zach Colburn
//'
//' @examples
//' # Function call
//' \dontrun{euclidean_linker_cpp(inputMatrix, critDist)}
//'
//' @import Rcpp
//'
//' @useDynLib Bioi, .registration = TRUE
//[[Rcpp::export(.euclidean_linker_cpp)]]
Rcpp::NumericVector euclidean_linker_cpp(
    Rcpp::NumericMatrix input,
    double critDist
) {
  // Initialize variables
  //
  // Get squared critical distance
  double scd = pow(critDist,2);

  // Initialize output group vector
  Rcpp::LogicalVector inGroup(input.nrow(),false);
  Rcpp::NumericVector output(input.nrow(),0.0);

  float gn=0;
  for(float i=0;i<output.length();i++){
    output(i)=gn;
    gn++;
  }

  // Get the number of dimensions to evaluate
  float nDim=input.ncol();

  // Get the number of points to evaluate
  float np=input.nrow();

  // Iterate through points
  //
  // fp = first point
  // sp = second point
  // np = number of points to evaluate
  // nDim = number of dimensions to evaluate
  for(float fp=0;fp<np-1;fp++){
    for(float sp=fp+1;sp<np;sp++){
      // Get the squared distance between the two points
      double sd = 0;
      for(float d=0;d<nDim;d++){
        sd += pow(input(fp,d)-input(sp,d),2);
        if(sd > scd){
          break;
        }
      }
      // If the squared distance between the two points is less than or equal to
      // the squared critical distance (scd, see above), then perform group
      // assignment.
      if(sd <= scd){
        if(!inGroup(fp) && !inGroup(sp)){
          inGroup(fp)=true;
          inGroup(sp)=true;
          gn++;
          output(fp)=gn;
          output(sp)=gn;
        }
        else if(!inGroup(fp) && inGroup(sp)){
          output(fp)=output(sp);
          inGroup(fp)=true;
        }
        else if(inGroup(fp) && !inGroup(sp)){
          output(sp)=output(fp);
          inGroup(sp)=true;
        }
        else if(output(fp) != output(sp)){
          for(float i=0;i<np;i++){
            if(output(i)==output(sp)){
              output(i)=output(fp);
            }
          }
        }
      }

    }
  }

  return output;
}
