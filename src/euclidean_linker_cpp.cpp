#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>


// #include <fstream>
// #include <string>
// #include <iostream>


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
//' @param use_prog_bar A logical indicating whether a progress bar should be
//' used. This must be set to false when running in parallel.
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
    double critDist,
    bool use_prog_bar = true
) {
  //std::ofstream logf("log.txt");
  //logf << "Creating log file.\n";

  // Initialize variables.
  //
  // Get squared critical distance.
  double scd = pow(critDist,2);
  //logf << "scd: " << scd << "\n";

  // Get the number of dimensions to evaluate.
  int nDim=input.ncol();
  //logf << "nDim: " << nDim << "\n";

  // Get the number of points to evaluate.
  int np=input.nrow();
  //logf << "np: " << np << "\n";

  // Initialize group_array (g_array).
  int* g_array = new int[np];
  for(int i=0;i<np;i++){
    g_array[i]=0;
  }

  // Initialize first_in_group_array (fig_array).
  int* fig_array = new int[np];
  for(int i=0;i<np;i++){
    fig_array[i]=-1;
  }

  // Create progress bar
  int num_bar_elements=50;
  int current_num_bars = 0;
  if(use_prog_bar){
    Rcpp::Rcout << "||";
    for(int i=0;i<num_bar_elements;i++){
      Rcpp::Rcout << "=";
    }
    Rcpp::Rcout << "||\n";
  }

  // Iterate through points.
  //
  // fp = first point
  // sp = second point
  // np = number of points to evaluate
  // nDim = number of dimensions to evaluate
  // gn = group number
  //logf << "Entering main loop\n";
  if(use_prog_bar){
    Rcpp::Rcout << "||";
  }
  int gn=1;
  for(int fp=0;fp<np-1;fp++){
    Rcpp::checkUserInterrupt();
    //logf << "============= New fp =============\n";
    //logf << "fp: " << fp << "\n";
    int sp=fp+1;
    if(sp == np){break;}
    // Compare the fp to multiple sp, but only as many as make sense to compare.
    // Figure out which is the maximum sp to compare.
    if(sp != np){
      int max_sp=sp;
      //logf << "Starting 'while' 1.\n";
      for(int i=sp;i<np;i++){
        if(!(abs(input(fp,0)-input(max_sp,0)) <= critDist)){
          break;
        }else{
          max_sp++;
        }
      }
      //logf << "max_sp: " << max_sp << "\n";
      //logf << "max_sp-sp+1: " << max_sp-sp+1 << "\n";

      // Initialize all elements of sp_indices to the corresponding sp.
      int* sp_indices = new int[max_sp-sp+1];
      //logf << "sp_indices created.\n";
      for(int i=0;i<(max_sp-sp+1);i++){
        sp_indices[i]=sp+i;
      }
      //logf << "sp_indices populated.\n";
      // Send to log
      // for(int k=0;k<(max_sp-sp+1);k++){
      //   //logf << sp_indices[k] << " ";
      // }
      //logf << "\n";

      // Determine how many sp are valid for this fp.
      int num_valid_sp=0;
      for(int i=0;i<(max_sp-sp+1);i++){
        if(sp_indices[i] == -1){break;}
        num_valid_sp++;
      }
      //logf << "num_valid_sp: " << num_valid_sp << "\n";

      // If there are no valid sp, then set a group number for this fp and skip
      // to the next fp.
      if(num_valid_sp == 0){
        //logf << "num_valid_sp is 0.\n";
        g_array[fp]=gn;
        gn++;
        delete[] sp_indices;
        continue;
        //logf << "Finished num_valid_sp is 0 section.\n";
      }

      // Determine the index of the max_sp for this fp.
      int max_sp_index = sp_indices[num_valid_sp-1];

      // For each valid sp, determine whether it is close enough to fp for
      // linkage. If it is then link it.
      //
      // csp = current second point
      int i=0;
      //logf << "Entering main 'while'.\n";
      while(
        (i < (max_sp-sp+1)) &&
          (sp_indices[i] != -1)
      ){
        //logf << "Inside main 'while'.\n";
        // Get squared distance between fp and the given sp in sp_indices.
        int csp=sp_indices[i];
        //logf << "csp: " << csp << "\n";
        double sd = 0;
        for(int d=0;d<nDim;d++){
          sd += pow(input(fp,d)-input(csp,d), 2);
          if(sd > scd){
            break;
          }
        }
        if(sd > scd){
          //logf << "Point not linked.\n";
          i++;
          continue;
        }
        // Perform linkage.
        if(sd <= scd){
          //logf << "Performing linkage.\n";
          if((g_array[fp] == 0) && (g_array[csp] == 0)){
            //logf << "Starting type 1 linkage.\n";
            g_array[fp]=gn;
            g_array[csp]=gn;
            //logf << "Updated group numbers.\n";
            fig_array[gn]=fp;
            //logf << "Updated fig_array.\n";
            gn++;
            //logf << "Completed type 1 linkage.\n";
          }
          else if((g_array[csp] != 0) && (g_array[fp] == 0)){
            //logf << "Starting type 2 linkage.\n";
            g_array[fp]=g_array[csp];
            //logf << "Completed type 2 linkage.\n";
          }
          else if((g_array[fp] != 0) && (g_array[csp] == 0)){
            //logf << "Starting type 3 linkage.\n";
            g_array[csp]=g_array[fp];
            //logf << "Completed type 3 linkage.\n";
          }
          else if(g_array[fp] != g_array[csp]){
            //logf << "Starting type 4 linkage.\n";
            int fp_gn=fig_array[fp];
            int sp_gn=fig_array[csp];
            int start=fp_gn;
            int old_g=g_array[csp];
            int new_g=g_array[fp];
            if(sp_gn < start){
              start=sp_gn;
              old_g=g_array[fp];
              new_g=g_array[csp];
            }
            for(int j=start;j<=max_sp_index;j++){
              if(g_array[j] == old_g){
                g_array[j]=new_g;
              }
            }
            //logf << "Completed type 4 linkage.\n";
          }
        }

        // Move to the next index in sp_indices.
        i++;
      }

      delete[] sp_indices;

    }else{
      g_array[fp]=gn;
    }

    // Account for the case in which no links were made. Assign fp to the next
    // group number.
    if(g_array[fp] == 0){
      g_array[fp]=gn;
      gn++;
    }

    // Progress bar.
    if(use_prog_bar){
      int bars_needed = ((float)(fp)) / ((float)(np)) * num_bar_elements;
      if((bars_needed - current_num_bars) != 0){
        for(int iter=0;iter<(bars_needed-current_num_bars);iter++){
          Rcpp::Rcout << ".";
          current_num_bars++;
        }
      }
    }
  }

  // Account for possibility that the last point is not in a group.
  if(g_array[np-1] == 0){
    g_array[np-1]=gn;
  }

  // Finish filling progress bar.
  if(use_prog_bar){
    if(current_num_bars < num_bar_elements){
      for(int i=0;i<num_bar_elements-current_num_bars;i++){
        Rcpp::Rcout << ".";
      }
    }
    Rcpp::Rcout << "||";
    Rcpp::Rcout << "\n";
  }

  //logf << "Escaped main loop!";

  // Initialize output group vector.
  Rcpp::NumericVector output(input.nrow());

  // Set values of output.
  for(int i=0;i<np;i++){
    output(i)=g_array[i];
  }

  // Clean up.
  delete[] g_array;
  delete[] fig_array;

  // Return output.
  return output;
}
