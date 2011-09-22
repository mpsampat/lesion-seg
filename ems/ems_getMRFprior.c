/* ems_getMRFprior.c    Koen Van Leemput - August 17, 2001 */

#include <math.h>
#include "mex.h"
  
static void getMRFprior(unsigned char *classification, unsigned char *playing, 
                        unsigned int DIM[3], 
                        unsigned int nrOfClasses, unsigned char *atlas, 
                        unsigned int nrOfPriors, unsigned int *lkp, 
                        double *G, double *H, 
                        unsigned char *MRFprior) 
  {
  double *u, *v, *voxelprior, sumvoxelprior;
  unsigned int colSize, planeSize, blockSize, ind, plane, row, col,
    class, neighbourprior, Gind, testind, maxCol, maxRow, maxPlane, prior,
    *nrOfClassesPerPrior;


  /* How many classes exist that use the same prior? */
  nrOfClassesPerPrior = mxCalloc(nrOfPriors, sizeof(unsigned int));
  for (class=0; class<nrOfClasses; class++)
    {
    nrOfClassesPerPrior[lkp[class]]++;
    }


  /* Some initialization */
  u = mxCalloc(nrOfPriors, sizeof(double));
  v = mxCalloc(nrOfPriors, sizeof(double));
  voxelprior = mxCalloc(nrOfPriors, sizeof(double));
  colSize = DIM[0];
  planeSize = DIM[0]*DIM[1];
  blockSize = DIM[0]*DIM[1]*DIM[2];
  maxRow = DIM[0]-1;
  maxCol = DIM[1]-1;
  maxPlane = DIM[2]-1;

  /* Now loop over all voxels */
  ind = 0;
  for (plane=0; plane<DIM[2]; plane++)
    {
    for (col=0; col<DIM[1]; col++)
      {
      for (row=0; row<DIM[0]; row++)
        {
        if (playing[ind])
          {
          /* Set neighbour counts to 0 */
          for (prior=0; prior<nrOfPriors; prior++)
            {
            u[prior] = 0;
            v[prior] = 0;
            }
          
          /* Look at west neighbour */
          if (col!=0)
            {
            testind = ind - colSize;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                u[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }
          
          /* Look at east neighbour */
          if (col!=maxCol)
            {
            testind = ind + colSize;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                u[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }
          
          /* Look at north neighbour */
          if (row!=0)
            {
            testind = ind - 1;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                u[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }
          
          /* Look at south neighbour */
          if (row!=maxRow)
            {
            testind = ind + 1;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                u[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }
          
          /* Look at top neighbour */
          if (plane!=0)
            {
            testind = ind - planeSize;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                v[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }
          
          /* Look at bottom neighbour */
          if (plane!=maxPlane)
            {
            testind = ind + planeSize;
            if (playing[testind])
              {
              for (class=0; class<nrOfClasses; class++)
                {
                v[lkp[class]] += classification[testind];
                testind += blockSize;
                }
              }
            }



          /* With the calculated neighbour counts, calculate prior */
          for (prior=0; prior<nrOfPriors; prior++)
            {
            voxelprior[prior] = 0;
            }
          
          Gind=0;
          for (neighbourprior=0; neighbourprior<nrOfPriors; neighbourprior++)
            {
            for (prior=0; prior<nrOfPriors; prior++)
              {
              voxelprior[prior] += G[Gind+prior] * 
                                   u[neighbourprior] + 
                                   H[Gind+prior] * 
                                   v[neighbourprior];
              }
            Gind += nrOfPriors;
            }


          sumvoxelprior = 0;
          testind = ind;
          for (prior=0; prior<nrOfPriors; prior++)
            {
            voxelprior[prior] = exp(-voxelprior[prior]/255) * atlas[testind];
            sumvoxelprior += voxelprior[prior];
            testind += blockSize;
            }
          sumvoxelprior /= 255;
          sumvoxelprior += 1E-15;
          testind = ind;
          for (class=0; class<nrOfClasses; class++)
            {
            MRFprior[testind] = (unsigned char) (voxelprior[lkp[class]]
                                        /sumvoxelprior 
                                        /((double) nrOfClassesPerPrior[lkp[class]]) + .5);
            testind += blockSize;
            }
          
 
          }
        ind++;
        }
      }
    }
    

  return;
  }

  
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray* prhs[])

  { 
  /* MRFprior = ems_getMRFprior(classification, playing, atlas, lkp, G, H) */
  
  unsigned char *classification, *playing, *atlas; 
  unsigned int DIM[3], *lkp; 
  unsigned int nrOfClasses, nrOfPriors;
  double *G, *H, *tmplkp; 
  unsigned char *MRFprior;
  
  int number_of_dims, i; 
  const int  *dim_array;


  /* Check if the user wants to see how to use the software */
  if ((nlhs==0) & (nrhs==0))
    {
    mexPrintf("\nUsage:\n\n");
    mexPrintf("MRFprior = ems_getMRFprior(classification, playing, atlas, lkp, G, H)\n\n");
    return;
    }


  /* Check for proper input and output arguments */
  if (nrhs != 6) { 
  mexErrMsgTxt("Six input arguments required."); 
  } else if (nlhs > 1) {
  mexErrMsgTxt("Too many output arguments."); 
    } 
 

  number_of_dims = mxGetNumberOfDimensions(prhs[0]);
  if ((number_of_dims != 4) || !mxIsUint8(prhs[0])) {
  mexErrMsgTxt("Classification must be a 4-dimensional uint8 matrix.");   
  }


  dim_array = mxGetDimensions(prhs[0]);
  DIM[0] = dim_array[0];
  DIM[1] = dim_array[1];
  DIM[2] = dim_array[2];
  nrOfClasses = dim_array[3];

  if ((mxGetNumberOfDimensions(prhs[1]) != 3) || !mxIsUint8(prhs[1])) {
  mexErrMsgTxt("Playing must be a 3-dimensional uint8 matrix.");   
  }

  for (i=0; i<3; i++) {
  if (DIM[i] != mxGetDimensions(prhs[1])[i]) {
  mexErrMsgTxt("Size of playing does not match size of classification.");   
  }
  }

  if ((mxGetNumberOfDimensions(prhs[2]) != 4) || !mxIsUint8(prhs[2])) {
  mexErrMsgTxt("Atlas must be a 4-dimensional uint8 matrix.");   
  }

  for (i=0; i<3; i++) {
  if (DIM[i] != mxGetDimensions(prhs[2])[i]) {
  mexErrMsgTxt("Size of atlas does not match size of classification.");   
  }
  }

  nrOfPriors = mxGetDimensions(prhs[2])[3];

  if ((mxGetNumberOfDimensions(prhs[3]) != 2) || (mxGetM(prhs[3]) != 1) 
    || (mxGetN(prhs[3]) != nrOfClasses)) {
  mexErrMsgTxt("lkp must be [1xm] vector that matches size of classification.");   
  }
  tmplkp = mxGetData(prhs[3]);
  lkp = mxCalloc(nrOfClasses, sizeof(unsigned int));
  for (i=0; i<nrOfClasses; i++) {
  lkp[i] = (unsigned int) (tmplkp[i] + .5) - 1;
  if (lkp[i]>nrOfPriors-1)
    {
    mexErrMsgTxt("lkp contains element(s) that are too high for atlas.");
    }
  }


  if ((mxGetM(prhs[4])!=nrOfPriors) || (mxGetN(prhs[4])!=nrOfPriors)) {
  mexErrMsgTxt("Matrix G does not match size of classification.");   
  }

  if ((mxGetM(prhs[5])!=nrOfPriors) || (mxGetN(prhs[5])!=nrOfPriors)) {
  mexErrMsgTxt("Matrix H does not match size of classification.");   
  }
  

  /* Get input argmuments */
  classification = mxGetData(prhs[0]);
  playing = mxGetData(prhs[1]);
  atlas = mxGetData(prhs[2]);
  G = mxGetPr(prhs[4]);
  H = mxGetPr(prhs[5]);


  /* Create a matrix for the return argument */ 
  plhs[0] = mxCreateNumericArray(4, dim_array, mxUINT8_CLASS, mxREAL);
  MRFprior = mxGetData(plhs[0]);

  
  /* Do the actual computations in a subroutine */
  getMRFprior(classification, playing, DIM, nrOfClasses, atlas, nrOfPriors,
              lkp, G, H, MRFprior);
  
  return;
  
  }

