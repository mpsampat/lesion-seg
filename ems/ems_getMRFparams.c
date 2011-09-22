/* ems_getMRFparams.c    Koen Van Leemput - August 17, 2001 */

#include <math.h>
#include "mex.h"

#define THRES1OLD 1E-4
#define THRES2OLD 1E-10
#define THRES1 1E-7
#define THRES2 1E-13


static void buildMRFhistogram(unsigned char *classification, 
                              unsigned int DIM[3], 
                              unsigned int nrOfClasses, 
                              unsigned int *lkp,
                              unsigned int *sampleInd, 
                              unsigned int nrOfSamples, 
                              double ***histogramPtr, 
                              unsigned char ***uListPtr,
                              unsigned char ***vListPtr,
                              unsigned int *nrOfusPtr, 
                              unsigned int *nrOfvsPtr,
                              unsigned int *nrOfPriorsPtr)
  {
  unsigned int *powersOf10, nrOfus, nrOfvs, *uReverseList, *vReverseList,
    nrOfNeighbourhoods, uReverseListNr, vReverseListNr, row, col, plane,
    nrOfPlanes, nrOfRows, nrOfCols, uNr, vNr, wInd, eInd, sInd, nInd, tInd, 
    bInd, ind, cLabel, wLabel, eLabel, nLabel, sLabel, tLabel, bLabel, colSize,
    planeSize, blockSize, neighbourhoodNr, neighbourhoodNrOffset, sample,
    offset, class, lastMessagePlane, nrOfPriors, *knownPrior, prior;
  mxArray *lhs[10], *rhs[10];
  double rest, *cClassification, *wClassification, *eClassification, 
    *sClassification, *nClassification, *tClassification, *bClassification,
    tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

  double **histogram;
  unsigned char **uList, **vList;


  mexPrintf("    Building neighborhood configuration histogram:      ");


  /* How many different priors are there? */
  nrOfPriors = 0;
  knownPrior = mxCalloc(nrOfClasses, sizeof(unsigned int));
  for (class=0; class<nrOfClasses; class++)
    {
    if (!knownPrior[class])
      {
      for (ind=0; ind<nrOfClasses; ind++)
        {
        if (lkp[ind]==lkp[class])
          {
          knownPrior[ind] = 1;
          }
        }
      nrOfPriors++;
      }
    }
  *nrOfPriorsPtr = nrOfPriors;



  /* Let's assign a follow-up number to each neighbourhood u 
     (4 neighbours in plane) */
  powersOf10 = mxCalloc(nrOfPriors, sizeof(unsigned int));
  for (prior=0; prior<nrOfPriors; prior++)
    {
    powersOf10[prior] = pow(10,nrOfPriors-1-prior);
    }
  

  rhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(rhs[0]) = nrOfPriors;
  rhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(rhs[1]) = 4;
  mexCallMATLAB(1,lhs,2,rhs,"ems_binBalls");

  nrOfus = mxGetM(lhs[0]);
  *nrOfusPtr = nrOfus;
  uList = mxCalloc(nrOfPriors, sizeof(unsigned char*));
  for (prior=0; prior<nrOfPriors; prior++)
    {
    uList[prior] = mxCalloc(nrOfus, sizeof(unsigned char));
    for (uNr=0; uNr<nrOfus; uNr++)
      {
      uList[prior][uNr] = (unsigned char) mxGetPr(lhs[0])[uNr+prior*nrOfus];
      }
    }
  *uListPtr = uList;

  uReverseList = mxCalloc(4*pow(10,nrOfPriors-1), sizeof(unsigned int));
  for (uNr=0; uNr<nrOfus; uNr++)
    {
    uReverseListNr = 0;
    for (prior=0; prior<nrOfPriors; prior++)
      {
      uReverseListNr += uList[prior][uNr]*powersOf10[prior];
      }
    uReverseList[uReverseListNr-1] = uNr;
    }

  /* Let's assign a follow-up number to each neighbourhood v 
     (2 neighbours out of plane) */
  rhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(rhs[0]) = nrOfPriors;
  rhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(rhs[1]) = 2;
  mexCallMATLAB(1,lhs,2,rhs,"ems_binBalls");

  nrOfvs = mxGetM(lhs[0]);
  *nrOfvsPtr = nrOfvs;
  vList = mxCalloc(nrOfPriors, sizeof(unsigned char*));
  for (prior=0; prior<nrOfPriors; prior++)
    {
    vList[prior] = mxCalloc(nrOfvs, sizeof(unsigned char));
    for (vNr=0; vNr<nrOfvs; vNr++)
      {
      vList[prior][vNr] = (unsigned char) mxGetPr(lhs[0])[vNr+prior*nrOfvs];
      }
    }
  *vListPtr = vList;
  
  vReverseList = mxCalloc(2*pow(10,nrOfPriors-1), sizeof(unsigned int));
  for (vNr=0; vNr<nrOfvs; vNr++)
    {
    vReverseListNr = 0;
    for (prior=0; prior<nrOfPriors; prior++)
      {
      vReverseListNr += vList[prior][vNr]*powersOf10[prior];
      }
    vReverseList[vReverseListNr-1] = vNr;
    }


  /* Initialize histogram neighbourhood vs. label to zero */
  nrOfNeighbourhoods = nrOfus * nrOfvs;
  histogram = mxCalloc(nrOfPriors, sizeof(double*));
  for (prior=0; prior<nrOfPriors; prior++)
    {
    histogram[prior] = mxCalloc(nrOfNeighbourhoods, sizeof(double));
    }
  *histogramPtr = histogram;

  /* Calculate histogram */
  colSize = DIM[0];
  planeSize = DIM[0]*DIM[1];
  blockSize = DIM[0]*DIM[1]*DIM[2];
  nrOfPlanes = DIM[2];
  nrOfCols = DIM[1];
  nrOfRows = DIM[0];
  wClassification = mxCalloc(nrOfPriors, sizeof(double));
  eClassification = mxCalloc(nrOfPriors, sizeof(double));
  sClassification = mxCalloc(nrOfPriors, sizeof(double));
  nClassification = mxCalloc(nrOfPriors, sizeof(double));
  tClassification = mxCalloc(nrOfPriors, sizeof(double));
  bClassification = mxCalloc(nrOfPriors, sizeof(double));
  cClassification = mxCalloc(nrOfPriors, sizeof(double));

  lastMessagePlane = 0;
  for (sample=0; sample<nrOfSamples; sample++)
    {
    ind = sampleInd[sample];

    /* Determine 3D position of this sample */
    rest = ind;
    plane = (unsigned int) (rest/((double) planeSize));
    rest -= plane*planeSize;
    col = (unsigned int) (rest/((double) colSize));
    rest -= col*colSize;
    row = rest;


    if ((plane!=0) & (plane!=nrOfPlanes-1) & (col!=0) & (col!=nrOfCols-1)
      & (row!=0) & (row!=nrOfRows-1))
      {

      /* What is classification of neigbours? */
      for (prior=0; prior<nrOfPriors; prior++)
        {
        wClassification[prior] = 0;
        eClassification[prior] = 0;
        sClassification[prior] = 0;
        nClassification[prior] = 0;
        tClassification[prior] = 0;
        bClassification[prior] = 0;
        cClassification[prior] = 0;
        }


      wInd = ind - colSize;
      eInd = ind + colSize;
      sInd = ind + 1;
      nInd = ind - 1;
      tInd = ind - planeSize;
      bInd = ind + planeSize;
      offset = 0;
      for (class=0; class<nrOfClasses; class++)
        {
        wClassification[lkp[class]] += classification[wInd + offset];
        eClassification[lkp[class]] += classification[eInd + offset];
        sClassification[lkp[class]] += classification[sInd + offset];
        nClassification[lkp[class]] += classification[nInd + offset];
        tClassification[lkp[class]] += classification[tInd + offset];
        bClassification[lkp[class]] += classification[bInd + offset];
        cClassification[lkp[class]] += classification[ind + offset];
        offset += blockSize;
        }



      /* Calculate probability for every possible MRF configuration 
       and add this probability to the MRF configuration histogram */
      for (tLabel=0; tLabel<nrOfPriors; tLabel++)
        {
        tmp1 = tClassification[tLabel];
        for (bLabel=0; bLabel<nrOfPriors; bLabel++)
          {
          vReverseListNr = powersOf10[tLabel] + powersOf10[bLabel];
          vNr = vReverseList[vReverseListNr-1];
          neighbourhoodNrOffset = vNr*nrOfus;
          tmp2 = bClassification[bLabel] * tmp1;
          for (wLabel=0; wLabel<nrOfPriors; wLabel++)
            {
            tmp3 = wClassification[wLabel] * tmp2;
            for (eLabel=0; eLabel<nrOfPriors; eLabel++)
              {
              tmp4 = eClassification[eLabel] * tmp3;
              for (sLabel=0; sLabel<nrOfPriors; sLabel++)
                {
                tmp5 = sClassification[sLabel] * tmp4;
                for (nLabel=0; nLabel<nrOfPriors; nLabel++)
                  {
                  tmp6 = nClassification[nLabel] * tmp5;
                  uReverseListNr = powersOf10[wLabel] + powersOf10[eLabel] + 
                                   powersOf10[sLabel] + powersOf10[nLabel];
                  uNr = uReverseList[uReverseListNr-1];
                  neighbourhoodNr = neighbourhoodNrOffset + uNr;
                  for (cLabel=0; cLabel<nrOfPriors; cLabel++)
                    {
                    histogram[cLabel][neighbourhoodNr] 
                      += tmp6 * cClassification[cLabel];
                    }
                  }
                }
              }
            }
          }
        }



      } /* End test if this sample falls well within ROI */

    if (plane!=lastMessagePlane)
      {
      lastMessagePlane = plane;
      for (prior=0; prior<5; prior++)
        {
        mexPrintf("\b");
        }
      mexPrintf("%-3d %%", 
                (unsigned int) (((double) (plane+1))/((double) DIM[2])*100));
#if 1
      fflush(stdout);
#endif
      }


    } /* End loop over samples */


    for (prior=0; prior<5; prior++)
      {
      mexPrintf("\b");
      }
    mexPrintf("%-3d %%", 100);

  mexPrintf("\n");

  return;
  }






static void MRFhistogram2params(double **histogram, unsigned char **uList, 
                                unsigned char **vList, unsigned int nrOfus, 
                                unsigned int nrOfvs, unsigned int nrOfPriors, 
                                unsigned int nrOfSamples, unsigned int constraintType,
                                double **Gptr, double **Hptr)
  {
  
  unsigned int nrOfNeigbourhoods, prior, neighbourhoodNr, nrOfZeros,
    nrOfEquations, equation, j, k, l, nrOfMRFparams, nrOfParamsInG,
    equationNr, uNr, vNr, nrOfNeighbourhoods, row, col, ind, indG, indH, indG2, indH2;
  double *lhs, *rhs, toAdd, threshold1, threshold2, *G, *H, sumG, sumH;
  mxArray *lhsMatrix, *rhsMatrix, *plhs[10], *prhs[10];
  
  nrOfNeighbourhoods = nrOfus * nrOfvs;

  /* Count number of zero entries in histogram */
  nrOfZeros=0;
  for (prior=0; prior<nrOfPriors; prior++)
    {
    for (neighbourhoodNr=0; neighbourhoodNr<nrOfNeighbourhoods;
         neighbourhoodNr++)
      {
      if (histogram[prior][neighbourhoodNr]==0)
        {
        nrOfZeros++;
        }
      }
    }  

#if 0
  mexPrintf("    nr of zero entries in histogram = %d\n", nrOfZeros);
#endif


  /* Set up set of linear equations */
  nrOfEquations = nrOfNeighbourhoods * 
                  (unsigned int) 
                  ((double)(nrOfPriors-1) * ((double) nrOfPriors) / 2);
  nrOfParamsInG = pow(nrOfPriors,2);
  nrOfMRFparams = 2*pow(nrOfPriors,2);


  lhsMatrix = mxCreateDoubleMatrix(nrOfEquations, nrOfMRFparams, mxREAL);
  rhsMatrix = mxCreateDoubleMatrix(nrOfEquations, 1, mxREAL);
  lhs = mxGetPr(lhsMatrix);
  rhs = mxGetPr(rhsMatrix);

  threshold1 = pow(255,7) * THRES1 * ((double) nrOfSamples);
  threshold2 = pow(255,7) * THRES2 * ((double) nrOfSamples); 


  equationNr=0;
  for (vNr=0; vNr<nrOfvs; vNr++)
    {
    for (uNr=0; uNr<nrOfus; uNr++)
      {
      neighbourhoodNr = uNr + nrOfus*vNr;
      for (j=0; j<nrOfPriors; j++)
        {
        for (k=j+1; k<nrOfPriors; k++)
          {
          if ((histogram[j][neighbourhoodNr]>threshold1) | 
              (histogram[k][neighbourhoodNr]>threshold1))
            {

            /* What is rhs of the current equation */
            if (histogram[j][neighbourhoodNr]<threshold2)
              {
              rhs[equationNr] = log(threshold2 /
                                   (histogram[k][neighbourhoodNr]));
              }
            else if (histogram[k][neighbourhoodNr]<threshold2)
              {
              rhs[equationNr] = log((histogram[j][neighbourhoodNr]) /
                                    threshold2);
              }
            else
              {
              rhs[equationNr] = log((histogram[j][neighbourhoodNr]) /
                                   (histogram[k][neighbourhoodNr]));
              }

            /* What is lhs of the current equation */
            for (l=0; l<nrOfPriors; l++)
              {
              lhs[equationNr + 
                 (j*nrOfPriors+l)*nrOfEquations] = -uList[l][uNr];
              lhs[equationNr + 
                 (k*nrOfPriors+l)*nrOfEquations] = uList[l][uNr];
              lhs[equationNr + 
                 (j*nrOfPriors+l+nrOfParamsInG)*nrOfEquations] = 
                -vList[l][vNr];
              lhs[equationNr + 
                 (k*nrOfPriors+l+nrOfParamsInG)*nrOfEquations] = 
                vList[l][vNr];
              }
            equationNr++;
            } /* End test if histogram can be trusted for this specific MRF 
               configuration */
          }
        } /* End loop over all pairs of central voxel labels */
      } 
    } /* End loop over all possible neigbourhoods */





  /* Dependent on value of constraintType, add constraints that some elements of G and H must be 
     the same */
  if (constraintType==1)
    {
    /* adjacent classes have the same prior probability in each other's neighborhood */
    for (col=0; col<nrOfPriors; col++)
      {
      indG = (col + col*nrOfPriors) * nrOfEquations;
      indH = indG + nrOfParamsInG*nrOfEquations;

      if (col==0)
        {
        for (equationNr=0; equationNr<nrOfEquations; equationNr++)
          {
          lhs[indG] += lhs[indG + nrOfPriors*nrOfEquations];
          lhs[indG] /= 2;
          lhs[indG + nrOfPriors*nrOfEquations] = lhs[indG];
          indG++;

          lhs[indH] += lhs[indH + nrOfPriors*nrOfEquations];
          lhs[indH] /= 2;
          lhs[indH + nrOfPriors*nrOfEquations] = lhs[indH];
          indH++;
          }
        }
      else if (col==nrOfPriors-1)
        {
        for (equationNr=0; equationNr<nrOfEquations; equationNr++)
          {
          lhs[indG] += lhs[indG - nrOfPriors*nrOfEquations];
          lhs[indG] /= 2;
          lhs[indG - nrOfPriors*nrOfEquations] = lhs[indG];
          indG++;

          lhs[indH] += lhs[indH - nrOfPriors*nrOfEquations];
          lhs[indH] /= 2;
          lhs[indH - nrOfPriors*nrOfEquations] = lhs[indH];
          indH++;
          }
        }
      else
        {
        for (equationNr=0; equationNr<nrOfEquations; equationNr++)
          {
          lhs[indG] += lhs[indG - nrOfPriors*nrOfEquations] +
                      lhs[indG + nrOfPriors*nrOfEquations];
          lhs[indG] /= 3;
          lhs[indG - nrOfPriors*nrOfEquations] = lhs[indG];
          lhs[indG + nrOfPriors*nrOfEquations] = lhs[indG];
          indG++;

          lhs[indH] += lhs[indH - nrOfPriors*nrOfEquations] +
                      lhs[indH + nrOfPriors*nrOfEquations];
          lhs[indH] /= 3;
          lhs[indH - nrOfPriors*nrOfEquations] = lhs[indH];
          lhs[indH + nrOfPriors*nrOfEquations] = lhs[indH];
          indH++;
          }
        }

      }
    }
  else if (constraintType==2)
    {
    /* Transition matrices are symmetric */
    for (row=0; row<nrOfPriors; row++)
      {
      for (col=row+1; col<nrOfPriors; col++)
        {
        indG = (col + row*nrOfPriors) * nrOfEquations;
        indH = indG + nrOfParamsInG*nrOfEquations;

        indG2 = (row + col*nrOfPriors) * nrOfEquations;
        indH2 = indG2 + nrOfParamsInG*nrOfEquations;
        
        for (equationNr=0; equationNr<nrOfEquations; equationNr++)
          {
          lhs[indG] += lhs[indG2];
          lhs[indG] /= 2;
          lhs[indG2] = lhs[indG];
          indG++;
          indG2++;
          
          lhs[indH] += lhs[indH2];
          lhs[indH] /= 2;
          lhs[indH2] = lhs[indH];
          indH++;
          indH2++;
          }

        }
      
      indG = (row + row*nrOfPriors) * nrOfEquations;
      indH = indG + nrOfParamsInG*nrOfEquations;
      for (equationNr=0; equationNr<nrOfEquations; equationNr++)
        {
        lhs[indG] = 0;
        indG++;
        
        lhs[indH] = 0;
        indH++;
        }
      }
    
    }
  else if (constraintType==3)
    {
    /* All non-diagonal elements are the same */
    /* mexPrintf("Testerdetest!\n"); */


    /* First, let's set diagonal elements to 0 */
    for (col=0; col<nrOfPriors; col++)
      {
      indG = (col + col*nrOfPriors) * nrOfEquations;
      indH = indG + nrOfParamsInG*nrOfEquations;

      for (equationNr=0; equationNr<nrOfEquations; equationNr++)
        {
        lhs[indG] = 0;
        indG++;
        
        lhs[indH] = 0;
        indH++;
        }
      
      }

    /* Now add contraint that all non-zero elements must be equal */
    for (equationNr=0; equationNr<nrOfEquations; equationNr++)
      {
      sumG = 0;
      sumH = 0;
      for (row=0; row<nrOfPriors; row++)
        {
        for (col=0; col<nrOfPriors; col++)
          {
          if ((row>col) || (row<col))
            {
            indG = (col+row*nrOfPriors) * nrOfEquations + equationNr;
            indH = indG + nrOfParamsInG*nrOfEquations;
            sumG += lhs[indG];
            sumH += lhs[indH];
            }
          }
        }
      sumG /= (double) ((nrOfPriors)*(nrOfPriors-1));
      sumH /= (double) ((nrOfPriors)*(nrOfPriors-1));

      for (row=0; row<nrOfPriors; row++)
        {
        for (col=0; col<nrOfPriors; col++)
          {
          if ((row>col) || (row<col))
            {
            indG = (col+row*nrOfPriors) * nrOfEquations + equationNr;
            indH = indG + nrOfParamsInG*nrOfEquations;
            lhs[indG] = sumG;
            lhs[indH] = sumH;
            }
          }
        }
      }


    }






  /* Solve least-squares fit */
#if 0
  mxSetName(lhsMatrix, "lhs");
  mexPutArray(lhsMatrix, "base");
  mxSetName(rhsMatrix, "rhs");
  mexPutArray(rhsMatrix, "base");
#endif

  prhs[0] = lhsMatrix;
  mexCallMATLAB(1,plhs,1,prhs,"pinv");
  prhs[0] = plhs[0];
  prhs[1] = rhsMatrix;
  mexCallMATLAB(1,plhs,2,prhs,"*");

  G = mxCalloc(nrOfPriors*nrOfPriors, sizeof(double));
  H = mxCalloc(nrOfPriors*nrOfPriors, sizeof(double));
  *Gptr = G;
  *Hptr = H;

  ind=0;
  for (row=0; row<nrOfPriors; row++)
    {
    for (col=0; col<nrOfPriors; col++)
      {
      G[row + col*nrOfPriors] = (mxGetPr(plhs[0]))[ind];
      H[row + col*nrOfPriors] = (mxGetPr(plhs[0]))[ind+nrOfParamsInG];
      ind++;
      }
    }


  /* Add constant to every column of G and H so that diagonal elements are
     zero */
  for (col=0; col<nrOfPriors; col++)
    {
    toAdd = -G[col + col*nrOfPriors];
    for (row=0; row<nrOfPriors; row++)
      {
      G[row + col*nrOfPriors] += toAdd;
      }
    toAdd = -H[col + col*nrOfPriors];
    for (row=0; row<nrOfPriors; row++)
      {
      H[row + col*nrOfPriors] += toAdd;
      }
    }

  return;
  }




static void getMRFparams(unsigned char *classification, unsigned int DIM[3], 
                         unsigned int nrOfClasses, unsigned int *lkp, 
                         unsigned int *sampleInd, 
                         unsigned int nrOfSamples, unsigned int constraintType,
                         mxArray **GmatrixPtr, mxArray **HmatrixPtr)
  {
  double **histogram = NULL, *G = NULL, *H = NULL;
  unsigned char **uList = NULL, **vList = NULL;
  unsigned int nrOfus, nrOfvs, nrOfPriors, ind;


  buildMRFhistogram(classification, DIM, nrOfClasses, lkp, 
                    sampleInd, nrOfSamples, 
                    &histogram, &uList, &vList, &nrOfus, &nrOfvs, &nrOfPriors);

  MRFhistogram2params(histogram, uList, vList, nrOfus, nrOfvs, nrOfPriors, 
                      nrOfSamples, constraintType, &G, &H);

  *GmatrixPtr = mxCreateDoubleMatrix(nrOfPriors, nrOfPriors, mxREAL);
  *HmatrixPtr = mxCreateDoubleMatrix(nrOfPriors, nrOfPriors, mxREAL);
  for (ind=0; ind<nrOfPriors*nrOfPriors; ind++)
    {
    (mxGetPr(*GmatrixPtr))[ind] = G[ind];
    (mxGetPr(*HmatrixPtr))[ind] = H[ind];
    }

  return;
  }



  
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray* prhs[])

  { 
  /* [G, H] = ems_getMRFparams(classification, lkp, sampleInd, constraintType) */
  
  unsigned char *classification; 
  unsigned int DIM[3]; 
  unsigned int nrOfClasses, *sampleInd, nrOfSamples, constraintType, *lkp;
  
  
  int number_of_dims, i; 
  const int  *dim_array;
  unsigned int nrOf3DElements;
  double *tmpsampleInd;

  /* Check if the user wants to see how to use the software */
  if ((nlhs==0) & (nrhs==0))
    {
    mexPrintf("\nUsage:\n\n");
    mexPrintf("[G, H] = ems_getMRFparams(classification, lkp, sampleInd, constraintType)\n\n");
    return;
    }


  /* Check for proper input and output arguments */
  if (nrhs != 4) { 
  mexErrMsgTxt("Four input arguments required."); 
  } else if (nlhs != 2) {
  mexErrMsgTxt("Two output arguments required."); 
    } 
 

  number_of_dims = mxGetNumberOfDimensions(prhs[0]);
  if ((number_of_dims != 4) || !mxIsUint8(prhs[0])) {
  mexErrMsgTxt("classification must be a 4-dimensional uint8 matrix.");   
  }


  dim_array = mxGetDimensions(prhs[0]);
  DIM[0] = dim_array[0];
  DIM[1] = dim_array[1];
  DIM[2] = dim_array[2];
  nrOfClasses = dim_array[3];


  if ((mxGetNumberOfDimensions(prhs[1]) != 2) || (mxGetM(prhs[1]) != 1) 
    || (mxGetN(prhs[1]) != nrOfClasses)) {
  mexErrMsgTxt("lkp must be [1xm] vector that matches size of classification.");   
  }




  if ((mxGetNumberOfDimensions(prhs[2]) != 2) || (mxGetN(prhs[2]) != 1)) 
    {
    mexErrMsgTxt("sampleInd must be a [mx1] vector.");
    }

  if ((mxGetNumberOfDimensions(prhs[3]) != 2) || (mxGetN(prhs[3]) != 1) ||
    (mxGetM(prhs[3]) != 1)) 
    {
    mexErrMsgTxt("constraintType must be a scalar.");
    }

  

  /* Get input argmuments */
  classification = mxGetData(prhs[0]);
  lkp = mxCalloc(nrOfClasses, sizeof(unsigned int));
  for (i=0; i<nrOfClasses; i++) 
    {
    lkp[i] = (unsigned int) ((mxGetPr(prhs[1]))[i] + .5) - 1;
    }
  nrOfSamples = mxGetM(prhs[2]);
  nrOf3DElements = DIM[0]*DIM[1]*DIM[2];
  tmpsampleInd = mxGetData(prhs[2]);
  sampleInd = mxCalloc(nrOfSamples, sizeof(unsigned int));
  for (i=0; i<nrOfSamples; i++)
    {
    sampleInd[i] = (unsigned int) (tmpsampleInd[i]-1); 
    if (sampleInd[i]>nrOf3DElements-1)
      {
      mexErrMsgTxt("sampleInd contains indices that are outside.");
      }
    }
  constraintType = mxGetPr(prhs[3])[0];


  /* Create matrices for the return arguments */ 
  plhs[0] = NULL;
  plhs[1] = NULL;

  
  /* Do the actual computations in a subroutine */
  getMRFparams(classification, DIM, nrOfClasses, lkp, sampleInd, nrOfSamples, 
               constraintType, &plhs[0], &plhs[1]);
  
  return;
  
  }

