/*
 * sjlt.c - Sparse Johnson-Lindenstrauss Transform
 *
 * Creates a random sparse Johnson-Lindenstrauss projection matrix. 
 * The columns are independent and each column has exactly s non-zero
 * entries. All non-zero entries are independent Rademacher random
 * variables. Details can be found in [1]. 
 *
 * The calling syntax is:
 *
 *		projection = sjlt(rows, columns, sparsity)
 *
 * This is a MEX file for MATLAB.
 * 
 * Depending on your compiler, you can compile the function using
 * one of the following calls:
 * $ mex CXXFLAGS='$CXXFLAGS -std=c++0x' COPTIMFLAGS='-O3 -DNDEBUG'  -largeArrayDims sjlt.cpp
 * or
 * $ mex CXXFLAGS='$CXXFLAGS -std=c++11' COPTIMFLAGS='-O3 -DNDEBUG'  -largeArrayDims sjlt.cpp
 *
 * Author: Tobias Pohlen <tobias.pohlen@rwth-aachen.de>
 *
 * References:
 * [1] Jean Bourgain, Sjoerd Dirksen, and Jelani Nelson. "Toward a Unified 
 *     Theory of Sparse Dimensionality Reduction in Euclidean Space", 
 *     Symposium on Theory of Computing, 2015. 
 */

#include "mex.h"
#include <random>

std::random_device rd;
std::mt19937 g(rd());

// We use this in order to generate rademacher random variables
std::uniform_int_distribution<int> rademacherDist(0, 1);
inline int rademacher()
{
    return 2*rademacherDist(g) - 1;
}

/* Tries to extract an integer from arg */
mwSize getIntegerScalar(const mxArray* arg)
{
    if (mxGetNumberOfElements(arg) == 1)
    {
        return mxGetScalar(arg);
    }
    else
    {
        mexErrMsgTxt("Integer scalar is not of size == [1 1].\n");
    }
}

/* Returns an integer from arg or 0 if the integer is negative */
mwSize getNonNegativeIntegerScalar(const mxArray* arg)
{
    int res = getIntegerScalar(arg);
    if (res < 0)
    {
        return 0;
    }
    else
    {
        return res;
    }
}

/* Shuffles the array randomly */
void shuffle(
    mwSize* array, 
    mwSize size, 
    std::uniform_int_distribution<mwSize> & indexDistribution)
{
    for (mwSize i = 0; i < size; i++)
    {
        std::swap(array[i], array[indexDistribution(g)]);
    }
}

/* Creates a sparse Johnson Lindenstrauss Transform of size numRows x numCols 
 * of sparsity. 
 */
void createSJLT(
    mwSize sparsity, 
    mwSize numRows, 
    mwSize numCols, 
    double *entries,
    mwSize* rowIndices, 
    mwSize* colIndices)
{
    // Create an array of row indices to shuffle. We use this in order
    // to draw random rows without replacement
    std::uniform_int_distribution<mwSize> rowDist(0, numRows-1);
    mwSize* rowCache = (mwSize*) malloc(numRows*sizeof(mwSize));
    for (mwSize i = 0; i < numRows; i++)
    {
        rowCache[i] = i;
    }
    
    // Fill the column indices and the entries (remember that the entries are
    // just independent rademacher random variables)
    mwSize colOffset = 0;
    for (mwSize c = 0; c < numCols; c++)
    {
        // Shuffle the row indices
        shuffle(rowCache, sparsity, rowDist);

        for (mwSize s = 0; s < sparsity; s++)
        {
            entries[colOffset+s] = rademacher();
            rowIndices[colOffset+s] = rowCache[s];
        }

        colIndices[c] = c*sparsity;

        colOffset += sparsity;
    }

    colIndices[numCols] = numCols*sparsity;

    free(rowCache);
}

/*
 * This is the function called by MATLAB. 
 */
void mexFunction(
    int numLeftHandSide, 
    mxArray *pointerLeftHandSide[],
    int numRightHandSide, 
    const mxArray *pointerRightHandSide[])
{
    // Inputs:
    // 1. number of rows
    // 2. number of columns
    // 3. sparsity (number of non-zeros per column)
    if(numRightHandSide != 3)
    {
        mexErrMsgIdAndTxt(
            "arrayProduct:numRightHandSide",
            "Three inputs required.");
    }
    
    // Outputs:
    // 1. SJLT matrix
    if (numLeftHandSide != 1)
    {
        mexErrMsgIdAndTxt(
            "arrayProduct:numLeftHandSide",
            "One output required.");
    }
    
    // Read the inputs
    int numRows = getNonNegativeIntegerScalar(pointerRightHandSide[0]);
    int numCols = getNonNegativeIntegerScalar(pointerRightHandSide[1]);
    int sparsity = getNonNegativeIntegerScalar(pointerRightHandSide[2]);
    
    // The sparsity cannot be higher than the number of rows
    if (sparsity > numRows)
    {
        sparsity = numRows;
    }
    
    // Create the outputs
    pointerLeftHandSide[0] = mxCreateSparse(numRows,numCols,numCols*sparsity,mxREAL);
    
    // Create the transformation
    createSJLT(
        sparsity, 
        numRows,
        numCols, 
        mxGetPr(pointerLeftHandSide[0]), 
        mxGetIr(pointerLeftHandSide[0]), 
        mxGetJc(pointerLeftHandSide[0]));
}