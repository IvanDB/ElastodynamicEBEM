#include "commonRoutines.cuh"

__global__ void kernelKboundary(double* __restrict__ matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double const4PiDeltaT,
                          const double* __restrict__ stdIntW, const double* __restrict__ stdIntNx, const double* __restrict__ stdIntNy, const double* __restrict__ stdIntNz, const int numPointExt,
                          const double* __restrict__ stdExtW, const double* __restrict__ stdExtNx, const double* __restrict__ stdExtNy, const double* __restrict__ stdExtNz,
                          const double* __restrict__ vertsT, const double* __restrict__ areeT, const double* __restrict__ normT, const int offsetZ, const int numBlocks, const double maxLen) 
{
    //Skip diagonal blocks
    if(blockIdx.x == blockIdx.y)
       return;

    //Time instant check
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Spatial condition check
    if(isKBlockNull_costant(vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
       return;

    //Shared memory initialization
    extern __shared__ double matrixSubBlock[][3][3];
    const size_t sharedBaseInd = threadIdx.x;
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Source triangle verteces
    double extVerts[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            extVerts[i][j] = vertsT[9*blockIdx.x + 3*j + i];  

    //Field triangle verteces
    double intVerts[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            intVerts[i][j] = vertsT[9*blockIdx.y + 3*j + i];

    //Field triangle vector components
    double intEdgeVectors[3][3];
    for(size_t j = 0; j < 3; ++j)
    {
        intEdgeVectors[0][j] = intVerts[1][j] - intVerts[0][j];
        intEdgeVectors[1][j] = intVerts[2][j] - intVerts[1][j];
        intEdgeVectors[2][j] = intVerts[0][j] - intVerts[2][j];
    }

    //Field triangle edge lengths
    const double intEdgeLengths[] = {norm2(intEdgeVectors[0]), 
                                     norm2(intEdgeVectors[1]), 
                                     norm2(intEdgeVectors[2])};
    
    //Point factors
    const double startFactor = (double) threadIdx.x / blockDim.x;
    const double endFactor = (double) (threadIdx.x + 1) / blockDim.x;

    //Eta costants
    const double etaCoeffs[4] = {-1, 3, -3, 1};
    const double etaValues[4] = {-2, -1, 0, 1};

    //Loop over source nodes
    for(size_t l = 0; l < numPointExt; ++l)
    {
        //Source weight
        const double extWeight = stdIntW[l] * areeT[blockIdx.x];

        //Source node
        const double stdNodeE[3] = {stdIntNx[l], stdIntNy[l], stdIntNz[l]};

        const double extNode[] = {stdNodeE[0] * extVerts[0][0] + stdNodeE[1] * extVerts[1][0] + stdNodeE[2] * extVerts[2][0],
                                  stdNodeE[0] * extVerts[0][1] + stdNodeE[1] * extVerts[1][1] + stdNodeE[2] * extVerts[2][1],
                                  stdNodeE[0] * extVerts[0][2] + stdNodeE[1] * extVerts[1][2] + stdNodeE[2] * extVerts[2][2]};
            
        //Loop over the field triangle edges
        for(size_t m = 0; m < 3; ++m)
        {
            const double segmentLength = intEdgeLengths[m] / blockDim.x;

            //Space vectors
            double pointStart[3];
            for(size_t i = 0; i < 3; ++i)
                pointStart[i] = extNode[i] - (intVerts[m][i] + startFactor * intEdgeVectors[m][i]);
                    
            double pointEnd[3];
            for(size_t i = 0; i < 3; ++i)
                pointEnd[i] = extNode[i] - (intVerts[m][i] + endFactor * intEdgeVectors[m][i]);

            const double normStart = norm2(pointStart); 
            const double normEnd = norm2(pointEnd); 
            
            //Edge tangent versor
            double tangentVersor[3];
            for(size_t i = 0; i < 3; ++i)
                tangentVersor[i] = intEdgeVectors[m][i] / intEdgeLengths[m];
    
            //Loop over eta factor
            for(size_t k = 0; k < 4; ++k)
            {
                const double timeInstant = deltaT * (offsetZ + double(blockIdx.z) + etaValues[k]);

                if(timeInstant <= 0)
                   continue;

                double tempValues[3][3] = {0.};
                
                //Add KR sub-kernel at start and end point
                nuKR_costant(tempValues, tangentVersor, pointStart, normStart, timeInstant, velP, velS, lambda, mu, rho);
                nuKR_costant(tempValues, tangentVersor, pointEnd, normEnd, timeInstant, velP, velS, lambda, mu, rho);

                //Add to shared memory
                for(size_t i = 0; i < 3; ++i)
                    for(size_t j = 0; j < 3; ++j)
                        matrixSubBlock[sharedBaseInd][i][j] += extWeight * etaCoeffs[k] * (tempValues[i][j] * segmentLength / 2); 
            }
        }
    }

    __syncthreads();

    //Loop to reduce in x dimension (x size is assumed to be a power of 2)
    size_t xDim = blockDim.x / 2;
    while(xDim > 0)
    {
        size_t sharedOffInd = threadIdx.x + xDim;
        if(threadIdx.x < xDim)
            for(size_t i = 0; i < 3; ++i)
                for(size_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    //Save in global memory
    if(threadIdx.x == 0)
        for(size_t i = 0; i < 3; ++i)
            for(size_t j = 0; j < 3; ++j)
            {
                const size_t ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                //matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiDeltaT;
            }
}



