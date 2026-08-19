#include "commonRoutines.cuh"

__global__ void kernelK(double* __restrict__ matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double const4PiDeltaT,
                          const double* __restrict__ stdExtW, const double* __restrict__ stdExtNx, const double* __restrict__ stdExtNy, const double* __restrict__ stdExtNz, const uint32_t numExtNodes,
                          const double* __restrict__ stdIntW, const double* __restrict__ stdIntNx, const double* __restrict__ stdIntNy, const double* __restrict__ stdIntNz,
                          const double* __restrict__ vertsT, const double* __restrict__ areeT, const double* __restrict__ normT, const int* __restrict__ indSMmatrix, const double* __restrict__ matCoeff, const double* __restrict__ vetCoeff, 
                          const int offsetZ, const int numBlocks, const double* __restrict__ nodesMesh, const double maxLen)  
{
    //Time instant check
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Spatial condition check
    if(isKBlockNull_linear(nodesMesh, vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
        return;

    //Shared memory initialization
    extern __shared__ double matrixSubBlock[][3][3];
    const uint32_t sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    for(uint32_t i = 0; i < 3; ++i)
        for(uint32_t j = 0; j < 3; ++j)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Source triangle verteces
    double extVerts[3][3];
    for(uint32_t i = 0; i < 3; ++i)
        for(uint32_t j = 0; j < 3; ++j)
            extVerts[i][j] = vertsT[9*blockIdx.x + 3*j + i];   

    //Base index in flag matrix
    const uint32_t baseFlagIdx = blockIdx.y * gridDim.x;

    //Eta costants
    const double etaCoeffs[4] = {-1, 3, -3, 1};
    const double etaValues[4] = {-2, -1, 0, 1};

    //Loop over field triangles
    for(uint32_t m = 0; m < gridDim.x; ++m)
    { 
        //Skip singular cases
        if(m == blockIdx.x)
            continue;

        //Vertex-triangle index
        const uint32_t indSMcurr = indSMmatrix[baseFlagIdx + m];
        if(indSMcurr == 0)
            continue;

        //Internal normal vector
        double intNormVect[3];
        for(uint32_t i = 0; i < 3; ++i)
            intNormVect[i] = normT[3*m + i];

        //Basis matrix and vector coefficient 
        double baseCoeffMatrix[3][3];
        for(uint32_t i = 0; i < 3; ++i)
            for(uint32_t j = 0; j < 3; ++j)
               baseCoeffMatrix[i][j] = matCoeff[9*m + 3*j + i];

        double baseCoeffVector[3];
        for(uint32_t i = 0; i < 3; ++i)
            baseCoeffVector[i] = vetCoeff[3*m + i];

        //Field weight
        const double intWeight = stdIntW[threadIdx.y] * areeT[m];

        //Field verteces
        double intVerts[3][3];
        for(uint32_t i = 0; i < 3; ++i)
            for(uint32_t j = 0; j < 3; ++j)
                intVerts[i][j] = vertsT[9*m + 3*j + i];

        //Field node
        const double stdNodeI[] = {stdIntNx[threadIdx.x*blockDim.y + threadIdx.y],
                                   stdIntNy[threadIdx.x*blockDim.y + threadIdx.y],
                                   stdIntNz[threadIdx.x*blockDim.y + threadIdx.y]};

        const double intNode[3] = {stdNodeI[0] * intVerts[0][0] + stdNodeI[1] * intVerts[1][0] + stdNodeI[2] * intVerts[2][0],
                                   stdNodeI[0] * intVerts[0][1] + stdNodeI[1] * intVerts[1][1] + stdNodeI[2] * intVerts[2][1],
                                   stdNodeI[0] * intVerts[0][2] + stdNodeI[1] * intVerts[1][2] + stdNodeI[2] * intVerts[2][2]};

        //Loop over source nodes
        for(uint32_t l = 0; l < numExtNodes; ++l)
        {
            //Source weight
            const double extWeight = stdExtW[l] * areeT[blockIdx.x];
    
            //Source node 
            const double stdNodeE[] = {stdExtNx[l], stdExtNy[l], stdExtNz[l]};

            const double extNode[] = {stdNodeE[0] * extVerts[0][0] + stdNodeE[1] * extVerts[1][0] + stdNodeE[2] * extVerts[2][0],
                                      stdNodeE[0] * extVerts[0][1] + stdNodeE[1] * extVerts[1][1] + stdNodeE[2] * extVerts[2][1],
                                      stdNodeE[0] * extVerts[0][2] + stdNodeE[1] * extVerts[1][2] + stdNodeE[2] * extVerts[2][2]};

            //Space vector
            const double point[] = {extNode[0] - intNode[0],
                                    extNode[1] - intNode[1],
                                    extNode[2] - intNode[2]};

            const double pointNorm = norm2(point);
    
            //Loop over eta values
            for(uint32_t k = 0; k < 4; ++k)
            {
                const double timeInstant = deltaT * (offsetZ + double(blockIdx.z) + etaValues[k]);

                if(timeInstant <= 0)
                   continue;

                double tempValues[3][3] = {0.};

                //Add KR sub-kernel
                nuKR_linear(tempValues, point, pointNorm, timeInstant, velP, velS, lambda, mu, rho, intNormVect, baseCoeffMatrix, indSMcurr);

                //Add KL and KT sub-kernels
                nuKLT(tempValues, point, intNormVect, pointNorm, timeInstant, intNode, baseCoeffMatrix, baseCoeffVector, indSMcurr, velP, velS, lambda, mu);
                
                //Add to shared memory
                for(uint32_t i = 0; i < 3; ++i)
                    for(uint32_t j = 0; j < 3; ++j)
                        matrixSubBlock[sharedBaseInd][i][j] += extWeight * intWeight * etaCoeffs[k] * tempValues[i][j];
            }
        }
    }

    __syncthreads();

    //First step to bring y dimension to power of 2
    uint32_t yDim = 1 << (31 - __clz(blockDim.y));
    uint32_t sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;

    if(threadIdx.y + yDim < blockDim.y)
        for(uint32_t i = 0; i < 3; ++i)
            for(uint32_t j = 0; j < 3; ++j)
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

    __syncthreads();

    //Loop to reduce in y dimension
    yDim /= 2;
    while(yDim > 0)
    {
        sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;
        if(threadIdx.y < yDim)
            for(uint32_t i = 0; i < 3; ++i)
                for(uint32_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        yDim /= 2;
    }

    //Loop to reduce in x dimension (x size is assumed to be a power of 2)
    uint32_t xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;
        if((threadIdx.x < xDim) && (threadIdx.y == 0))
            for(uint32_t i = 0; i < 3; ++i)
                for(uint32_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    //Save in global memory
    if(threadIdx.x == 0 && threadIdx.y == 0)
        for(uint32_t i = 0; i < 3; ++i)
            for(uint32_t j = 0; j < 3; ++j)
            {
                const uint32_t ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                //matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(fabs(matrixSubBlock[0][i][j]) > 1e-14)
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiDeltaT;
            }
    
}






