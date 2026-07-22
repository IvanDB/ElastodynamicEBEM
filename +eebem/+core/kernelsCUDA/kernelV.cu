#include "commonRoutines.cuh"

__global__ void kernelV(double* __restrict__ matrix, const double deltaT, const double velP, const double velS, const double const4PiRho,
                        const double* __restrict__ stdExtW, const double* __restrict__ stdExtNx, const double* __restrict__ stdExtNy, const double* __restrict__ stdExtNz, const int numPointExt,
                        const double* __restrict__ stdIntW, const double* __restrict__ stdIntNx, const double* __restrict__ stdIntNy, const double* __restrict__ stdIntNz,
                        const double* __restrict__ vertsT, const double* __restrict__ areeT, const int offsetZ, const int numBlocks)  
{
    //Skip diagonal blocks
    if(blockIdx.x == blockIdx.y)
        return;

    //Time instant check
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Spatial condition check
    if(isVBlockNull(vertsT, deltaT, velP, velS, offsetZ + blockIdx.z))
        return;

    //Shared memory initialization
    extern __shared__ double matrixSubBlock[][3][3];
    const size_t sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Field triangle verteces 
    double intVerts[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            intVerts[i][j] = vertsT[9*blockIdx.y + 3*j + i];
    
    //Field node
    const size_t intNodeIdx = threadIdx.x*blockDim.y + threadIdx.y;
    const double stdNodeI[] = {stdIntNx[intNodeIdx], stdIntNy[intNodeIdx], stdIntNz[intNodeIdx]};

    const double intNode[] = {stdNodeI[0] * intVerts[0][0] + stdNodeI[1] * intVerts[1][0] + stdNodeI[2] * intVerts[2][0],
                              stdNodeI[0] * intVerts[0][1] + stdNodeI[1] * intVerts[1][1] + stdNodeI[2] * intVerts[2][1],
                              stdNodeI[0] * intVerts[0][2] + stdNodeI[1] * intVerts[1][2] + stdNodeI[2] * intVerts[2][2]};
    
    //Field weight
    const double intWeight = stdIntW[threadIdx.y] * areeT[blockIdx.y];

    //Source triangle verteces
    double extVerts[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            extVerts[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    //Eta costants
    const double etaCoeffs[3] = {1, -2, 1};
    const double etaValues[3] = {-1, 0, 1};

    //Loop over source nodes
    for(size_t l = 0; l < numPointExt; ++l)
    {
        const double extWeight = stdExtW[l] * areeT[blockIdx.x];
      
        const double stdNodeE[] = {stdExtNx[l], stdExtNy[l], stdExtNz[l]};
        const double extNode[3] = {stdNodeE[0] * extVerts[0][0] + stdNodeE[1] * extVerts[1][0] + stdNodeE[2] * extVerts[2][0],
                                   stdNodeE[0] * extVerts[0][1] + stdNodeE[1] * extVerts[1][1] + stdNodeE[2] * extVerts[2][1],
                                   stdNodeE[0] * extVerts[0][2] + stdNodeE[1] * extVerts[1][2] + stdNodeE[2] * extVerts[2][2]};
        
        //Space vector
        const double point[3] = {intNode[0] - extNode[0], 
                                 intNode[1] - extNode[1], 
                                 intNode[2] - extNode[2]};

        const double pointNorm = norm2(point);
    
        //Loop over eta factor
        for(size_t k = 0; k < 3; ++k)
        {
            const double timeInstant = deltaT * (offsetZ + double(blockIdx.z) + etaValues[k]);
            
            if(timeInstant < 0)
                continue;
            
            double tempValues[3][3] = {0.};
            nuV(tempValues, point, pointNorm, timeInstant, velP, velS);
            
            //Add to shared memory
            for(size_t i = 0; i < 3; ++i)
                for(size_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += extWeight * intWeight * etaCoeffs[k] * tempValues[i][j];
    
        }
    }

    __syncthreads();
    
    size_t sharedOffInd;
    
    //First step to reduce y dimension to power of 2
    size_t yDim = pow(2.0, (int) floor(log2((float) blockDim.y)));
    sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;

    if(threadIdx.y + yDim < blockDim.y)
        for(size_t i = 0; i < 3; ++i)
            for(size_t j = 0; j < 3; ++j)
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
    
    __syncthreads();

    //Loop to reduce in y dimension
    yDim /= 2;
    while(yDim > 0)
    {
        sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;
        if(threadIdx.y < yDim)
            for(size_t i = 0; i < 3; ++i)
                for(size_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        yDim /= 2;
    }

    //Loop to reduce in x dimension (x size is assumed to be a power of 2)
    size_t xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;
        if(threadIdx.x < xDim)
            for(size_t i = 0; i < 3; ++i)
                for(size_t j = 0; j < 3; ++j)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    //Save in global memory
    if(threadIdx.x == 0 && threadIdx.y == 0)
        for(size_t i = 0; i < 3; ++i)
            for(size_t j = 0; j < 3; ++j)
            {
                const size_t ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                // matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiRho;
            }

}
