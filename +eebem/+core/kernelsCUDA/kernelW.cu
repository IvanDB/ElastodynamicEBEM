static __device__ inline bool isBlockNull_W(const double *nodesMesh, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen);

static __device__ inline double dot(const double vettA[3], const double vettB[3]);
static __device__ inline void cross(const double vettA[3], const double vettB[3], double vettC[3]);

static __device__ inline double VFunction(const double triangleCoeffs[3], const double triangleMatrix[3][3], const int vertexInd, const double xi[3]);
static __device__ inline double VTildeFunction(const double triangleCoeffs[3], const double triangleMatrix[3][3], const int vertexInd, const double x[3]);

static __device__ inline void kernelCalc(double kernel[3][3], const double t, const double VTildeX, const double VXi,
                                         const double VAlphaS[3], const double VTildeAlphaS[3], const double rVector[3],
                                         const double r, const double n[3], const double v[3], const double pi, const double mu, const double lambda,
                                         const double rho, const double cS, const double cP);

static __device__ inline double kernelTCalc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                        const double rVector[3], const double r, const double n[3], const double v[3], const double pi, const double mu,
                                           const double lambda, const double rho, const double cS, const double cP);
static __device__ inline double kernelRCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3],
                                           const double rVector[3], const double r, const double pi, const double mu, const double lambda,
                                           const double rho, const double cS, const double cP);

static __device__ inline double kernelTDelta1Calc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                                 const double rVector[3], const double r, const double n[3], const double v[3], const double pi, const double mu,
                                                 const double lambda, const double cS, const double cP);
static __device__ inline double kernelTDelta2Calc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                                 const double r, const double n[3], const double v[3], const double pi, const double mu,
                                                 const double lambda, const double rho, const double cS, const double cP);
static __device__ inline double kernelTHCalc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                            const double rVector[3], const double r, const double n[3], const double v[3], const double pi,
                                            const double mu, const double lambda, const double cS, const double cP);

static __device__ inline double kernelRDeltaCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3], const double rVector[3],
                                                const double r, const double pi, const double mu, const double lambda, const double rho, const double cS, const double cP);
static __device__ inline double kernelRHCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3],
                                            const  double rVector[3], const double r, const double pi, const double mu, const double lambda, const double rho,
                                            const double cS, const double cP);

__global__ void kernelW(double *matrix, const double deltaT, const double cP, const double cS, const double lambda, const double mu, const double rho, const double pi,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const double *normT, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff, 
                          const int offsetZ, const int numBlocks, const double *nodesMesh, const double maxLen, const int MaxTrianglesPerNode, const int *TrianglesPerNode)
{
    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
    {
        return;
    }
    
    // Check which 3x3 blocks are null
    if (isBlockNull_W(nodesMesh, deltaT, cP, cS, offsetZ + blockIdx.z, maxLen))
    {
        return;
    }
    


    extern __shared__ double matrixSubBlock[][3][3]; 
    int i, j, k, l;
    
    


    const unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    // Inizializzazione shared memory
    #pragma unroll
    for (i = 0; i < 3; i++)
    {
        #pragma unroll
        for (j = 0; j < 3; j++)
        {
            matrixSubBlock[sharedBaseInd][i][j] = 0;
        }
    }


    int outerNode = blockIdx.x;
    int innerNode = blockIdx.y;
    int outerIndex;
    int innerIndex;
    int currentOuterTriangleIndex;
    int currentInnerTriangleIndex;
    int currentVertexNumberOfOutNode;
    int currentVertexNumberOfInNode;

    // We have a matrix TrianglesPerNode, such that, each row represents the node with the corresponding index
    //(row 1 is the first node and row 2 is the second node etc etc), and in each row, the columns tell us the indexes of the 
    //triangles that have that node as a vertex (for node s, the s-th row would be \mathcal{T}_s)
    //Since not all nodes have the same number of triangles touching it, we will compute beforehand the maximum nuber of triangles that
    //touch a node for the entire mesh: MaxTrianglesPerNode. This means that TrianglesPerNode will be a NumNodes x MaxTrianglesPerNode
    //size matrix. if node s has less than MaxTrianglesPerNode touching it, then the remaining spaces in the row will be filled with zeros
    //However, note that TrianglesPerNode will be passed as a one dimentional array, and so the indexing needs to be adjusted.

    for (outerIndex = 0; outerIndex < MaxTrianglesPerNode; outerIndex++) //cycle on all triangles that have outer node as a vertex
    {
        currentOuterTriangleIndex = TrianglesPerNode[outerIndex*gridDim.x + outerNode]; //index of current outer triangle
        
        if (currentOuterTriangleIndex == 0) // if 0 this means we have already checked all triangles: exit the loop
        {
            break;
        }
        for (innerIndex = 0; innerIndex < MaxTrianglesPerNode; innerIndex++) //cycle on all triangles that have inner node as a vertex
        {
            currentInnerTriangleIndex = TrianglesPerNode[innerIndex*gridDim.y + innerNode]; //index of current inner triangle
                
            if (currentInnerTriangleIndex == 0)
            {                           // if 0 we have already checked all triangles
                break;
            }
            if (currentInnerTriangleIndex == currentOuterTriangleIndex)
            {
                continue;                                               // Singular integrations are executed on CPU
            }    
                // To find wich vertex the current nodes correspond to in regards to the currend inner and outer triangles, we use
                // indSMatrix. In this matrix each row corresponds to a node, and each column corresponds to a triangle that touches that Node, and it in in correspondance
                // with TrianglesPerNode. So the s-th row tells us exactly which vertex is node s, revalive to the trianglse that are present in the s-th row of TrianglesPerNode
             
            currentVertexNumberOfOutNode = indSMmatrix[outerIndex*gridDim.x + outerNode]; // this tells us which vertex is outerNode
                                                                                                                              // with regards to the current outer triangle
                                                                                                                              // so we can extract the correct shape function
    
            currentVertexNumberOfInNode = indSMmatrix[innerIndex*gridDim.y + innerNode]; // this tells us which vertex is innerNode
                                                                                                                              // with regards to the current inner triangle
                                                                                                                              // so we can extract the correct shape function
                


            // Now we extract the current outer triangle vertexes coordinates
            double currentOuterVertexes[3][3];
            #pragma unroll
            for (i=0; i<3; i++)
            {
                #pragma unroll
                for(j=0; j<3; j++)
                {
                    currentOuterVertexes[i][j] = vertsT[9*(currentOuterTriangleIndex - 1) + 3*j + i];// remember to pass vertsT the right way
                }
            }


            // Now we extract the current inner triangle vertexes coordinates
            double currentInnerVertexes[3][3];
            #pragma unroll
            for (i=0; i<3; i++)
            {
                #pragma unroll
                for(j=0; j<3; j++)
                {
                    currentInnerVertexes[i][j] = vertsT[9*(currentInnerTriangleIndex - 1) + 3*j + i];
                }
            }
            
            // Extract the current outer triangle normal vector
            double currentOuterNormal[3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                currentOuterNormal[i] = normT[3*(currentOuterTriangleIndex - 1) + i];
            }


            // Extract the current inner triangle normal vector
            double currentInnerNormal[3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                currentInnerNormal[i] = normT[3*(currentInnerTriangleIndex - 1) + i];
            }
                      
            
            //Extraction of outer triangle matrix and piece of shape function that only depends on the triangle 
            double currentOuterTriangleMatrix[3][3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                #pragma unroll
                for (j = 0; j < 3; j++)
                {
                    currentOuterTriangleMatrix[i][j] = matCoeff[9*(currentOuterTriangleIndex - 1) + 3*j + i];
                }
            }

            double currentOuterTriangleCoeffs[3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                currentOuterTriangleCoeffs[i] = vetCoeff[3*(currentOuterTriangleIndex - 1) + i];
            }



             //Extraction of inner triangle matrix and piece of shape function that only depends on the triangle 
            double currentInnerTriangleMatrix[3][3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                #pragma unroll
                for (j = 0; j < 3; j++)
                {
                    currentInnerTriangleMatrix[i][j] = matCoeff[9*(currentInnerTriangleIndex - 1) + 3*j + i];
                }
            }

            double currentInnerTriangleCoeffs[3];
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                currentInnerTriangleCoeffs[i] = vetCoeff[3*(currentInnerTriangleIndex - 1) + i];
            }
            

            double currentOuterVVector[3];
            double tempTriangleMatrixRow[3];
            #pragma unroll
            for (i=0; i<3; i++)
            {
                tempTriangleMatrixRow[i] = currentOuterTriangleMatrix[currentVertexNumberOfOutNode-1][i];
            }
            cross(currentOuterNormal, tempTriangleMatrixRow, currentOuterVVector);// This is \bm{V_{\tilde{\alpha}}^{\tilde{s}}}


            double currentInnerVVector[3];
            #pragma unroll
            for (i=0; i<3; i++)
            {
                tempTriangleMatrixRow[i] = currentInnerTriangleMatrix[currentVertexNumberOfInNode-1][i];
            }
            cross(currentInnerNormal, tempTriangleMatrixRow, currentInnerVVector); // This is \bm{V_{\alpha}^{s}}
            
            

            // "current" standard Composite Gauss-Hammer quadrature node
            double standardGaussNode[3];
            
            standardGaussNode[0] = stdGHCnx[threadIdx.x*blockDim.y + threadIdx.y];
            standardGaussNode[1] = stdGHCny[threadIdx.x*blockDim.y + threadIdx.y];
            standardGaussNode[2] = stdGHCnz[threadIdx.x*blockDim.y + threadIdx.y];

            // mapping "current" standard quadrature node onto "current" inner quadrature nodo on the actual current inner triangle
            double InnerGaussNode[3];
            InnerGaussNode[0] = standardGaussNode[0] * currentInnerVertexes[0][0] + standardGaussNode[1] * currentInnerVertexes[1][0] + standardGaussNode[2] * currentInnerVertexes[2][0];
            InnerGaussNode[1] = standardGaussNode[0] * currentInnerVertexes[0][1] + standardGaussNode[1] * currentInnerVertexes[1][1] + standardGaussNode[2] * currentInnerVertexes[2][1];
            InnerGaussNode[2] = standardGaussNode[0] * currentInnerVertexes[0][2] + standardGaussNode[1] * currentInnerVertexes[1][2] + standardGaussNode[2] * currentInnerVertexes[2][2];
                
            // wheight of "current" inner gauss quadrature node
            double innerGaussWeight = stdGHCw[threadIdx.y] * areeT[currentInnerTriangleIndex-1];
            
            // compute shape function on inner node
            double currentVXi;
            currentVXi = VFunction(currentInnerTriangleCoeffs, currentInnerTriangleMatrix, currentVertexNumberOfInNode, InnerGaussNode); 
            
            // Loop on outer composite Gauss-Hammer Nodes on current outer triangle
            
            for (l = 0; l < numPointExt; l++)
            {
                // Reading curren standard outer coposite Gauss_Hammer node
                standardGaussNode[0] = stdGHnx[l];
                standardGaussNode[1] = stdGHny[l];
                standardGaussNode[2] = stdGHnz[l];

                // Mapping current node on current triangle

                double currentOuterGaussNode[3];
                currentOuterGaussNode[0] = standardGaussNode[0] * currentOuterVertexes[0][0] + standardGaussNode[1] * currentOuterVertexes[1][0] + standardGaussNode[2] * currentOuterVertexes[2][0];
                currentOuterGaussNode[1] = standardGaussNode[0] * currentOuterVertexes[0][1] + standardGaussNode[1] * currentOuterVertexes[1][1] + standardGaussNode[2] * currentOuterVertexes[2][1];
                currentOuterGaussNode[2] = standardGaussNode[0] * currentOuterVertexes[0][2] + standardGaussNode[1] * currentOuterVertexes[1][2] + standardGaussNode[2] * currentOuterVertexes[2][2];



                // Getting current outer composite Gauss-Hammer weight
                double currentOuterGaussWeight = stdGHw[l] * areeT[currentOuterTriangleIndex-1];

                // compute shape function on current outer node
                double currentVTildeX;
                currentVTildeX = VTildeFunction(currentOuterTriangleCoeffs, currentOuterTriangleMatrix, currentVertexNumberOfOutNode, currentOuterGaussNode);
                
                // r = x - \xi
                double current_rVector[3];
                current_rVector[0] = currentOuterGaussNode[0] - InnerGaussNode[0];
                current_rVector[1] = currentOuterGaussNode[1] - InnerGaussNode[1];
                current_rVector[2] = currentOuterGaussNode[2] - InnerGaussNode[2];

                // r norm
                double current_rNorm = sqrt(current_rVector[0]*current_rVector[0] + current_rVector[1]*current_rVector[1] + current_rVector[2]*current_rVector[2]);

                double timeInstant = 0;
                double kernelValues[3][3];
                const int kernelCoeffs[4] = {1, -3, 3, -1};
                const int timeShift[4] = {-2, -1, 0, 1};
                #pragma unroll
                for (k = 0; k < 4; k++)
                {
                    timeInstant = deltaT * (offsetZ + double(blockIdx.z) +timeShift[k]); // (l+\eta)\Delta t
                        
                    // skip negative time
                    if(timeInstant <= 0)
                    {
                        continue;
                    }

                    // initializing kernel components
                    #pragma unroll
                    for (i = 0; i < 3; i++)
                    {
                        #pragma unroll
                        for (j = 0; j < 3; j++)
                        {
                            kernelValues[i][j] = 0;
                        }
                    }
                    // void function that will fill kernelValues
                    kernelCalc(kernelValues, timeInstant, currentVTildeX, currentVXi, currentInnerVVector, currentOuterVVector, current_rVector, current_rNorm,
                               currentInnerNormal, currentOuterNormal, pi, mu, lambda, rho, cS, cP);
                
                    //Somma pesata dei valori del nucleo alla shared memory
                    #pragma unroll
                    for (i = 0; i < 3; i++)
                    {
                        #pragma unroll
                        for (j = 0; j < 3; j++)
                        {
                            matrixSubBlock[sharedBaseInd][i][j] += currentOuterGaussWeight * innerGaussWeight * kernelCoeffs[k]/(deltaT*deltaT) * kernelValues[i][j];
                        }
                    }
                }
            }       
        }
    }


    __syncthreads();

    // Starting Parallel Reduction
    unsigned int xDim, yDim, sharedOffInd;

    // Dimention along y is not guaranteed to be a power of two
    // First reduction to the biggest power of two smaller than y
    yDim = pow(2.0, (int) floor(log2((float) blockDim.y)));
    sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;

    if(threadIdx.y + yDim < blockDim.y)
    {
        #pragma unroll
        for (i = 0; i < 3; i++)
        {
            #pragma unroll
            for (j = 0; j < 3; j++)
            {
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
            }
        }
    }

    __syncthreads();

    //Iteration along y
    yDim /= 2;
    while(yDim > 0)
    {
        sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;
        if(threadIdx.y < yDim)
        {
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                #pragma unroll
                for (j = 0; j < 3; j++)
                {
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
                }
            }
        }

        __syncthreads();
        yDim /= 2;
    }

    //Iteration along x (always a power of two)
    xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;
        if(threadIdx.x < xDim && threadIdx.y == 0)
        {
            #pragma unroll
            for (i = 0; i < 3; i++)
            {
                #pragma unroll
                for (j = 0; j < 3; j++)
                {
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
                }
            }
        }

        __syncthreads();
        xDim /= 2;
    }

    // Saving on global memory
    unsigned long ind;
    if(threadIdx.x == 0 && threadIdx.y == 0)
    {
        #pragma unroll
        for (i = 0; i < 3; i++)
        {
            #pragma unroll
            for (j = 0; j < 3; j++)
            {
                ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                //matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(fabs(matrixSubBlock[0][i][j]) > 1e-14)
                {
                    matrix[ind] += matrixSubBlock[0][i][j];
                }
            }
        }
    }


    



}


static __device__ inline bool isBlockNull_W(const double *nodesMesh, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
{
    
    double outerNode[3];
    double innerNode[3];

    // blockIdx.x = outerNode (\tilde{s}), blockIdx.y = innerNode (s)
    #pragma unroll
    for (int i = 0; i < 3; i++) 
    {
        outerNode[i] = nodesMesh[3 * blockIdx.x + i];
        innerNode[i]  = nodesMesh[3 * blockIdx.y + i];
    }

    // \tilde{s} - s
    double vettDist[3];
    vettDist[0] = outerNode[0] - innerNode[0];
    vettDist[1] = outerNode[1] - innerNode[1];
    vettDist[2] = outerNode[2] - innerNode[2];
    
    // node - to - node distance
    double d_NN = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]);

    // max and min distance estimate between shape function domains
    double distMax_Patch = d_NN + 2.0 * maxLen;
    double distMin_Patch = d_NN - 2.0 * maxLen;

    // S wave already passed
    bool wavePassed = ((indTemp - 2) * cS * deltaT > distMax_Patch);
    
    // P wave not reached
    bool waveNotReached = ((indTemp + 1) * cP * deltaT < distMin_Patch);

    // if either one is true skip block
    return wavePassed || waveNotReached;
}






// dot product
static __device__ inline double dot(const double vettA[3], const double vettB[3])
{
    return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];
}


// cross product
static __device__ inline void cross(const double vettA[3], const double vettB[3], double vettC[3])
{
    vettC[0] = vettA[1]*vettB[2] - vettA[2]*vettB[1];
    vettC[1] = vettA[2]*vettB[0] - vettA[0]*vettB[2];
    vettC[2] = vettA[0]*vettB[1] - vettA[1]*vettB[0];
}

// VFunction to evaluate V_{\alpha}^{s}(\xi)

static __device__ inline double VFunction(const double triangleCoeffs[3], const double triangleMatrix[3][3], const int vertexInd, const double xi[3])
{
    double eVec[3] = {0.0, 0.0, 0.0};
    eVec[vertexInd - 1] = 1.0;
    double shapeVec[3];

    shapeVec[0] = triangleCoeffs[0] + triangleMatrix[0][0] * xi[0] + triangleMatrix[0][1] * xi[1] + triangleMatrix[0][2] * xi[2];
    shapeVec[1] = triangleCoeffs[1] + triangleMatrix[1][0] * xi[0] + triangleMatrix[1][1] * xi[1] + triangleMatrix[1][2] * xi[2];
    shapeVec[2] = triangleCoeffs[2] + triangleMatrix[2][0] * xi[0] + triangleMatrix[2][1] * xi[1] + triangleMatrix[2][2] * xi[2];

    return dot(shapeVec, eVec);
}


static __device__ inline double VTildeFunction(const double triangleCoeffs[3], const double triangleMatrix[3][3], const int vertexInd, const double x[3])
{
    double eVec[3] = {0.0, 0.0, 0.0};
    eVec[vertexInd - 1] = 1.0;
    double shapeVec[3];

    shapeVec[0] = triangleCoeffs[0] + triangleMatrix[0][0] * x[0] + triangleMatrix[0][1] * x[1] + triangleMatrix[0][2] * x[2];
    shapeVec[1] = triangleCoeffs[1] + triangleMatrix[1][0] * x[0] + triangleMatrix[1][1] * x[1] + triangleMatrix[1][2] * x[2];
    shapeVec[2] = triangleCoeffs[2] + triangleMatrix[2][0] * x[0] + triangleMatrix[2][1] * x[1] + triangleMatrix[2][2] * x[2];

    return dot(shapeVec, eVec);
}




//function arguments:  kernelValues, timeInstant, currentVTildeX, currentVXi, currentInnerVVector, currentOuterVVector, current_rVector, current_rNorm,
                               // currentInnerNormal, currentOuterNormal, pi, mu, lambda, rho, cS, cP

static __device__ inline void kernelCalc(double kernel[3][3], const double t, const double VTildeX, const double VXi,
                                         const double VAlphaS[3], const double VTildeAlphaS[3], const double rVector[3],
                                         const double r, const double n[3], const double v[3], const double pi, const double mu, const double lambda,
                                         const double rho, const double cS, const double cP)
{
    int i, k;
    // declare and initialize kernel T components
    double kernelT[3][3];
    #pragma unroll
    for (i = 0; i < 3; i++)
    {
        #pragma unroll
        for (k = 0; k < 3; k++)
        {
            kernelT[i][k] = 0;
        }
    }

    // declare and initialize kernel R components
    double kernelR[3][3];

    #pragma unroll
    for (i = 0; i < 3; i++)
    {
        #pragma unroll
        for (k = 0; k < 3; k++)
        {
            kernelR[i][k] = 0;
        }
    }

    // compute and add kernels
    #pragma unroll
    for (i = 0; i < 3; i++)
    {
        #pragma unroll
        for(k = 0; k < 3; k++)
        {
            kernelT[i][k] = kernelTCalc(i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, rho, cS, cP);
            kernelR[i][k] = kernelRCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP);
            kernel[i][k] = kernelT[i][k] + kernelR[i][k];
        }
    }
}



// function arguments: i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, rho, cS, cP
static __device__ inline double kernelTCalc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                        const double rVector[3], const double r, const double n[3], const double v[3], const double pi, const double mu,
                                           const double lambda, const double rho, const double cS, const double cP)
{

    // compute components
    double kernelTDelta1_ik = kernelTDelta1Calc(i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, cS, cP);
    double kernelTDelta2_ik = kernelTDelta2Calc(i, k, t, VTildeX, VXi, r, n, v, pi, mu, lambda, rho, cS, cP);
    double kernelTH_ik = kernelTHCalc(i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, cS, cP);

    //add together and return

    return kernelTDelta1_ik + kernelTDelta2_ik + kernelTH_ik;
}


// function arguments: i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP
static __device__ inline double kernelRCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3],
                                           const double rVector[3], const double r, const double pi, const double mu, const double lambda,
                                           const double rho, const double cS, const double cP)
{

    // compute components
    double kernelRDelta_ik = kernelRDeltaCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP);
    double kernelRH_ik = kernelRHCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP);
    // add together and return
    return kernelRDelta_ik + kernelRH_ik;
}

// function arguments: i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, cS, cP
static __device__ inline double kernelTDelta1Calc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                                 const double rVector[3], const double r, const double n[3], const double v[3], const double pi, const double mu,
                                                 const double lambda, const double cS, const double cP)
{
   
    // Heavyside functions and Kronecker delta
    double H_P = (t - (r / cP) > 0.0) ? 1.0 : 0.0;
    double H_S = (t - (r / cS) > 0.0) ? 1.0 : 0.0;
    double delta_ik = (i == k) ? 1.0 : 0.0;

    // double time integrals
    double IDeltaS1 = H_S/(cS*cS);
    double IDeltaP1 = H_P/(cP*cP);

    // scalar products
    double Rv = dot(rVector,v);
    double Rn = dot(rVector,n); 
    double Nv = dot(n,v);

    // constant 
    double constTerm = 1.0/(2.0*pi*(lambda+mu));

    return VTildeX*VXi*constTerm*(lambda*lambda*((n[k]*v[i])/(2.0*r)) + lambda*mu*((rVector[i]*n[k]*Rv + rVector[k]*v[i]*Rn)/(r*r*r)) + 
                                                    mu*mu*((rVector[k]*n[i]*Rv + rVector[i]*v[k]*Rn + delta_ik*Rv*Rn + rVector[i]*rVector[k]*Nv)/(2.0*r*r*r)))*(IDeltaS1-IDeltaP1);
}

// function arguments: i, k, t, VTildeX, VXi, r, n, v, pi, mu, lambda, rho, cS, cP
static __device__ inline double kernelTDelta2Calc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                                 const double r, const double n[3], const double v[3], const double pi, const double mu,
                                                 const double lambda, const double rho, const double cS, const double cP)
{
    // Heavyside functions and Kronecker delta
    double H_P = (t - (r/cP) > 0.0) ? 1.0 : 0.0;
    double H_S = (t - (r/cS) > 0.0) ? 1.0 : 0.0;
    double delta_ik = (i == k) ? 1.0 : 0.0;

    // double time integrals
    double IDeltaS2 = H_S/(cS*cS*cS*cS); 
    double IDeltaP2 = H_P/(cP*cP*cP*cP);

    // scalar products
    double Nv = dot(n,v);

    //constant term
    double constTerm = -(lambda*mu + 2.0*mu*mu)/(4.0*pi*rho*(lambda+mu));

    return VTildeX*VXi*constTerm*((lambda*n[k]*v[i] + mu*n[i]*v[k] + mu*delta_ik*(Nv))/(r))*(IDeltaS2-IDeltaP2);
}

// function arguments i, k, t, VTildeX, VXi, rVector, r, n, v, pi, mu, lambda, cS, cP
static __device__ inline double kernelTHCalc(const int i, const int k, const double t, const double VTildeX, const double VXi,
                                            const double rVector[3], const double r, const double n[3], const double v[3], const double pi,
                                            const double mu, const double lambda, const double cS, const double cP)
{
    // Heavyside functions and Kronecker delta
    double H_P = (t - (r/cP) > 0.0) ? 1.0 : 0.0;
    double H_S = (t - (r/cS) > 0.0) ? 1.0 : 0.0;
    double delta_ik = (i == k) ? 1.0 : 0.0;

    // double time integrals
    double IHS = H_S*(((t*t) - ((r/cS)*(r/cS)))/2.0);
    double IHP = H_P*(((t*t) - ((r/cP)*(r/cP)))/2.0);

    // scalar products
    double Rn = dot(rVector,n);
    double Rv = dot(rVector,v);
    double Nv = dot(n,v);

    // constant term
    double constTerm = mu/(2.0*pi*(lambda+mu));

    return VTildeX*VXi*constTerm*(3.0*lambda*((rVector[k]*v[i]*Rn + rVector[i]*n[k]*Rv)/(r*r*r*r*r)) - 2.0*lambda*((n[k]*v[i])/(r*r*r)) 
                                  + 3.0*mu*((rVector[k]*n[i]*Rv + rVector[i]*v[k]*Rn + delta_ik*Rv*Rn + rVector[i]*rVector[k]*Nv)/(2.0*r*r*r*r*r)) - mu*((n[i]*v[k] + delta_ik*Nv)/(r*r*r)))*(IHS-IHP);
}

// function arguments: i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP
static __device__ inline double kernelRDeltaCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3], const double rVector[3],
                                                const double r, const double pi, const double mu, const double lambda, const double rho, const double cS, const double cP)
{
    // Heavyside functions and Kronecker delta
    double H_P = (t - (r/cP) > 0.0) ? 1.0 : 0.0;
    double H_S = (t - (r/cS) > 0.0) ? 1.0 : 0.0;
    double delta_ik = (i == k) ? 1.0 : 0.0;

    // double time integrals
    double IDeltaRGS = H_S*(((t-(r/cS))*(t-(r/cS)))/(2*cS*cS));
    double IDeltaRGP = H_P*(((t-(r/cP))*(t-(r/cP)))/(2*cP*cP));

    // dot and cross products
    double A[3];
    double B[3];

    cross(rVector, VTildeAlphaS, A); 
    cross(rVector, VAlphaS, B); 

    double Vdot = dot(VTildeAlphaS, VAlphaS); 

    // constant term 
    double constTerm = -((mu*mu)/(4.0*pi*rho*(lambda+mu)));


    return constTerm*((2.0*lambda / (r*r*r)) * A[i] * B[k] + ((lambda + 2.0*mu) / (r*r*r)) * B[i] * A[k] + (lambda + 2.0*mu) * Vdot * (delta_ik/r - (rVector[i]*rVector[k])/(r*r*r)))*(IDeltaRGS - IDeltaRGP);

}

// function arguments: i, k, t, VAlphaS, VTildeAlphaS, rVector, r, pi, mu, lambda, rho, cS, cP
static __device__ inline double kernelRHCalc(const int i, const int k, const double t, const double VAlphaS[3], const double VTildeAlphaS[3],
                                            const  double rVector[3], const double r, const double pi, const double mu, const double lambda, const double rho,
                                            const double cS, const double cP)
{
    // Heavyside functions and Kronecker delta
    double H_P = (t - (r/cP) > 0.0) ? 1.0 : 0.0;
    double H_S = (t - (r/cS) > 0.0) ? 1.0 : 0.0;
    double delta_ik = (i == k) ? 1.0 : 0.0;

    // double time integrals
    double IHRGS = H_S*(((t-(r/cS))*(t-(r/cS))*(t-(r/cS))*(t-(r/cS)))/(24.0) + r*((t-(r/cS))*(t-(r/cS))*(t-(r/cS)))/(6.0*cS));
    double IHRGP = H_P*(((t-(r/cP))*(t-(r/cP))*(t-(r/cP))*(t-(r/cP)))/(24.0) + r*((t-(r/cP))*(t-(r/cP))*(t-(r/cP)))/(6.0*cP));


    // dot and cross products
    double A[3];
    double B[3];

    cross(rVector, VTildeAlphaS, A); 
    cross(rVector, VAlphaS, B); 

    double Vdot = dot(VTildeAlphaS, VAlphaS); 

    // constant term
    double constTerm = -((mu*mu)/(4.0*pi*rho*(lambda+mu)));

    // kernel terms
    double term_r5 = (3.0/(r*r*r*r*r)) * ( 2.0*lambda * A[i] * B[k] + (lambda + 2.0*mu) * B[i] * A[k] - (lambda + 2.0*mu) * Vdot * rVector[i] * rVector[k]);

    double term_r3 = (1/(r*r*r)) * ( -2.0*lambda * delta_ik * Vdot + 2.0*lambda * VAlphaS[i] * VTildeAlphaS[k] + (lambda + 2.0*mu) * VTildeAlphaS[i] * VAlphaS[k]);

    return constTerm * (term_r5 + term_r3) * (IHRGS - IHRGP);
}