__global__ void kernelWNodeNode(double *matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double DeltaT,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const double *normT, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff, 
                          const int offsetZ, const int numBlocks, const double *nodesMesh, const double maxLen)
{
    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
        return;
    

    //Place-holder per quellache sarà isBlockNull di questo kernel. Da trovare gli argomenti necessari!!
    //if (isBlockNull(nodesMesh, vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
        //return;


    extern __shared__ double matrixSubBlock[][3][3]; 
    int k, i, j, l, m, q, alpha, alphaTilde;



    const unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;



    int outerNode = blockIdx.x;
    int innerNode = blockIdx.y;

    //Suppose we have a matrix TrianglesPerNode, such that, each row represents the node with the corresponding index
    //(row 1 is the first node and row 2 is the second node etc etc), and in each row, the columns tell us the indexes of the 
    //triangles that have that node as a vertex (for node s, the s-th row would be \mathcal{T}_s)
    //Since not all nodes have the same number of triangles touching it, we will compute beforehand the maximum nuber of triangles that
    //touch a node for the entire mesh: MaxTrianglesPerNode. This means that TrianglesPerNode will be a NumNodes x MaxTrianglesPerNode
    //size matrix. if node s has less than MaxTrianglesPerNode touching it, then the remaining spaces in the row will be filled with zeros
    //However, note that TrianglesPerNode will be passed as a one dimentional array, and so the indexing needs to be adjusted.

    for (outerIndex = 0; outerIndex < MaxTrianglesPerNode; outerIndex++) //cycle on all triangles that have outer node as a vertex
    {
        currentOuterTriangleIndex = TrianglesPerNode[outerNode*MaxTrianglesPerNode + outerIndex]; //index of current outer triangle
        
        if (currentOuterTriangleIndex == 0) // if 0 this means we have already checked all triangles: exit the loop
            continue;

            for (innerIndex = 0; innerIndex < MaxTrianglesPerNode; innerIndex++) //cycle on all triangles that have inner node as a vertex
            {
                currentInnerTriangleIndex = TrianglesPerNode[innerNode*MaxTrianglesPerNode + innerIndex]; //index of current inner triangle
                
                if (currentInnerTriangleIndex == 0) // if 0 we have already checked all triangles
                    continue;
                
                // To find wich vertex the current nodes correspond to in regards to the currend inner and outer triangles, we use
                // indSMatrix. In this matrix each row corresponds to a node, and each column corresponds to a triangle. If node "s" is the
                // "n-th" vertex of triangle "m" then element (s,m) of indSMatrix will be n. This matrix is passed to the GPU as a 1-D array
                // which means that indexes have to be adjusted
             
                vertexNumberOfOutNodeInCurrOutTriang = indSMatrix[outerNode*NumTriangles + currentOuterTriangleIndex]; // this tells us which vertex is outerNode
                                                                                                                              // with regards to the current outer triangle
                                                                                                                              // so we can extract the correct shape function
    
                vertexNumberOfInnerNodeInCurrInnerTriang = indSMatrix[innerNode*NumTriangles + currentInnerTriangleIndex]; // this tells us which vertex is innerNode
                                                                                                                              // with regards to the current inner triangle
                                                                                                                              // so we can extract the correct shape function
                //here we will extract the  vertexes for the current inner and outer triangles

                    
                // here we will map the current gauss-hammer quadrature node on the inner and outer triangles using
                // threads and maybe another loop, and compute the weight

                //Here we will extract the current inner and outer triangles normal vectors

                //Here we will extract the current inner and outer values for the shape functions

                //Here we will compute the integral kernel 

                //Here we willimplement the quadrature procedure

                //Here we will sum all the results into matrixSubBlock
            }
    }


    



}