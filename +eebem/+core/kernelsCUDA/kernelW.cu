__global__ void kernelK(double *matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double DeltaT,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const double *normT, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff, 
                          const int offsetZ, const int numBlocks, const double *nodesMesh, const double maxLen)
{   
    //Evito i blocchi diagonali
    if(blockIdx.x == blockIdx.y)
        return;
    
    
    
    //Solito controllo per non sforare con il tempo
    if(offsetZ + blockIdx.z >= numBlocks)
        return;
    
    //Place-holder per quellache sarà isBlockNull di questo kernel. Da trovare gli argomenti necessari!!
    //if (isBlockNull(nodesMesh, vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
        //return;


    extern __shared__ double matrixSubBlock[][9][9]; 
    int k, i, j, l, m, q;

    const unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;

    // Inizializzazione shared memory
    for (i = 0; i < 9; i++)//Anche qui 9x9 forse
        for (j = 0; j < 9; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Dichiarazione e inizializzazione tensore di permutazione E
    const int tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};
    
   //Dichiarazione variabili
    double vertInnerTriangle[3][3];
    double vertOuterTriangle[3][3];
    double temporaryStandardNode[3];
    
    //Calcolo peso nodo iterno GHC corrente
    double currentInnerWeight = stdGHCw[threadIdx.y] * areeT[blockIdx.y];

    // Lettura vertici triangolo interno e esterno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertInnerTriangle[i][j] = vertsT[9*blockIdx.y + 3*j + i];
            vertOuterTriangle[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    //Lettura delle normali correnti del triangolo interno ed esterno
    double normInt[3], normExt[3];

    for (i = 0; i < 3; i++) 
        normInt[i] = normT[3*blockIdx.y + i];
        normExt[i] = normT[3*blockIdx.x + i];
    
    
    // Lettura coordinate nodo interno GHC corrente su triangolo stardard
    temporaryStandardNode[0] = stdGHCnx[threadIdx.x*blockDim.y + threadIdx.y];
    temporaryStandardNode[1] = stdGHCny[threadIdx.x*blockDim.y + threadIdx.y];
    temporaryStandardNode[2] = stdGHCnz[threadIdx.x*blockDim.y + threadIdx.y];

    //Mappaggio nodo GHC corrente su triangolo interno
    double currentInnerNode[3];
    currentInnerNode[0] = temporaryStandardNode[0] * vertInnerTriangle[0][0] + temporaryStandardNode[1] * vertInnerTriangle[1][0] + temporaryStandardNode[2] * vertInnerTriangle[2][0];
    currentInnerNode[1] = temporaryStandardNode[0] * vertInnerTriangle[0][1] + temporaryStandardNode[1] * vertInnerTriangle[1][1] + temporaryStandardNode[2] * vertInnerTriangle[2][1];
    currentInnerNode[2] = temporaryStandardNode[0] * vertInnerTriangle[0][2] + temporaryStandardNode[1] * vertInnerTriangle[1][2] + temporaryStandardNode[2] * vertInnerTriangle[2][2];


    //Ciclo sui numPointExt nodi del triangolo esterno
    for(l = 0; l < numPointExt; l++)
    {
        // Calcolo peso nodo GH triangolo esterno
        double currentExternWeight = stdGHw[l] * areeT[blockIdx.x];
        
        //Lettura nodo GHC esterno corrente su triangolo standard
        temporaryStandardNode[0] = stdGHnx[l];
        temporaryStandardNode[1] = stdGHny[l];
        temporaryStandardNode[2] = stdGHnz[l];

        //Mappaggio nodo GHC corrente su triangolo esterno
        double currentOuterNode[3];
        currentOuterNode[0] = temporaryStandardNode[0] * vertOuterTriangle[0][0] + temporaryStandardNode[1] * vertOuterTriangle[1][0] + temporaryStandardNode[2] * vertOuterTriangle[2][0];
        currentOuterNode[1] = temporaryStandardNode[0] * vertOuterTriangle[0][1] + temporaryStandardNode[1] * vertOuterTriangle[1][1] + temporaryStandardNode[2] * vertOuterTriangle[2][1];
        currentOuterNode[2] = temporaryStandardNode[0] * vertOuterTriangle[0][2] + temporaryStandardNode[1] * vertOuterTriangle[1][2] + temporaryStandardNode[2] * vertOuterTriangle[2][2];
    
        //Calcolo coordinate vettore differenza
        double point[3];
        rDifference[0] = currentOuterNode[0] - currentInnerNode[0];
        rDifference[1] = currentOuterNode[1] - currentInnerNode[1];
        rDifference[2] = currentOuterNode[2] - currentInnerNode[2];

        //Calcolo norma vettore differenza
        double rNorm = sqrt(rDifference[0]*rDifference[0] + rDifference[1]*rDifference[1] + rDifference[2]*rDifference[2]);
}