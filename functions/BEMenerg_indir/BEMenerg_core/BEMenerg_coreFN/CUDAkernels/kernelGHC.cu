__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, const double cP, const double cS);
__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp);


__global__ void kernelGHC(double *matrix, const double deltaT, const double velP, const double velS, const double const4PiRho,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const int offsetZ, const int numBlocks)  
{
    //Evito i blocchi diagonali
    if(blockIdx.x == blockIdx.y)
        return;

    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Controllo condizione teorica necessità di calcolo
    if (isBlockNull(vertsT, deltaT, velP, velS, offsetZ + blockIdx.z))
        return;

    extern __shared__ double matrixSubBlock[][3][3];
    int k, i, j, l;

    const unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Dichiarazione variabili
    double vertsTemp[3][3];
    double nodoTemp[3];
    
    //Calcolo peso nodo GHC corrente
    double pesoInt = stdGHCw[threadIdx.y] * areeT[blockIdx.y];

    // Lettura vertici triangolo interno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTemp[i][j] = vertsT[9*blockIdx.y + 3*j + i];
    
    // Lettura coordinate nodo GHC corrente su triangolo stardard
    nodoTemp[0] = stdGHCnx[threadIdx.x*blockDim.y + threadIdx.y];
    nodoTemp[1] = stdGHCny[threadIdx.x*blockDim.y + threadIdx.y];
    nodoTemp[2] = stdGHCnz[threadIdx.x*blockDim.y + threadIdx.y];

    //Mappaggio nodo GHC corrente su triangolo interno
    double nodoGHCcurr[3];
    nodoGHCcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
    nodoGHCcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
    nodoGHCcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];

    //Lettura vertici triangolo esterno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTemp[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    //Ciclo sui numPointExt nodi del triangolo esterno
    for(l = 0; l < numPointExt; l++)
    {
        // Calcolo peso nodo GH triangolo esterno
        double pesoExt = stdGHw[l] * areeT[blockIdx.x];
        
        //Lettura nodo GH corrente su triangolo standard
        nodoTemp[0] = stdGHnx[l];
        nodoTemp[1] = stdGHny[l];
        nodoTemp[2] = stdGHnz[l];

        //Mappaggio nodo GH corrente su triangolo esterno
        double nodoGHcurr[3];
        nodoGHcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
        nodoGHcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
        nodoGHcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];
    
        //Calcolo coordinate vettore differenza
        double point[3];
        point[0] = nodoGHCcurr[0] - nodoGHcurr[0];
        point[1] = nodoGHCcurr[1] - nodoGHcurr[1];
        point[2] = nodoGHCcurr[2] - nodoGHcurr[2];

        //Calcolo norma vettore differenza
        double pointNorm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
    
        //Inizializzazione variabili
        double istTemp = 0;
        double tempValues[3][3];
        double coeffNucleo[3] = {1, -2, 1};
        double coeffTemp[3] = {-1, 0, 1};

        //Ciclo sui 3 istanti temporali
        for (k = 0; k < 3; k++)
        {
            // Calcolo istante temporale corrente
            istTemp = deltaT * (offsetZ + blockIdx.z + coeffTemp[k]);
            
            //Check necessità di calcolo
            if(istTemp < 0)
                continue;
            
            //Calcolo delle 9 componenti del nucleo
            nucleo(tempValues, point, pointNorm, istTemp, velP, velS);
            
            //Somma pesata dei valori del nucleo alla shared memory
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += pesoExt * pesoInt * coeffNucleo[k] * tempValues[i][j];
    
        }
    }

    //Sync prima di inziare la riduzione
    __syncthreads();
    
    unsigned int xDim, yDim, sharedOffInd;
    
    //Primo step per scendere a potenza di 2 in y
    yDim = pow(2.0, (int) floor(log2((float) blockDim.y)));
    sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;

    if(threadIdx.y + yDim < blockDim.y)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
    
    __syncthreads();

    //Iterazione per scendere ad 1 in y
    yDim /= 2;
    while(yDim > 0)
    {
        sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;
        if(threadIdx.y < yDim)
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        yDim /= 2;
    }

    //Iterazione per scendere ad 1 in x (x è sempre potenza di 2)
    xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;
        if(threadIdx.x < xDim)
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    //Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0 && threadIdx.y == 0)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                // matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiRho;
            }

    
}

__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp)
{
    int i, j;

    double vertsS[3][3];
    double vertsF[3][3];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];
            vertsF[i][j] = vertsT[9*blockIdx.y + 3*j + i];
        }

    double distMax = 0;
    double distCurr;
    
    double vettDist[3];

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            vettDist[0] = vertsS[i][0] - vertsF[j][0];
            vettDist[1] = vertsS[i][1] - vertsF[j][1];
            vettDist[2] = vertsS[i][2] - vertsF[j][2];
            distCurr = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]);
            if (distCurr > distMax)
                distMax = distCurr;
        }
    
    return ((indTemp - 1) * cS * deltaT > distMax);
}


__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, const double cP, const double cS)
{
    nu[0][0] = ((x[0] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[0] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[0] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[0][1] = ((x[0] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[0] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[0] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[0][2] = ((x[0] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[0] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[0] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));
    nu[1][0] = ((x[1] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[1] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[1] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[1][1] = ((x[1] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[1] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[1] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[1][2] = ((x[1] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[1] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[1] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][0] = ((x[2] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[2] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[2] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][1] = ((x[2] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[2] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[2] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][2] = ((x[2] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[2] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[2] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));


}