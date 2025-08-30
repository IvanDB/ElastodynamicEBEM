__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, double const cP, double const cS);
__device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double currDiffTemp);


__global__ void kernelPP(double *matrix, const double velP, const double velS, const double const4PiRho,
                         const double *sourcePoint, const double *diffTemp,
                         const double *stdPPw, const double *stdPPnx, const double *stdPPny, const double *stdPPnz, const int numSubPoint,
                         const double *vertsT, const double *areeT, const double *centerT, const double *maxLenT)  
{
    extern __shared__ double matrixSubBlock[][3][3];
    unsigned int i, j, l;

    //Controllo condizione teorica necessità di calcolo
    if (isBlockNull(centerT, maxLenT, sourcePoint, velP, diffTemp[blockIdx.y]))
        return;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[threadIdx.x][i][j] = 0;

    //Dichiarazione variabili
    double vertsTemp[3][3];

    // Lettura vertici triangolo
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTemp[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    //Lettura istante temporale
    double currT;
    currT = diffTemp[blockIdx.y];

    //Ciclo sui numSubPoint nodi di ciascuna sottoregione del triangolo interno
    for(l = 0; l < numSubPoint; l++)
    {
        //Dichiarazione variabili
        double nodoTemp[3];

        // Lettura coordinate nodo GHC corrente su triangolo stardard
        nodoTemp[0] = stdPPnx[3*threadIdx.x + l];
        nodoTemp[1] = stdPPny[3*threadIdx.x + l];
        nodoTemp[2] = stdPPnz[3*threadIdx.x + l];

        //Mappaggio nodo GHC corrente su triangolo interno
        double nodoPPcurr[3];
        nodoPPcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
        nodoPPcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
        nodoPPcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];
        
        //Calcolo peso nodo GHC corrente
        double pesoPP = stdPPw[l] * areeT[blockIdx.x];

        //Calcolo coordinate vettore differenza
        double point[3];
        point[0] = nodoPPcurr[0] - sourcePoint[0];
        point[1] = nodoPPcurr[1] - sourcePoint[1];
        point[2] = nodoPPcurr[2] - sourcePoint[2];

        //Calcolo norma vettore differenza
        double pointNorm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);

        //Matrice contenente il risultato del singolo nodo
        double tempValues[3][3];
        
        //Calcolo delle 9 componenti del nucleo
        nucleo(tempValues, point, pointNorm, currT, velP, velS);

        //Somma pesata dei valori del nucleo alla shared memory
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixSubBlock[threadIdx.x][i][j] += pesoPP * tempValues[i][j];
    }
    
    //Sync prima di inziare la riduzione
    __syncthreads();

    //Iterazione per scendere ad 1
    unsigned int xDim, sharedOffInd, sharedBaseInd;
    
    sharedBaseInd = threadIdx.x;
    xDim = blockDim.x / 2;
    while(xDim > 0)
    {
        sharedOffInd = threadIdx.x + xDim;
        if(threadIdx.x < xDim)
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    //Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                ind = 3*gridDim.y*(3*blockIdx.x + j) + 3*blockIdx.y + i;
                //matrix[3*blockIdx.y + i][3*blockIdx.x + j]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiRho;
            }    
}




__device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double currDiffTemp)
{
    unsigned short i;

    double vettDist[3];
    for (i = 0; i < 3; i++)
        vettDist[i] = sourcePoint[i] - centerT[3*blockIdx.x + i];

    const double distMin = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) - maxLenT[blockIdx.x];
       
    return (currDiffTemp * cP < distMin);
}




__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, double const cP, double const cS)
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