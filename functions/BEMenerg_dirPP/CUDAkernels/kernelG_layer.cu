__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, double const cP, double const cS);
__device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double cS, const double diffTempMin, const double diffTempMax, const unsigned int indBlock);

__global__ void kernelPPnewGlayer(double *uRawG, const double velP, const double velS, const double const4PiRho,
                                 const double *sourcePoints, const double *diffTemp, const double *vettLin,
                                 const int numSubRegion, const int numSubPoint, const double *stdPPw, const double *stdPPnx, const double *stdPPny, const double *stdPPnz,
                                 const int numTriang, const double *vertsT, const double *areeT, const double *centerT, const double *maxLenT)  
{
    extern __shared__ double uGshared[][3];
    unsigned int i, j, l, k, p, n;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        uGshared[threadIdx.x][i] = 0;

    //Lettura punto sorgente
    const double sourcePoint[3] = {sourcePoints[3*blockIdx.x + 0], sourcePoints[3*blockIdx.x + 1], sourcePoints[3*blockIdx.x + 2]};
    
    //Ciclo sui blocchi
    for(p = 0; p < numTriang; p++)
    {
        //Check condizione teorica blocco nullo
        if (isBlockNull(centerT, maxLenT, sourcePoint, velP, velS, diffTemp[threadIdx.x + 1], diffTemp[threadIdx.x], p))
            continue;

        //Inizializzazione blocco 3x3
        double matrixBlock[3][3] = {0};
    
        // Lettura vertici triangolo corrente
        double vertsTemp[3][3];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                vertsTemp[i][j] = vertsT[9*p + 3*j + i];
    
        // Inizializzazione istanti temporali
        const int coeffTemp[2] = {-1, 0};
        const int coeffNucleo[2] = {1, -1};
    
        //Ciclo sugli istanti temporali
        for(k = 0; k < 2 ; k++)
        {
            //Lettura istante temporale
            double currT;
            currT = diffTemp[threadIdx.x + 1 + coeffTemp[k]];
    
            if(currT <= 0)
                continue;

            //Ciclo sulle numSubRegion sottoregioni
            for(n = 0; n < numSubRegion; n++)
            {
                //Ciclo sui numSubPoint nodi di ciascuna sottoregione
                for(l = 0; l < numSubPoint; l++)
                {
                    // Lettura coordinate nodo GHC corrente su triangolo stardard
                    double nodoTemp[3];
                    nodoTemp[0] = stdPPnx[n*numSubPoint + l];
                    nodoTemp[1] = stdPPny[n*numSubPoint + l];
                    nodoTemp[2] = stdPPnz[n*numSubPoint + l];
            
                    //Mappaggio nodo GHC corrente su triangolo interno
                    double nodoPPcurr[3];
                    nodoPPcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
                    nodoPPcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
                    nodoPPcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];
                    
                    //Calcolo peso nodo GHC corrente
                    double pesoPP = stdPPw[l] * areeT[p];
            
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
                            matrixBlock[i][j] += coeffNucleo[k] * pesoPP * tempValues[i][j];
                }
            }
        }
    
        //Applicazione coefficiente costante e pulizia numerica
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixBlock[i][j] = (abs(matrixBlock[i][j]) > pow(10.0, -14)) * matrixBlock[i][j] / const4PiRho;
    
        //Lettura componenti della densità/del dato relativi al blocco corrente
        double vettLinCurr[3];
        vettLinCurr[0] = vettLin[(3*numTriang)*threadIdx.x + 3*p + 0];
        vettLinCurr[1] = vettLin[(3*numTriang)*threadIdx.x + 3*p + 1];
        vettLinCurr[2] = vettLin[(3*numTriang)*threadIdx.x + 3*p + 2];
    
        //Calcolo prodotto righe colonne blocco matriciale - blocco densità/dato
        uGshared[threadIdx.x][0] += matrixBlock[0][0] * vettLinCurr[0] + matrixBlock[0][1] * vettLinCurr[1] + matrixBlock[0][2] * vettLinCurr[2];
        uGshared[threadIdx.x][1] += matrixBlock[1][0] * vettLinCurr[0] + matrixBlock[1][1] * vettLinCurr[1] + matrixBlock[1][2] * vettLinCurr[2];
        uGshared[threadIdx.x][2] += matrixBlock[2][0] * vettLinCurr[0] + matrixBlock[2][1] * vettLinCurr[1] + matrixBlock[2][2] * vettLinCurr[2];
    }

    //Sync prima di inziare la riduzione
    __syncthreads();

    //Iterazione per scendere a potenza di 2
    unsigned int xDim;
    xDim = pow(2.0, (int) floor(log2((float) blockDim.x)));
    if(threadIdx.x + xDim < blockDim.x)
        for (i = 0; i < 3; i++)
            uGshared[threadIdx.x][i] += uGshared[threadIdx.x + xDim][i];

    __syncthreads();

    //Iterazioni per scendere ad 1
    xDim /= 2;
    while(xDim > 0)
    {
        if(threadIdx.x < xDim)
            for (i = 0; i < 3; i++)
                uGshared[threadIdx.x][i] += uGshared[threadIdx.x + xDim][i];

        __syncthreads();
        xDim /= 2;
    }

    //Salvataggio dati in memoria globale
    if(threadIdx.x == 0)
        for (i = 0; i < 3; i++)
            uRawG[3*blockIdx.x + i] += uGshared[threadIdx.x][i];
}




__device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double cS, const double diffTempMin, const double diffTempMax, const unsigned int indBlock)
{
    unsigned short i;

    double vettDist[3];
    for (i = 0; i < 3; i++)
        vettDist[i] = sourcePoint[i] - centerT[3*indBlock + i];

    const double distMin = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) - maxLenT[indBlock];
    const double distMax = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) + maxLenT[indBlock];
       
    return (diffTempMax * cP < distMin) || (diffTempMin * cS > distMax);
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