__device__ void nucleo(double nu[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
                       const double cP, const double cS, const double lambda, const double mu);
__device__ bool isSubBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double cS, const double diffTempMin, const double diffTempMax, const unsigned int indTriang);

__device__ double dotProd3D(const double vettA[3], const double vettB[3]);
__device__ double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM);

__global__ void kernelPPnewPlayer(double *uRawP, const double velP, const double velS, const double lambda, const double mu, const double const4PiRhoDeltaT,
                         const double *sourcePoints, const double *diffTemp, const double *vettLin,
                         const int numSubRegion, const int numSubPoint, const double *stdPPw, const double *stdPPnx, const double *stdPPny, const double *stdPPnz,
                         const int numTriang, const double *vertsT, const double *areeT, const double *normT, const double *centerT, const double *maxLenT, 
                         const int numNodes, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff)  
{
    extern __shared__ double uPshared[][3];
    unsigned int i, j, l, k, m, n, p;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        uPshared[threadIdx.x][i] = 0;

    //Lettura punto sorgente
    const double sourcePoint[3] = {sourcePoints[3*blockIdx.x + 0], sourcePoints[3*blockIdx.x + 1], sourcePoints[3*blockIdx.x + 2]};
    
    //Ciclo sui blocchi
    for(p = 0; p < numNodes; p++)
    {
        //Inizializzazione blocco 3x3
        double matrixBlock[3][3] = {0};

        //Inizializzazione indice matrice delle flag
        const int baseIndFlag = p * numTriang;
        
        //Ciclo sui triangoli di campo
        for(m = 0; m < numTriang; m++)
        {
            
            //Estrazione indice vertice-triangolo corrente
            int indSMcurr = indSMmatrix[baseIndFlag + m];
    
            //Check indice corrente
            if(indSMcurr == 0)
                continue;

            //Check condizione teorica sottoblocco nullo
            if (isSubBlockNull(centerT, maxLenT, sourcePoint, velP, velS, diffTemp[threadIdx.x + 2], diffTemp[threadIdx.x], m))
                continue;
    
            // Lettura vertici triangolo corrente
            double vertsTemp[3][3];
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    vertsTemp[i][j] = vertsT[9*m + 3*j + i];

            // Lettura vettore normale corrente
            double normIntCurr[3];
            for (i = 0; i < 3; i++)
                normIntCurr[i] = normT[3*m + i];
    
            //Estrazione matrice e vettore dei coefficienti correnti
            double matCoeffCurr[3][3];
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                   matCoeffCurr[i][j] = matCoeff[9*m + 3*j + i];
    
            double vetCoeffCurr[3];
            for (i = 0; i < 3; i++)
                vetCoeffCurr[i] = vetCoeff[3*m + i];
            
            // Inizializzazione istanti temporali
            const int coeffTemp[3] = {-1, 0, 1};
            const int coeffNucleo[3] = {1, -2, 1};

            //Ciclo sugli istanti temporali
            for(k = 0; k < 3 ; k++)
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
                        double pesoPP = stdPPw[l] * areeT[m];
            
                        //Calcolo coordinate vettore differenza
                        double point[3];
                        point[0] = sourcePoint[0] - nodoPPcurr[0];
                        point[1] = sourcePoint[1] - nodoPPcurr[1];
                        point[2] = sourcePoint[2] - nodoPPcurr[2];
            
                        //Calcolo norma vettore differenza
                        double pointNorm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
                
                        //Matrice contenente il risultato del singolo nodo
                        double tempValues[3][3];
                    
                        //Calcolo delle 9 componenti del nucleo
                        nucleo(tempValues, point, normIntCurr, pointNorm, currT, nodoPPcurr, matCoeffCurr, vetCoeffCurr, indSMcurr, velP, velS, lambda, mu);
            
                        //Somma pesata dei valori del nucleo alla shared memory
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 3; j++)
                                matrixBlock[i][j] += coeffNucleo[k] * pesoPP * tempValues[i][j];
                    }
                }
            }
        }

        //Applicazione coefficiente costante e pulizia numerica
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixBlock[i][j] = (abs(matrixBlock[i][j]) > pow(10.0, -14)) * matrixBlock[i][j] / const4PiRhoDeltaT;
    
        //Lettura componenti della densità/del dato relativi al blocco corrente
        double vettLinCurr[3];
        vettLinCurr[0] = vettLin[(3*numNodes)*threadIdx.x + 3*p + 0];
        vettLinCurr[1] = vettLin[(3*numNodes)*threadIdx.x + 3*p + 1];
        vettLinCurr[2] = vettLin[(3*numNodes)*threadIdx.x + 3*p + 2];
    
        //Calcolo prodotto righe colonne blocco matriciale - blocco densitàdato
        uPshared[threadIdx.x][0] += matrixBlock[0][0] * vettLinCurr[0] + matrixBlock[0][1] * vettLinCurr[1] + matrixBlock[0][2] * vettLinCurr[2];
        uPshared[threadIdx.x][1] += matrixBlock[1][0] * vettLinCurr[0] + matrixBlock[1][1] * vettLinCurr[1] + matrixBlock[1][2] * vettLinCurr[2];
        uPshared[threadIdx.x][2] += matrixBlock[2][0] * vettLinCurr[0] + matrixBlock[2][1] * vettLinCurr[1] + matrixBlock[2][2] * vettLinCurr[2];
    }
    
    //Sync prima di inziare la riduzione
    __syncthreads();

    //Iterazione per scendere a potenza di 2
    unsigned int xDim;
    xDim = pow(2.0, (int) floor(log2((float) blockDim.x)));
    if(threadIdx.x + xDim < blockDim.x)
        for (i = 0; i < 3; i++)
            uPshared[threadIdx.x][i] += uPshared[threadIdx.x + xDim][i];

    __syncthreads();

    //Iterazioni per scendere ad 1
    xDim /= 2;
    while(xDim > 0)
    {
        if(threadIdx.x < xDim)
            for (i = 0; i < 3; i++)
                uPshared[threadIdx.x][i] += uPshared[threadIdx.x + xDim][i];

        __syncthreads();
        xDim /= 2;
    }

    //Salvataggio dati in memoria globale
    if(threadIdx.x == 0)
        for (i = 0; i < 3; i++)
            uRawP[3*blockIdx.x + i] = uPshared[threadIdx.x][i];
}

__device__ bool isSubBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double cS, const double diffTempMin, const double diffTempMax, const unsigned int indTriang)
{
    unsigned short i;

    double vettDist[3];
    for (i = 0; i < 3; i++)
        vettDist[i] = sourcePoint[i] - centerT[3*indTriang + i];

    const double distMin = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) - maxLenT[indTriang];
    const double distMax = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) + maxLenT[indTriang];
       
    return (diffTempMax * cP < distMin) || (diffTempMin * cS > distMax);
}

__device__ double dotProd3D(const double vettA[3], const double vettB[3])
{
    return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];
}

__device__ double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM)
{
    double leftVector[3];
    leftVector[0] = vetCoeff[0] + matCoeff[0][0] * nodeInt[0] + matCoeff[0][1] * nodeInt[1] + matCoeff[0][2] * nodeInt[2];
    leftVector[1] = vetCoeff[1] + matCoeff[1][0] * nodeInt[0] + matCoeff[1][1] * nodeInt[1] + matCoeff[1][2] * nodeInt[2];
    leftVector[2] = vetCoeff[2] + matCoeff[2][0] * nodeInt[0] + matCoeff[2][1] * nodeInt[1] + matCoeff[2][2] * nodeInt[2];

    const double rightVector[3] = {indSM == 1, indSM == 2, indSM == 3};
    return dotProd3D(leftVector, rightVector);
}


__device__ void nucleo(double nu[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
                       const double cP, const double cS, const double lambda, const double mu)
{
    nu[0][0] = ((lambda * x[0] * n[0] / pow(r, 2)) + (2 * mu * x[0] * x[0] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[0] * n[0] / pow(r, 2)) + (((1 / pow(r, 2)) - (2 * x[0] * x[0] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[0] * n[0] / pow(r, 3)) - (2 * mu * x[0] * n[0] / pow(r, 3)) - (2 * mu * ((1 / pow(r, 3)) - (6 * x[0] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[0] * n[0] / pow(r, 3)) + (3 * x[0] * n[0] / pow(r, 3)) + (3 * ((1 / pow(r, 3)) - (4 * x[0] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[0] * n[0] / pow(r, 5)) + (x[0] * n[0] / pow(r, 5)) + (((1 / pow(r, 5)) - (5 * x[0] * x[0] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[0][1] = ((lambda * x[0] * n[1] / pow(r, 2)) + (2 * mu * x[1] * x[0] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[1] * n[0] / pow(r, 2)) + (- ((2 * x[1] * x[0] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[0] * n[1] / pow(r, 3)) - (2 * mu * x[1] * n[0] / pow(r, 3)) - (2 * mu * (- (6 * x[1] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[0] * n[1] / pow(r, 3)) + (3 * x[1] * n[0] / pow(r, 3)) + (3 * (- (4 * x[1] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[0] * n[1] / pow(r, 5)) + (x[1] * n[0] / pow(r, 5)) + ((- (5 * x[1] * x[0] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[0][2] = ((lambda * x[0] * n[2] / pow(r, 2)) + (2 * mu * x[2] * x[0] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[2] * n[0] / pow(r, 2)) + (- ((2 * x[2] * x[0] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[0] * n[2] / pow(r, 3)) - (2 * mu * x[2] * n[0] / pow(r, 3)) - (2 * mu * (- (6 * x[2] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[0] * n[2] / pow(r, 3)) + (3 * x[2] * n[0] / pow(r, 3)) + (3 * (- (4 * x[2] * x[0] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[0] * n[2] / pow(r, 5)) + (x[2] * n[0] / pow(r, 5)) + ((- (5 * x[2] * x[0] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[1][0] = ((lambda * x[1] * n[0] / pow(r, 2)) + (2 * mu * x[0] * x[1] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[0] * n[1] / pow(r, 2)) + (- ((2 * x[0] * x[1] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[1] * n[0] / pow(r, 3)) - (2 * mu * x[0] * n[1] / pow(r, 3)) - (2 * mu * (- (6 * x[0] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[1] * n[0] / pow(r, 3)) + (3 * x[0] * n[1] / pow(r, 3)) + (3 * (- (4 * x[0] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[1] * n[0] / pow(r, 5)) + (x[0] * n[1] / pow(r, 5)) + ((- (5 * x[0] * x[1] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[1][1] = ((lambda * x[1] * n[1] / pow(r, 2)) + (2 * mu * x[1] * x[1] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[1] * n[1] / pow(r, 2)) + (((1 / pow(r, 2)) - (2 * x[1] * x[1] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[1] * n[1] / pow(r, 3)) - (2 * mu * x[1] * n[1] / pow(r, 3)) - (2 * mu * ((1 / pow(r, 3)) - (6 * x[1] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[1] * n[1] / pow(r, 3)) + (3 * x[1] * n[1] / pow(r, 3)) + (3 * ((1 / pow(r, 3)) - (4 * x[1] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[1] * n[1] / pow(r, 5)) + (x[1] * n[1] / pow(r, 5)) + (((1 / pow(r, 5)) - (5 * x[1] * x[1] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[1][2] = ((lambda * x[1] * n[2] / pow(r, 2)) + (2 * mu * x[2] * x[1] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[2] * n[1] / pow(r, 2)) + (- ((2 * x[2] * x[1] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[1] * n[2] / pow(r, 3)) - (2 * mu * x[2] * n[1] / pow(r, 3)) - (2 * mu * (- (6 * x[2] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[1] * n[2] / pow(r, 3)) + (3 * x[2] * n[1] / pow(r, 3)) + (3 * (- (4 * x[2] * x[1] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[1] * n[2] / pow(r, 5)) + (x[2] * n[1] / pow(r, 5)) + ((- (5 * x[2] * x[1] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[2][0] = ((lambda * x[2] * n[0] / pow(r, 2)) + (2 * mu * x[0] * x[2] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[0] * n[2] / pow(r, 2)) + (- ((2 * x[0] * x[2] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[2] * n[0] / pow(r, 3)) - (2 * mu * x[0] * n[2] / pow(r, 3)) - (2 * mu * (- (6 * x[0] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[2] * n[0] / pow(r, 3)) + (3 * x[0] * n[2] / pow(r, 3)) + (3 * (- (4 * x[0] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[2] * n[0] / pow(r, 5)) + (x[0] * n[2] / pow(r, 5)) + ((- (5 * x[0] * x[2] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[2][1] = ((lambda * x[2] * n[1] / pow(r, 2)) + (2 * mu * x[1] * x[2] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[1] * n[2] / pow(r, 2)) + (- ((2 * x[1] * x[2] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[2] * n[1] / pow(r, 3)) - (2 * mu * x[1] * n[2] / pow(r, 3)) - (2 * mu * (- (6 * x[1] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[2] * n[1] / pow(r, 3)) + (3 * x[1] * n[2] / pow(r, 3)) + (3 * (- (4 * x[1] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[2] * n[1] / pow(r, 5)) + (x[1] * n[2] / pow(r, 5)) + ((- (5 * x[1] * x[2] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));

    nu[2][2] = ((lambda * x[2] * n[2] / pow(r, 2)) + (2 * mu * x[2] * x[2] / pow(r, 4) * dotProd3D(x, n))) * (((t - (r/cP)) > 0) / pow(cP, 3))
                    + mu * ((x[2] * n[2] / pow(r, 2)) + (((1 / pow(r, 2)) - (2 * x[2] * x[2] / pow(r, 4))) * dotProd3D(x, n))) * (((t - (r/cS)) > 0) / pow(cS, 3))
                    + (((lambda - 2 * mu) * x[2] * n[2] / pow(r, 3)) - (2 * mu * x[2] * n[2] / pow(r, 3)) - (2 * mu * ((1 / pow(r, 3)) - (6 * x[2] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cP)) * ((t - (r/cP)) > 0) / pow(cP, 2))
                    + mu * ((2 * x[2] * n[2] / pow(r, 3)) + (3 * x[2] * n[2] / pow(r, 3)) + (3 * ((1 / pow(r, 3)) - (4 * x[2] * x[2] / pow(r, 5))) * dotProd3D(x, n))) * ((t - (r/cS)) * ((t - (r/cS)) > 0) / pow(cS, 2))
                    - 6 * mu * ((x[2] * n[2] / pow(r, 5)) + (x[2] * n[2] / pow(r, 5)) + (((1 / pow(r, 5)) - (5 * x[2] * x[2] / pow(r, 7))) * dotProd3D(x, n))) 
                                                * (((t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0) - (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0)));
    
    // Calcolo funzione di base
    const double baseFunctionValue = baseFunctionSM(nodeInt, matCoeff, vetCoeff, indSM);

    unsigned int i, j;
    //Applicazione componente relativa alla funzione di base
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) 
            nu[i][j] *= baseFunctionValue;

}