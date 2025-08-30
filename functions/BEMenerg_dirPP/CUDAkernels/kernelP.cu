__device__ void nucleo(double nu[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
                       const double cP, const double cS, const double lambda, const double mu);
// __device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double currDiffTemp);
__device__ double dotProd3D(const double vettA[3], const double vettB[3]);
__device__ double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM);

__global__ void kernelPPNEWP(double *matrix, const double velP, const double velS, const double lambda, const double mu, const double const4PiRhoDeltaT, const double numT,
                         const double *sourcePoint, const double *diffTemp,
                         const double *stdPPw, const double *stdPPnx, const double *stdPPny, const double *stdPPnz, const int numSubPoint,
                         const double *vertsT, const double *areeT, const double *normT, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff)  
{
    extern __shared__ double matrixSubBlock[][3][3];
    unsigned int i, j, l, k, m;

    //Controllo condizione teorica necessità di calcolo
    //if (isBlockNull(centerT, maxLenT, sourcePoint, velP, diffTemp[blockIdx.y]))
    //    return;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[threadIdx.x][i][j] = 0;

    //Dichiarazione variabili
    double vertsTemp[3][3];

    // Inizializzazione istanti temporali
    const int coeffTemp[3] = {-1, 0, 1};
    const int coeffNucleo[3] = {1, -2, 1};

    //Ciclo sugli istanti temporali
    for(k = 0; k < 3 ; k++)
    {
        //Lettura istante temporale
        double currT;
        currT = diffTemp[blockIdx.y + 1 + coeffTemp[k]];

        if(currT <= 0)
            continue;

        //Inizializzazione indice matrice delle flag
        const int baseIndFlag = blockIdx.x * numT;
    
        //Ciclo sui triangoli di campo
        for(m = 0; m < numT; m++)
        { 
            //Estrazione indice vertice-triangolo corrente
            int indSMcurr = indSMmatrix[baseIndFlag + m];
    
            //Check indice associato al blocco
            if(indSMcurr == 0)
                continue;
    
            //Estrazione vettore normale corrente
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
    
            // Lettura vertici triangolo
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    vertsTemp[i][j] = vertsT[9*m + 3*j + i];
    
            //Ciclo sui numSubPoint nodi di ciascuna sottoregione del triangolo interno
            for(l = 0; l < numSubPoint; l++)
            {
                //Dichiarazione variabili
                double nodoTemp[3];
        
                // Lettura coordinate nodo GHC corrente su triangolo stardard
                nodoTemp[0] = stdPPnx[numSubPoint*threadIdx.x + l];
                nodoTemp[1] = stdPPny[numSubPoint*threadIdx.x + l];
                nodoTemp[2] = stdPPnz[numSubPoint*threadIdx.x + l];
        
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
                        matrixSubBlock[threadIdx.x][i][j] += coeffNucleo[k] * pesoPP * tempValues[i][j];
            }
        }
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
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiRhoDeltaT;
            }    
}

// __device__ bool isBlockNull(const double *centerT, const double *maxLenT, const double *sourcePoint, const double cP, const double currDiffTemp)
// {
//     unsigned short i;
// 
//     double vettDist[3];
//     for (i = 0; i < 3; i++)
//         vettDist[i] = sourcePoint[i] - centerT[3*blockIdx.x + i];
// 
//     const double distMin = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]) - maxLenT[blockIdx.x];
// 
//     return (currDiffTemp * cP < distMin);
// }

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