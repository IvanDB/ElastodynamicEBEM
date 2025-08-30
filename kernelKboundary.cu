__device__ void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j);
__device__ void nucleoKR(double nuKR[3][3], const double tau[3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3]);

__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen);
__device__ double dotProd3D(const double vettA[3], const double vettB[3]);

__global__ void kernelKboundary(double *matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double const4PiDeltaT,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const double *normT, const int offsetZ, const int numBlocks, const double maxLen) 
{
    //Evito i blocchi singolari (calcolati su CPU)
    if(blockIdx.x == blockIdx.y)
       return;

    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Controllo condizione teorica necessità di calcolo del singolo blocco
    if (isBlockNull(vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
       return;

    extern __shared__ double matrixSubBlock[][3][3];
    int k, i, j, l, m;

    const unsigned int sharedBaseInd = threadIdx.x;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Dichiarazione e inizializzazione tensore di permutazione E
    const int tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

    //Dichiarazione variabili
    double vertsExt[3][3];
    double vertsInt[3][3];
    double nodoTempExt[3];
    double vettLatiInt[3][3];
    double normLatiInt[3];

    //Lettura vertici triangolo esterno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsExt[i][j] = vertsT[9*blockIdx.x + 3*j + i];  

    // Lettura vertici triangolo interno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsInt[i][j] = vertsT[9*blockIdx.y + 3*j + i];

    //Calcolo componenti vettori lati
    for (j = 0; j < 3; j++)
    {
        vettLatiInt[0][j] = vertsInt[1][j] - vertsInt[0][j];
        vettLatiInt[1][j] = vertsInt[2][j] - vertsInt[1][j];
        vettLatiInt[2][j] = vertsInt[0][j] - vertsInt[2][j];
    }

    //Calcolo norme vettori lati
    for (i = 0; i < 3; i++)
        normLatiInt[i] = sqrt(vettLatiInt[i][0]*vettLatiInt[i][0] + vettLatiInt[i][1]*vettLatiInt[i][1] + vettLatiInt[i][2]*vettLatiInt[i][2]);

    //Calcolo fattori punti thread corrente
    const double startFactor = (double) threadIdx.x / blockDim.x;
    const double endFactor = (double) (threadIdx.x + 1) / blockDim.x;

    //Ciclo sui numPointExt nodi del triangolo esterno
    for(l = 0; l < numPointExt; l++)
    {
        // Calcolo peso nodo GH triangolo esterno
        double pesoExt = stdGHw[l] * areeT[blockIdx.x];

        //Lettura nodo GH corrente su triangolo standard
        nodoTempExt[0] = stdGHnx[l];
        nodoTempExt[1] = stdGHny[l];
        nodoTempExt[2] = stdGHnz[l];

        //Mappaggio nodo GH corrente su triangolo esterno
        double nodoGHcurr[3];
        nodoGHcurr[0] = nodoTempExt[0] * vertsExt[0][0] + nodoTempExt[1] * vertsExt[1][0] + nodoTempExt[2] * vertsExt[2][0];
        nodoGHcurr[1] = nodoTempExt[0] * vertsExt[0][1] + nodoTempExt[1] * vertsExt[1][1] + nodoTempExt[2] * vertsExt[2][1];
        nodoGHcurr[2] = nodoTempExt[0] * vertsExt[0][2] + nodoTempExt[1] * vertsExt[1][2] + nodoTempExt[2] * vertsExt[2][2];

        //Ciclo sui tre lati del triangolo interno
        for (m = 0; m < 3; m++)
        {
            double normSegDec = normLatiInt[m] / blockDim.x;

            double pointStart[3];
            for (i = 0; i < 3; i++)
                pointStart[i] = nodoGHcurr[i] - (vertsInt[m][i] + startFactor * vettLatiInt[m][i]);
                    
            double pointEnd[3];
            for (i = 0; i < 3; i++)
                pointEnd[i] = nodoGHcurr[i] - (vertsInt[m][i] + endFactor * vettLatiInt[m][i]);

            //Calcolo norma vettore differenza
            double normStart = sqrt(pointStart[0]*pointStart[0] + pointStart[1]*pointStart[1] + pointStart[2]*pointStart[2]); 
            double normEnd = sqrt(pointEnd[0]*pointEnd[0] + pointEnd[1]*pointEnd[1] + pointEnd[2]*pointEnd[2]); 
            
            //Riscalamento vettore tangente
            double vettTangCurr[3];
            for (i = 0; i < 3; i++)
                vettTangCurr[i] = vettLatiInt[m][i] / normLatiInt[m];

            //Inizializzazione variabili
            double istTemp = 0;
            double tempValues[3][3];
            const int coeffNucleo[4] = {-1, 3, -3, 1};
            const int coeffTemp[4] = {-2, -1, 0, 1};
    
            //Ciclo sui 4 istanti temporali
            for (k = 0; k < 4; k++)
            { 
                // Calcolo istante temporale corrente
                istTemp = deltaT * (offsetZ + double(blockIdx.z) + coeffTemp[k]);

                //Check necessità di calcolo
                if(istTemp <= 0)
                   continue;

                //Inizializzazione componenti nucleo totale
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        tempValues[i][j] = 0;
                
                //Aggiunta delle delle componenti dei nuclei a Start ed End (formula dei trapezi)
                nucleoKR(tempValues, vettTangCurr, pointStart, normStart, istTemp, velP, velS, lambda, mu, rho, tensorE);
                nucleoKR(tempValues, vettTangCurr, pointEnd, normEnd, istTemp, velP, velS, lambda, mu, rho, tensorE);

                //Somma pesata dei valori del nucleo alla shared memory
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        matrixSubBlock[sharedBaseInd][i][j] += pesoExt * coeffNucleo[k] * (tempValues[i][j] * normSegDec / 2); 
            }
        }
    }

    //Sync prima di inziare la riduzione
    __syncthreads();

    unsigned int xDim, sharedOffInd;

    //Iterazione per scendere ad 1 in x (x supposto potenza di 2)
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

    // Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                //matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiDeltaT;
            }
}

__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
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
    double distMin = 100; //Mettere valore vicino a +Inf;
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
            if (distCurr < distMin)
                distMin = distCurr;
        }
    
    return ((indTemp - 2) * cS * deltaT > distMax) || ((indTemp + 1) * cP * deltaT < distMin - maxLen);
}

__device__ double dotProd3D(const double vettA[3], const double vettB[3])
{
    return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];
}

__device__ void nucleoKR(double nuKR[3][3], const double tau[3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3])
{
    int j, i, k;
    double nucleoCurr[3][3];

    // Inizializzazione nucleo temporaneo
    for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++) 
            nucleoCurr[i][k] = 0;

    // Ciclo sull'indice j dei nuclei nuKRj
    for(j = 0; j < 3; j++)
    {
        // Calcolo nucleo corrente
        nucleoKRj(nucleoCurr, x, r, t, cP, cS, lambda, mu, rho, tensorE, j);

        // Ciclo sulle 9 componenti
        for (i = 0; i < 3; i++)
            for (k = 0; k < 3; k++)
                nuKR[i][k] += nucleoCurr[i][k] * tau[j];
    }
}

__device__ void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j)
{
    nuKRj[0][0] = (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[0] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[0] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP)))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][0][2] * (3 * x[2] * x[1] / pow(r, 5)))
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[1] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[2] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[2] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[2] / pow(r, 3)))
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][1][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][1] = (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[1] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][2][1] * (3 * x[1] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][2] = (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[2] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));
}