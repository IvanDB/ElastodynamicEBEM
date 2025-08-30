__device__ void nucleo(double nu[3][3], const double x[3], const double r, const double t, double const cP, double const cS);


__global__ void kernelDQC3(double *matrix, const double deltaT, double const velP, double const velS, double const const4PiRho,
                           const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz,
                           const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz, const int numSubPoint,
                           const double *vertsT, const double *areeT)  
{
    //Evito i blocchi diagonali
    if(blockIdx.x == blockIdx.y)
        return;

    extern __shared__ double matrixSubBlock[][3][3];
    int k, i, j, l;

    unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;
    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Dichiarazione variabili
    double vertsTemp[3][3];
    double nodoTemp[3];

    // Calcolo peso nodo GH triangolo esterno
    double pesoExt = stdGHw[threadIdx.x] * areeT[blockIdx.x] * 2;

    // Lettura vertici triangolo esterno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTemp[i][j] = vertsT[9*blockIdx.x + 3*j + i];
    
    //Lettura nodo GH corrente su triangolo standard
    nodoTemp[0] = stdGHnx[threadIdx.x];
    nodoTemp[1] = stdGHny[threadIdx.x];
    nodoTemp[2] = stdGHnz[threadIdx.x];

    //Mappaggio nodo GH corrente su triangolo esterno
    double nodoGHcurr[3];
    nodoGHcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
    nodoGHcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
    nodoGHcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];

    //Lettura vertici triangolo interno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTemp[i][j] = vertsT[9*blockIdx.y + 3*j + i];

    //Ciclo sui numSubPoint nodi di ciascuna sottoregione del triangolo interno
    for(l = 0; l < numSubPoint; l++)
    {
        //Calcolo peso nodo GHC corrente
        double pesoInt = stdGHCw[l] * areeT[blockIdx.y] / (sqrt(3.0)/2);
   
        
        // Lettura coordinate nodo GHC corrente su triangolo stardard
        nodoTemp[0] = stdGHCnx[3*threadIdx.y + l];
        nodoTemp[1] = stdGHCny[3*threadIdx.y + l];
        nodoTemp[2] = stdGHCnz[3*threadIdx.y + l];

        //Mappaggio nodo GHC corrente su triangolo interno
        double nodoGHCcurr[3];
        nodoGHCcurr[0] = nodoTemp[0] * vertsTemp[0][0] + nodoTemp[1] * vertsTemp[1][0] + nodoTemp[2] * vertsTemp[2][0];
        nodoGHCcurr[1] = nodoTemp[0] * vertsTemp[0][1] + nodoTemp[1] * vertsTemp[1][1] + nodoTemp[2] * vertsTemp[2][1];
        nodoGHCcurr[2] = nodoTemp[0] * vertsTemp[0][2] + nodoTemp[1] * vertsTemp[1][2] + nodoTemp[2] * vertsTemp[2][2];
    
    
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
            istTemp = deltaT * (blockIdx.z + coeffTemp[k]);
            
            //Check necessità di calcolo
            if(istTemp < 0)
                continue;
            
            //Calcolo delle 9 componenti del nucleo
            nucleo(tempValues, point, pointNorm, istTemp, velP, velS);
            
            //Somma pesata dei valori dle nucleo alla shared memory
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += pesoExt * pesoInt * coeffNucleo[k] * tempValues[i][j];
    
        }
    }

    //Sync prima di inziare la riduzione
    __syncthreads();

    //Primo step per scendere a potenza di 2 in x
    unsigned int xDim, yDim, sharedOffInd;
    xDim = pow(2.0, (int) floor(log2((float) blockDim.x)));
    sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;

    if(threadIdx.x + xDim < blockDim.x)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];
    
    __syncthreads();

    //Iterazione per scendere ad 1 in x
    xDim /= 2;
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

    //Iterazione per scendere ad 1 in y (y è già potenza di 2)
    yDim = blockDim.y/2;
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

    //Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0 && threadIdx.y == 0)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                // matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                matrix[ind] += matrixSubBlock[0][i][j] / const4PiRho;
            }

    
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