clc
clear
close all
reset(gpuDevice);
addpath(genpath("./functions"))
addpath(genpath("./sorgentiGPU"))

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");

%Definizione 
glbIndexFigures = 10;

%%
input_file_name = 'input_sphereUniform_mid.txt';
numGHext = 12;
numGHCint = 16;
pb_param = readInputFile(input_file_name);
numSlice = pb_param.Nt;

domainMesh = readSpaceMesh(pb_param.domain_type, pb_param.lev);
%glbIndexFigures = plotMesh(domainMesh, glbIndexFigures);


[gha, ghw] = GaussHammer_base(28);
[ghcn, ghcw] = GaussHammerComposite(numGHCint, 3);

matrixIn = gpuArray(zeros((3*domainMesh.number_triangles)^2 * numSlice, 1));
deltaT = pb_param.T_fin / pb_param.Nt;
stdGHw = squeeze(ghw(numGHext, 1:numGHext));
stdGHnx = squeeze(gha(numGHext, 1:numGHext, 1));
stdGHny = squeeze(gha(numGHext, 1:numGHext, 2));
stdGHnz = squeeze(gha(numGHext, 1:numGHext, 3));
stdGHCw = ghcw;
stdGHCnx = squeeze(ghcn(:, 1));
stdGHCny = squeeze(ghcn(:, 2));
stdGHCnz = squeeze(ghcn(:, 3));
areeT = domainMesh.area;
vertsT = zeros(domainMesh.number_triangles*9, 1);
for i = 0 : (domainMesh.number_triangles - 1)
    vertsT(9*i + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(i+1, 1:3), :), [9 1]);
end

rhs = zeros(3*domainMesh.number_triangles, numSlice);
tnf = zeros(3*domainMesh.number_triangles, 1);
density = zeros(3*domainMesh.number_triangles, numSlice);

%% GPU SETUP
%kernel = parallel.gpu.CUDAKernel("kernelTEST.ptx", "kernelTEST.cu");
kernel = parallel.gpu.CUDAKernel("kernelDQC3.ptx", "kernelDQC3.cu");
kernel.GridSize = [domainMesh.number_triangles domainMesh.number_triangles numSlice];
kernel.ThreadBlockSize = [numGHext numGHCint 1];
kernel.SharedMemorySize = numGHext * numGHCint * 9 * 8;

%% RUN
clc
disp("Done - setup")
%pause(2)
%wait(gpuDevice);
disp("Done - sleep")
val = 1;
for temp = 1 : 10
    matrixOut{temp} = feval(kernel, matrixIn, deltaT, pb_param.velP, pb_param.velS, pi*pb_param.rho*val,  ...
                        stdGHw, stdGHnx, stdGHny, stdGHnz,...
                        stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, 3, ...
                        vertsT, areeT);
    val = 100*max(matrixOut{1});
    disp("Done - launch" + num2str(temp))
    pause(0.5)
end
for i = 1 : 0
    parfor j = 1 : domainMesh.number_triangles
        matrixSubBlockSA(:, :, j, i) = BEMenerg_calcMatrxSubBlockDiag_SA(j, pb_param, domainMesh, i-1, deltaT, gha, ghw);
    end
    rhs(:, i) = BEMenerg_calcTnBlock(rhs(:, i), i - 1, deltaT, pb_param, domainMesh, gha, ghw);
    disp(strcat("Done - CPUiter", num2str(i)));
end
disp("Done - CPU")
wait(gpuDevice);
disp("Done - GPU")
%%
matrixOut = reshape(matrixOut, [3*domainMesh.number_triangles 3*domainMesh.number_triangles numSlice]);


for i = 1 : numSlice
    for j = 1 : domainMesh.number_triangles
        indRC = j - 1;
        matrixOut(3*indRC + (1:3), 3*indRC + (1:3), i) = matrixSubBlockSA(:, :, j, i);
    end
end

matrixSist = matrixOut(:, :, 1);
[L, U, P] = lu(matrixSist);



for currInd = 1 : numSlice
    tnf = rhs(:, currInd);

    for ind_inner = 2 : currInd    
        tnf = tnf - matrixOut(:, :, ind_inner) * density(:, currInd - ind_inner + 1);
    end
    
    density(:, currInd) = U\(L\(P*tnf));
end


%%
glbIndexFigures = plotDensity(pb_param, domainMesh, density, glbIndexFigures);
