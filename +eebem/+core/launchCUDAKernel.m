function deviceMatrix = launchCUDAKernel(gpuID, kernel, varargin)
%LAUNCHCUDAKERNEL Summary of this function goes here
%   Detailed explanation goes here
arguments
    gpuID (1,1) parallel.gpu.GPUDevice
    kernel (1,1) parallel.gpu.CUDAKernel
end

arguments (Repeating)
    varargin
end

%Array di input
linearMatrixSixe = 9 * prod(kernel.GridSize);

deviceMatrix = zeros(linearMatrixSixe, 1, 'gpuArray');
wait(gpuID);

%Avvio computazione GPU
deviceMatrix = feval(kernel, deviceMatrix, varargin{:});

end