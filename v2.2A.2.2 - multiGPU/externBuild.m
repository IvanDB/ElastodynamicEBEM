function externBuild()

% CPU code
compilerInfo = mex.getCompilerConfigurations('C');
while(isempty(compilerInfo))
    mex -setup C
    compilerInfo = mex.getCompilerConfigurations('C');
end

cmd = "mex -v -R2018a COMPFLAGS='$COMPFLAGS /Ox' kernelV_test.c";
cmdout = evalc(cmd);
if(~contains(cmdout, "MEX completed successfully."))
    error("C compilation failed! Shell output:" + newline + cmdout)
end

cmd = "mex -v -R2018a COMPFLAGS='$COMPFLAGS /Ox' kernelK_test.c";
cmdout = evalc(cmd);
if(~contains(cmdout, "MEX completed successfully."))
    error("C compilation failed! Shell output:" + newline + cmdout)
end

% GPU code
scriptMSCV = ['"' compilerInfo.Details.CommandLineShell '" ' ...
                    compilerInfo.Details.CommandLineShellArg];
gpuInfo = gpuDevice;
gpuCC = gpuInfo.ComputeCapability;
gpuCC(gpuCC == '.') = [];

[status, cmdout] = system([scriptMSCV ' && nvcc -ptx -O3 -arch=compute_' gpuCC ' ./functions/BEMenerg_dir/CUDAkernels/kernelK.cu -Wno-deprecated-gpu-targets']);
if(status)
    error("CUDA compilation failed! Shell output:" + newline + cmdout)
end

[status, cmdout] = system([scriptMSCV ' && nvcc -ptx -O3 -arch=compute_' gpuCC ' ./functions/BEMenerg_dir/CUDAkernels/kernelV.cu -Wno-deprecated-gpu-targets']);
if(status)
    error("CUDA compilation failed! Shell output:" + newline + cmdout)
end
disp("Successful compilation")
end