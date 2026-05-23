clc; 
k = 4; p = 3;
%%
for N = 1 : 26
    clearvars -except N p k
    nIter_ful = floor(N / (k * p));
    nIter_tot = ceil(N / (k * p));
    
    nBlks_fin = mod(N, (k * p));
    
    blkSize_I = ceil(nBlks_fin / k);
    blkSize_F = floor(nBlks_fin / k);
    
    nBlocks_F = mod(-nBlks_fin, k);
    nBlocks_I = (k - nBlocks_F) * (nBlks_fin ~= 0);
    
    offsets = 0 : p : (k*p)*nIter_ful;
    
    offsets_b = blkSize_I * (1 : nBlocks_I);
    offsets_bf = offsets(end) + offsets_b;
    offsets = [offsets, offsets_bf];

    offsets_e = blkSize_F * (1 : nBlocks_F);
    offsets_ef = offsets(end) + offsets_e;
    
    offsets = [offsets, offsets_ef];
    disp("" + sprintf('%2d | ', N) + sprintf('%2d %2d %2d %2d - ', offsets))
end


 





