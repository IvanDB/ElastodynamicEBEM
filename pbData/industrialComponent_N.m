function g = industrialComponent_N(pbParam)
%industrialComponent_N undefined
%   undefined

arguments
    pbParam (1, 1) struct
end

g = @(x, t, n) [0; 0; sign(n(3)) * (abs(n(3)) > 0.5)];
end
