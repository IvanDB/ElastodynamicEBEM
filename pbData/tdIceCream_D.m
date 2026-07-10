function g = ConeBall_D(pbParam)

arguments
    pbParam (1, 1) struct
end

R = 0.5;
g = @(x, t) [0; 0; sin(t * (1 - ((x(1)^2 + x(2)^2) / R))).^5];

end