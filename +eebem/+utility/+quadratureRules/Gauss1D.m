function [nodes, weight] = Gauss1D(type, n, alpha, beta)
%GAUSS1D  Compute 1D Gaussian quadrature nodes and weights via the Golub-Welsch algorithm.
%   [NODES, WEIGHT] = GAUSS1D(TYPE, N, ALPHA, BETA) returns N nodes and weights
%   for one of six classical Gaussian quadrature families (selected by TYPE), by
%   building the corresponding three-term recurrence coefficients and computing
%   them as the eigenvalues/first eigenvector-components of the associated
%   Jacobi (tridiagonal) matrix, via the local helper EXALGOLPROCEDURE.
%
%   Input arguments:
%       TYPE  - (positive integer) 1 Gauss-Legendre on (-1,1), w(x)=1; 2
%               Gauss-Chebyshev 1st kind on (-1,1), w(x)=1/sqrt(1-x^2); 3
%               Gauss-Chebyshev 2nd kind on (-1,1), w(x)=sqrt(1-x^2); 4 Gauss-Hermite
%               on (-inf,inf), w(x)=exp(-x^2); 5 Gauss-Jacobi on (-1,1) (unimplemented,
%               see Notes); 6 Gauss-Laguerre on (0,inf), w(x)=exp(-x)*x^alpha.
%       N     - (positive integer) number of nodes/weights requested.
%       ALPHA - (double, optional, default 0) Jacobi/Laguerre exponent.
%       BETA  - (double, optional, default 0) Jacobi exponent (TYPE 5 only).
%
%   Output arguments:
%       NODES  - (1xN double) quadrature nodes.
%       WEIGHT - (1xN double) quadrature weights.
%
%   Notes:
%       TYPE = 5 (Gauss-Jacobi) is not implemented and errors
%       unconditionally ("Caso incompleto"). This project only uses TYPE =
%       1 (Gauss-Legendre), via GENERATEFINALG2DNODES and DOPPIOGAUSS1D.
%
%   See also GENERATEFINALG2DNODES, DOPPIOGAUSS1D

arguments
    type    (1, 1) double {mustBeInteger, mustBePositive}
    n       (1, 1) double {mustBeInteger, mustBePositive}
    alpha   (1, 1) double {} = 0
    beta    (1, 1) double {} = 0
end

arguments (Output)
    nodes  (1, :) double
    weight (1, :) double
end

import eebem.utility.quadratureRules.*

% CALCOLO COEFFICIENTI RELAZIONE A TRE TERMINI E 0-th MOMENT (?)
switch type
    case 1 %Gauss - Legendre -> I = (-1, 1) 
           %                   w(x) = 1
        mu = 2;
        a = zeros(1, n);
        tmp = 1 : (n-1);
        b = tmp ./ sqrt(4.*(tmp.^2) - 1);

    case 2 %Gauss - Chebyshev 1° tipo   -> I = (-1, 1) 
           %                              w(x) = 1 / sqrt(1 - x^2)
        mu = pi;
        a = zeros(1, n);
        b = [sqrt(0.5), 0.5 * ones(1, n-2)];

    case 3 %Gauss - Chebyshev 2° tipo   -> I = (-1, 1) 
           %                              w(x) = sqrt(1 - x^2)
        mu = pi / 2;
        a = zeros(1, n);
        b = 0.5 * ones(1, n-1);
    
    case 4 %Gauss - Hermite -> I = (-inf, +inf) 
           %                  w(x) = exp(-x^2)
        mu = sqrt(pi);
        a = zeros(1, n);
        b = sqrt((1 : (n-1)) ./2 );

    case 5 %Gauss - Jacobi   -> I = (-1, 1) 
           %                  w(x) = (1 - x)^alpha + (1 + x)^beta  (alpha, beta > -1)
        
        mu = (2 ^ (alpha + beta + 1)) * gamma(alpha + 1) .* gamma(beta + 1) / gamma(alpha + beta + 2); 

        error("Caso incompleto")

    case 6 %Gauss - Laguerre -> I = (0, +inf) 
           %                  w(x) = exp(-x) * x^alpha     (alpha > -1)
        mu = gamma(alpha + 1);
        tmp = 1 : n;
        a = (2 .* tmp) - 1 + alpha;
        tmp = 1 : (n-1);
        b = sqrt(tmp .* (tmp + alpha));

    otherwise
        error("Caso non implementato")
end

% CALCOLO NODI E PESI MEDIANTE exAlgolProcedure
% NOTA: This subroutine is a translation of an algol procedure,
%       num. math. 12, 377-383(1968) by Martin and Wilkinson,
%       as modified in num. math. 15, 450(1970) by Dubrulle.
%       handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
%       This is a modified version of the 'eispack' routine imtql2.

[nodes, ~, weight, ~] = exAlgolProcedure(n, a, b);
weight = (weight.^2) .* mu;

end

function [d, e, z, ierr] = exAlgolProcedure(n, d, e)
%EXALGOLPROCEDURE  Golub-Welsch eigen-decomposition of the Jacobi tridiagonal matrix for Gaussian quadrature.
% This subroutine is a translation of an algol procedure,
% num. math. 12, 377-383(1968) by Martin and Wilkinson,
% as modified in num. math. 15, 450(1970) by Dubrulle.
% handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
% This is a modified version of the 'eispack' routine imtql2.
% 
% This subroutine finds the eigenvalues and first components of the
% eigenvectors of a symmetric tridiagonal matrix by the implicit ql
% method.
% 
% on input:
% 
%    n is the order of the matrix;
% 
%    d contains the diagonal elements of the input matrix;
% 
%    e contains the subdiagonal elements of the input matrix
%      in its first n-1 positions.  e(n) is arbitrary;
% 
%    z contains the first row of the identity matrix.
% 
% on output:
% 
%    d contains the eigenvalues in ascending order.  if an
%      error exit is made, the eigenvalues are correct but
%      unordered for indices 1, 2, ..., ierr-1;
% 
%    e has been destroyed;
% 
%    z contains the first components of the orthonormal eigenvectors
%      of the symmetric tridiagonal matrix.  if an error exit is
%      made, z contains the eigenvectors associated with the stored
%      eigenvalues;
% 
%    ierr is set to
%      zero       for normal return,
%      j          if the j-th eigenvalue has not been
%                 determined after 30 iterations.

ierr = 0;

z = zeros(1, n); 
z(1) = 1;

if(n == 1) 
    return;
end

e(n) = 0;
for l = 1 : n
    j = 0;
    %   :::::::::: look for small sub-diagonal element ::::::::::
    while true
        for m = l : n
            if m == n
                break
            end
            if abs(e(m)) <= eps * (abs(d(m)) + abs(d(m+1)))
                break
            end
        end
        
        p = d(l);

        if m == l
            break;
        end

        if j == 30
            ierr = l;
            return
        end
        j = j + 1;

        %   :::::::::: form shift ::::::::::
        g = (d(l+1) - p) / (2 * e(l));
        r = hypot(g, 1.0);

        ggg = abs(r) * (1 - 2*(g < 0));
       
        g = d(m) - p + e(l) / (g + ggg);
        s = 1;
        c = 1;
        p = 0;
            
        for i = (m - 1) : -1 : l
            f = s * e(i);
            b = c * e(i);

            signFactor = sign(g) * (abs(g) > abs(f)) + sign(f) * (abs(g) <= abs(f));

            r = hypot(f, g);           
            e(i+1) = r * signFactor;
            s = f / r * signFactor;
            c = g / r * signFactor;

            g = d(i+1) - p;
            r = (d(i) - g) * s + 2.0 * c * b;
            p = s * r;
            d(i+1) = g + p;
            g = c * r - b;
            
            %   :::::::::: form first component of vector ::::::::::
            f = z(i+1);
            z(i+1) = s * z(i) + c * f;
            z(i) = c * z(i) - s * f;
        end

        d(l) = d(l) - p;
        e(l) = g;
        e(m) = 0.0;
    end
end 

% Order eigenvalues and eigenvectors
[d, idx] = sort(d);
z = z(idx);

end
