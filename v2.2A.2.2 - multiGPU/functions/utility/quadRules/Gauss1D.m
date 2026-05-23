function [nodes, weight] = Gauss1D(type, n, alpha, beta)
% INPUT 
%   - type: intero indicante la tipologia di nodi di Gauss da generare
%   - n: intero contenente il numero di nodi da generare
%   - alpha: parametro opzionale necessario per i nodi di Gauss-
%   - beta: parametro opzionale necessario per i nodi di Gauss-
%
% OUTPUT:
%   - nodes: vettore contenente le coordinate dei nodi di Gauss-* sul
%                   rispettivo intervallo di riferimento
%   - weight: vettore contenente i pesi di Gauss-* sul rispettivo 
%                   intervallo di riferimento



%% CALCOLO COEFFICIENTI RELAZIONE A TRE TERMINI E 0-th MOMENT (?)

% Vecchia chiamata
% [a, b, mu] = class(type, n, alpha, beta);

% Controllo tipologia nodi richiesta
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

%% CALCOLO NODI E PESI MEDIANTE ??
% NOTA: this subroutine is a translation of an algol procedure,
%       num. math. 12, 377-383(1968) by martin and wilkinson,
%       as modified in num. math. 15, 450(1970) by dubrulle.
%       handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
%       this is a modified version of the 'eispack' routine imtql2.

[nodes, e, weight, ier] = exAlgolProcedure(n, a, b); % -> controllare
weight = (weight.^2) .* mu;

end

