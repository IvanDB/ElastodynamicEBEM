function ris = BEMenerg_core_TnTest_classic(nucl_1, nucl_2, nucl_3, L)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione VETTORE contenente VALUTAZIONE PUNTUALE del DATO al BORDO
ris = zeros(3,1); 

ris(1) = integral2(nucl_1, -L, L, -L, L, 'Method', 'iterated', 'AbsTol', 1e-6, 'RelTol', 1e-5);
ris(2) = integral2(nucl_2, -L, L, -L, L, 'Method', 'iterated', 'AbsTol', 1e-6, 'RelTol', 1e-5);
ris(3) = integral2(nucl_3, -L, L, -L, L, 'Method', 'iterated', 'AbsTol', 1e-6, 'RelTol', 1e-5);
return