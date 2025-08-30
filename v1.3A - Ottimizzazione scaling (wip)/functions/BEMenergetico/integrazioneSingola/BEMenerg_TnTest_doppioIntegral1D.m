function ris = BEMenerg_TnTest_doppioIntegral1D(nucl_1, nucl_2, nucl_3, L)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione VETTORE contenente VALUTAZIONE PUNTUALE del DATO al BORDO
ris = zeros(3, 1);

%Calcolo valori mediante doppio integral 1D
ris(1) = integral(@(y1) integral(@(y2) nucl_1(y1, y2), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);
ris(2) = integral(@(y1) integral(@(y2) nucl_2(y1, y2), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);
ris(3) = integral(@(y1) integral(@(y2) nucl_3(y1, y2), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5), -L, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);

return