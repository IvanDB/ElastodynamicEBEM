function ris = BEMenerg_TnTest_doppioIntegral1DSpezzato(nucl_1, nucl_2, nucl_3, L, x1, x2)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione VETTORE contenente VALUTAZIONE PUNTUALE del DATO al BORDO
ris = zeros(3, 1);

%Calcolo valori mediante doppio integral 1D
temp_n1_p1 = @(y1) integral(@(y2) nucl_1(y1, y2), -L, x2, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n1_p2 = @(y1) integral(@(y2) nucl_1(y1, y2), x2, L, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n1 = @(y1) temp_n1_p1(y1) + temp_n1_p2(y1);

temp_n2_p1 = @(y1) integral(@(y2) nucl_2(y1, y2), -L, x2, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n2_p2 = @(y1) integral(@(y2) nucl_2(y1, y2), x2, L, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n2 = @(y1) temp_n2_p1(y1) + temp_n2_p2(y1);

temp_n3_p1 = @(y1) integral(@(y2) nucl_3(y1, y2), -L, x2, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n3_p2 = @(y1) integral(@(y2) nucl_3(y1, y2), x2, L, 'AbsTol', 1e-6, 'RelTol', 1e-5);
temp_n3 = @(y1) temp_n3_p1(y1) + temp_n3_p2(y1);

ris(1) = integral(@(y1) temp_n1(y1), -L, x1, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true) + integral(@(y1) temp_n1(y1), x1, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);
ris(2) = integral(@(y1) temp_n2(y1), -L, x1, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true) + integral(@(y1) temp_n2(y1), x1, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);
ris(3) = integral(@(y1) temp_n3(y1), -L, x1, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true) + integral(@(y1) temp_n3(y1), x1, L, 'AbsTol', 1e-6, 'RelTol', 1e-5, 'ArrayValued', true);

return