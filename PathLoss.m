function[PL] = PathLoss(d, fc)
%{
Calculando o PATH LOSS em DB

INPUT:
    d: vetor de dist�ncias entre as UE's e a BS [m]
    fc: frequ�ncia de portadora [GHz]

OUTPUT:
    PL: vetor de Path Loss
%}    
    
    % Path Loss
    PL = 36.7.*log10(d) + 22.7 + 26*log10(fc);

end