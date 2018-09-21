function[PL] = PathLoss(d, fc)
    % Calculando o PATH LOSS em DB, onde:
    % d � o vetor de dist�ncias dos UE's [m]
    % fc � a frequ�ncia de portadora  [GHz]
    
    % Path Loss
    PL = 36.7.*log10(d) + 22.7 + 26*log10(fc);

end