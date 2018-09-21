function[PL] = PathLoss(d, fc)
    % Calculando o PATH LOSS em DB, onde:
    % d é o vetor de distâncias dos UE's [m]
    % fc é a frequência de portadora  [GHz]
    
    % Path Loss
    PL = 36.7.*log10(d) + 22.7 + 26*log10(fc);

end