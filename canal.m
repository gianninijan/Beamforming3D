function h = canal(G, A_t, PL, vtShad, vtFastShad)
%{
Calculo dos coeficientes do canal

INPUT:
    G: ganho da antena BS [linear]
    A_t: padr�o de radia��o da antena total [linear]
    PL: Perda de caminho [linear]
    vtShad: vetor de desvanecimento log-normal
    vtFastShad: vetor de desvanecimento r�pido  

OUTPUT:
    h: vetor de coeficientes do canal entre a BS do setor e os UE's
%}
    
    % calculo dos coeficientes do canal
    h = sqrt(G.*A_t.*PL.*vtShad).*vtFastShad;     

end