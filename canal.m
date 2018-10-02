function h = canal(G, A_t, PL, vtShad, vtFastShad)
%{
Calculo dos coeficientes do canal

INPUT:
    G: ganho da antena BS [linear]
    A_t: padrão de radiação da antena total [linear]
    PL: Perda de caminho [linear]
    vtShad: vetor de desvanecimento log-normal
    vtFastShad: vetor de desvanecimento rápido  

OUTPUT:
    h: vetor de coeficientes do canal entre a BS do setor e os UE's
%}
    
    % calculo dos coeficientes do canal
    h = sqrt(G.*A_t.*PL.*vtShad).*vtFastShad;     

end