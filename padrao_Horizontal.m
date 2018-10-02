function A_h = padrao_Horizontal( vtPhi, phiSTd, phi3DBd, Am)
%{
Calculo do padr�o horizontal da antena para cada UE

INPUT:
    vtPhi: vetor de angulos azimutais dos UE's [radianos]
    phiSTd: angulo de boresight da antena, i.e, angulo de ganho maximo da antena na horizontal  [graus] 
    phi3DBd: largura de feixe em 3dB [graus]
    Am: Raz�o frente-tr�s [dB]
    
OUTPUT:
   A_h: padr�o de radia��o horizontal da antena [DB]
%}

    vtPhid = (180/pi)*vtPhi;    % angulos azimutais dos UE's em GRAUS
    
    A_h = -min((12.*(((vtPhid - phiSTd)/phi3DBd).^2)), Am);
    
end