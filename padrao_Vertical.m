function A_v = padrao_Vertical(vtThetas, vtAngDT, vtThetas3Db, SLA)
%{
Calculo do padr�o vertical da antena para cada UE

INPUT:
    vtThetas: vetor de angulos de eleva��o entre a BS e os UE's. [graus]
    vtAngDT: vetor com angulos de downtild da antena [graus]  
    vtThetas3Db: vetor com largura de feixe 3db (50% maximo) na vertical [graus] 
    SLA: limite de n�vel do l�bulo lateral [dB]

OUTPUT:
    A_v: padr�o de radia��o vertical da antena [dB]

%}
    
    % calculo do padr�o de radia��o vertical 
    A_v = -min(12.*(((vtThetas-vtAngDT)./vtThetas3Db).^2), SLA);  % [dB]
    
end