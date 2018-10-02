function A_v = padrao_Vertical(vtThetas, vtAngDT, vtThetas3Db, SLA)
%{
Calculo do padrão vertical da antena para cada UE

INPUT:
    vtThetas: vetor de angulos de elevação entre a BS e os UE's. [graus]
    vtAngDT: vetor com angulos de downtild da antena [graus]  
    vtThetas3Db: vetor com largura de feixe 3db (50% maximo) na vertical [graus] 
    SLA: limite de nível do lóbulo lateral [dB]

OUTPUT:
    A_v: padrão de radiação vertical da antena [dB]

%}
    
    % calculo do padrão de radiação vertical 
    A_v = -min(12.*(((vtThetas-vtAngDT)./vtThetas3Db).^2), SLA);  % [dB]
    
end