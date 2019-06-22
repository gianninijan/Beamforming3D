function [ mtAngsStering ] = geradorAngsDIRECAO( Bh, M, FatorSetor )
    % calcula o grupo de Angs. de DIREÇÃO p/ cada setor
    % 
    % IN:
    %    Bh - número de grupos horizontais
    %    M - número de células
    %    FatorSetor - Setorização da célula
    %
    % OUT:
    %     mtAngsStering - Matriz com angs de Direção de cada setor --> [Linha_X, Coluna_Y] = [Setor_X, Ang_Direção_grupo_Y]
    
    % ang. AZIMUTAL 
    angHorInicCelula = 0;     % inicial da célula em, GRAUS (º)    
    angHorFinalCelula = 360;  % final da célula em, GRAUS (º)
    
    % passo entre os angulos de DIREÇÃO do grupo
    passo = (angHorFinalCelula - angHorInicCelula)/(Bh*FatorSetor);
    
    % gerando os angs. de DIREÇÃO de UMA CÉLULA [em, GRAUS (º)]
    angsFiSt_gr = linspace(angHorInicCelula + (passo/2), angHorFinalCelula - (passo/2), Bh*FatorSetor);
    
    % gerando os angs. de DIREÇÃO de cada SETOR [em, GRAUS (º)]
    mtAngsStering = reshape(angsFiSt_gr, Bh, FatorSetor);    
    mtAngsStering = mtAngsStering';                          
    mtAngsStering = repmat(mtAngsStering, M, 1);     % [linha_X, coluna_Y] = [setor_X, Ang_Direção_grupo_Y] 
end
