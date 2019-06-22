function [ mtAngsStering ] = geradorAngsDIRECAO( Bh, M, FatorSetor )
    % calcula o grupo de Angs. de DIRE��O p/ cada setor
    % 
    % IN:
    %    Bh - n�mero de grupos horizontais
    %    M - n�mero de c�lulas
    %    FatorSetor - Setoriza��o da c�lula
    %
    % OUT:
    %     mtAngsStering - Matriz com angs de Dire��o de cada setor --> [Linha_X, Coluna_Y] = [Setor_X, Ang_Dire��o_grupo_Y]
    
    % ang. AZIMUTAL 
    angHorInicCelula = 0;     % inicial da c�lula em, GRAUS (�)    
    angHorFinalCelula = 360;  % final da c�lula em, GRAUS (�)
    
    % passo entre os angulos de DIRE��O do grupo
    passo = (angHorFinalCelula - angHorInicCelula)/(Bh*FatorSetor);
    
    % gerando os angs. de DIRE��O de UMA C�LULA [em, GRAUS (�)]
    angsFiSt_gr = linspace(angHorInicCelula + (passo/2), angHorFinalCelula - (passo/2), Bh*FatorSetor);
    
    % gerando os angs. de DIRE��O de cada SETOR [em, GRAUS (�)]
    mtAngsStering = reshape(angsFiSt_gr, Bh, FatorSetor);    
    mtAngsStering = mtAngsStering';                          
    mtAngsStering = repmat(mtAngsStering, M, 1);     % [linha_X, coluna_Y] = [setor_X, Ang_Dire��o_grupo_Y] 
end
