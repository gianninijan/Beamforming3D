function [ mtAngsInclinacao ] = geradorAngsINCLINACAO(Bv, S, angVertInicCelula, angVertFinalCelula)
    % calcula o grupo de Angs. de INCLINAÇÃO p/ cada setor
    %
    % IN:
    %    Bv - número de grupos verticais
    %    S - número de Setores
    %    FatorSetor - Setorização da célula
    %    angVertInicCelula - ang. de Elevação mínimo em, GRAUS (º)
    %    angVertFinalCelula - ang. de Elevação máximo em, GRAUS (º)
    %
    % OUT:
    %     mtAngsInclinacao - Matriz com angs de INCLINAÇÃO de cada setor --> [Linha_X, Coluna_Y] = [Setor_X, Ang_Direção_grupo_Y]
    
    
    % largura do setor Vertical
    larguraSetorVertical = (angVertFinalCelula - angVertInicCelula)/Bv;
    
    % gerando os angs. de DIREÇÃO de UMA CÉLULA [em, GRAUS (º)]
    angDownTild_gr = zeros(1, Bv);
    
    % CÓDIGO ANTERIOR
    % angTILD_Inicio = 4;
    %  passo = 2;
    %  for jj = 1:Bv, 
    %      angDownTild_gr(jj) = angTILD_Inicio + (jj - 1)*passo;
    %  end
    
    
    % CÓGIGO ATUAL SETORIZAÇÃO IGUALMENTE
    angDownTild_gr = [(angVertInicCelula + (larguraSetorVertical/2)): larguraSetorVertical: (angVertFinalCelula - (larguraSetorVertical/2))]; 
    
    % matriz de angulos de INCLINAÇÃO p/ cada SETOR;
    mtAngsInclinacao = repmat(angDownTild_gr, S, 1); % [linha, coluna] = [setor, grupos do ang. Inclinação (TILD) por setor] 
end

