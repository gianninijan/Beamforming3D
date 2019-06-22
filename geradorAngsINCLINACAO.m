function [ mtAngsInclinacao ] = geradorAngsINCLINACAO(Bv, S, angVertInicCelula, angVertFinalCelula)
    % calcula o grupo de Angs. de INCLINA��O p/ cada setor
    %
    % IN:
    %    Bv - n�mero de grupos verticais
    %    S - n�mero de Setores
    %    FatorSetor - Setoriza��o da c�lula
    %    angVertInicCelula - ang. de Eleva��o m�nimo em, GRAUS (�)
    %    angVertFinalCelula - ang. de Eleva��o m�ximo em, GRAUS (�)
    %
    % OUT:
    %     mtAngsInclinacao - Matriz com angs de INCLINA��O de cada setor --> [Linha_X, Coluna_Y] = [Setor_X, Ang_Dire��o_grupo_Y]
    
    
    % largura do setor Vertical
    larguraSetorVertical = (angVertFinalCelula - angVertInicCelula)/Bv;
    
    % gerando os angs. de DIRE��O de UMA C�LULA [em, GRAUS (�)]
    angDownTild_gr = zeros(1, Bv);
    
    % C�DIGO ANTERIOR
    % angTILD_Inicio = 4;
    %  passo = 2;
    %  for jj = 1:Bv, 
    %      angDownTild_gr(jj) = angTILD_Inicio + (jj - 1)*passo;
    %  end
    
    
    % C�GIGO ATUAL SETORIZA��O IGUALMENTE
    angDownTild_gr = [(angVertInicCelula + (larguraSetorVertical/2)): larguraSetorVertical: (angVertFinalCelula - (larguraSetorVertical/2))]; 
    
    % matriz de angulos de INCLINA��O p/ cada SETOR;
    mtAngsInclinacao = repmat(angDownTild_gr, S, 1); % [linha, coluna] = [setor, grupos do ang. Inclina��o (TILD) por setor] 
end

