function [ mtdifAngsVer_gr ] = calculadorDifAngsVert( S, numUE, mtThetaUE, mtUeSector, mtAngsTild_gr)
    % calcula a diferen�a entre o �ng. ELEVA��O de cada UE e o ang. de INCLINA��O da BS de cada SETOR
    %   
    % IN:
    %    S - n�mero de SETORES
    %    numUE - n�mero de USU�RIOS
    %    mtThetaUE - matriz com os angs. de ELEVA��o dos UE's, onde: [linha_X, coluna_Y] = [setor_X, UE_Y]
    %    mtUeSector - matriz com os indices dos UE's por setores (linhas da mt.)
    %    mtAngsTild_gr - matriz com os angs. de INCLINA��O de cada setor (linha da mt.)
    %
    % OUT:
    %    mtdifAngsVer_gr - mt. com as diferen�as dos �ngulos do plano VERTICAL,
    %    onde [linha_X, coluna_Y] = [Setor_X, UE_Y] = (| \theta_UE_Y - \thetaTILD_BS_X |)
    
    % la�o percorrendo cada UE
    for jj = 1:numUE,
        
        % [setor, inst]  = [ n� do setor que UE_jj pertence, instante de tempo q UE_jj est� ATIVO]  
        [setor, inst] = find(mtUeSector == jj);
    
        % ang. de ELEVA��O do UE_jj em rela��o a seu BS_setor, em
        angElevUE = mtThetaUE(setor, jj); % graus (�)
        
        % calculando a dif. entre o ang. ELEVA��O do UE_jj e ang. de INCLINA��O da BS_setor
        mtdifAngsVer_gr(setor, jj) = min(abs(angElevUE - mtAngsTild_gr(setor,:)));
        
        % la�o percorrendo os SETORES
        for s = 1:S,
        
            % se o setor_s � diferente do 'setor' do UE_jj (la�o externo), ent�o
            if s ~= setor, 
                
                % ang. ELEVA��O do UE_jj em rela��o a BS_s, em:
                angElevUeRel = mtThetaUE(s,jj);
                
                % ang. de ELEVA��O do UE_ATIVO no INSTANTE_DE_TEMPO 'inst' em rela��o BS_s 
                angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));
                
                % [ ... , indice] = [... , indice do ang. de INCLINA��O da BS_s]
                [diferenca, indice] = min(abs(angElevUeAtivo - mtAngsTild_gr(s,:)));
                
                % calculando a dif. entre o ang. ELEVA��O do UE_jj e ang. de INCLINA��O da BS_s, em
                mtdifAngsVer_gr(s, jj) = min(abs(angElevUeRel - mtAngsTild_gr(s, indice)));
            end
        end
    end
end

