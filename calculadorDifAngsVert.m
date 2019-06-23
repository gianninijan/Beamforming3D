function [ mtdifAngsVer_gr ] = calculadorDifAngsVert( S, numUE, mtThetaUE, mtUeSector, mtAngsTild_gr)
    % calcula a diferença entre o âng. ELEVAÇÃO de cada UE e o ang. de INCLINAÇÃO da BS de cada SETOR
    %   
    % IN:
    %    S - número de SETORES
    %    numUE - número de USUÁRIOS
    %    mtThetaUE - matriz com os angs. de ELEVAÇÃo dos UE's, onde: [linha_X, coluna_Y] = [setor_X, UE_Y]
    %    mtUeSector - matriz com os indices dos UE's por setores (linhas da mt.)
    %    mtAngsTild_gr - matriz com os angs. de INCLINAÇÃO de cada setor (linha da mt.)
    %
    % OUT:
    %    mtdifAngsVer_gr - mt. com as diferenças dos ângulos do plano VERTICAL,
    %    onde [linha_X, coluna_Y] = [Setor_X, UE_Y] = (| \theta_UE_Y - \thetaTILD_BS_X |)
    
    % laço percorrendo cada UE
    for jj = 1:numUE,
        
        % [setor, inst]  = [ nº do setor que UE_jj pertence, instante de tempo q UE_jj está ATIVO]  
        [setor, inst] = find(mtUeSector == jj);
    
        % ang. de ELEVAÇÃO do UE_jj em relação a seu BS_setor, em
        angElevUE = mtThetaUE(setor, jj); % graus (º)
        
        % calculando a dif. entre o ang. ELEVAÇÃO do UE_jj e ang. de INCLINAÇÃO da BS_setor
        mtdifAngsVer_gr(setor, jj) = min(abs(angElevUE - mtAngsTild_gr(setor,:)));
        
        % laço percorrendo os SETORES
        for s = 1:S,
        
            % se o setor_s é diferente do 'setor' do UE_jj (laço externo), então
            if s ~= setor, 
                
                % ang. ELEVAÇÃO do UE_jj em relação a BS_s, em:
                angElevUeRel = mtThetaUE(s,jj);
                
                % ang. de ELEVAÇÂO do UE_ATIVO no INSTANTE_DE_TEMPO 'inst' em relação BS_s 
                angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));
                
                % [ ... , indice] = [... , indice do ang. de INCLINAÇÃO da BS_s]
                [diferenca, indice] = min(abs(angElevUeAtivo - mtAngsTild_gr(s,:)));
                
                % calculando a dif. entre o ang. ELEVAÇÃO do UE_jj e ang. de INCLINAÇÂO da BS_s, em
                mtdifAngsVer_gr(s, jj) = min(abs(angElevUeRel - mtAngsTild_gr(s, indice)));
            end
        end
    end
end

