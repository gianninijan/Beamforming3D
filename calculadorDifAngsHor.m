function [ mtDifAngsHor_gr ] = calculadorDifAngsHor( S, numUE, vtUePos, vtBsSetor, mtUeSector, mtAngsStering )
    % 
    % calcula a diferença entre o âng. AZIMUTAL de cada UE e o ang. de DIREÇÃO da BS de cada SETOR
    % 
    % IN:
    %   S - número de SETORES
    %   numUE - número de USUÁRIOS
    %   vtUePos - vetor com a s posições dos UE's na horizontal
    %   mtUeSector - matriz com os indices dos UE's por setores (linhas da mt.)
    %   mtAngsStering - matriz com os angs. de DIREÇÃO de cada setor (linha da mt.)
    %
    % OUT:
    %     mtDifAngsHor_gr - mt. com as diferenças dos ângulos, onde:
    %     [linha_X, coluna_Y] = [Setor_X, UE_Y]
    
    % laço percorrendo cada UE
    for jj = 1:numUE, 
       
       % [setor, inst]  = [ nº do setor que UE_jj pertence, instante de tempo q UE_jj está ATIVO]    
       [setor, inst] = find(mtUeSector == jj);
       
       % posição do UE_jj em relação a sua BS_setor
       pos_x = real(vtUePos(jj)) - real(vtBsSetor(setor));
       pos_y = imag(vtUePos(jj)) - imag(vtBsSetor(setor));
       z_aux = pos_x + 1j*pos_y;
    
       % ang. AZIMUTAL do UE_jj em relação a sua BS_setor, em:
       anguloUE180 = (180/pi).*angle(z_aux);  % ~ [-180º, +180º]
       anguloUE360 = wrapTo360(anguloUE180);  % ~ [0º, +360º]
       
       % angs. de DIREÇÃO p/ BS_setor, em:
       angFiGrupo_360 = mtAngsStering(setor,:);     % ~ [0º, 360] 
       angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % ~ [-180º, +180º]
    
       % calculando a dif. entre o ang. AZIMUTAL do UE_jj e ang. de DIREÇÂO da BS_setor
       dif1 = min(abs(anguloUE360 - angFiGrupo_360));
       dif2 = min(abs(anguloUE180 - angFiGrupo_180));
       difAngulo = min(dif1, dif2);
       mtDifAngsHor_gr(setor, jj) = difAngulo;
    
       % laço percorrendo os SETORES 
       for s = 1:S,
            
           % se o setor_s é diferente do 'setor' do UE_jj (laço externo), então
           if s ~= setor,
               
               % posição do UE_jj em relação a célula de BS_s
               pos_xRel = real(vtUePos(jj)) - real(vtBsSetor(s));
               pos_yRel = imag(vtUePos(jj)) - imag(vtBsSetor(s));
               zRel = pos_xRel + 1j*pos_yRel;
               
               % ang. AZIMUTAL do UE_jj em relação a BS_s, em:
               anguloRelUE180 = (180/pi)*angle(zRel);       % ~ [-180º, +180º]
               anguloRelUE360 = wrapTo360(anguloRelUE180);  % ~ [0º, 360º]
               
               % posição do UE ATIVO no setor_s no INSTANTE de TEMPO 'inst'
               xBS_rel = real(vtUePos(mtUeSector(s, inst))) - real(vtBsSetor(s));
               yBS_rel = imag(vtUePos(mtUeSector(s, inst))) - imag(vtBsSetor(s));
               zBS_rel = xBS_rel + 1j*yBS_rel;
               
               % ang. AZIMUTAL do zBS_rel em relação a BS_s
               angAzimUE_ATIVO180 = (180/pi).*angle(zBS_rel); 
               angAzimUE_ATIVO360 = wrapTo360(angAzimUE_ATIVO180);
               
               % [ ... , indice] = [ ... , indice do ang. de DIREÇÃO da BS_s] 
               [diferenca, indice] = min(abs(angAzimUE_ATIVO360 - mtAngsStering(s,:)));
            
               % calculando a dif. entre o ang. AZIMUTAL do UE_jj e ang. de DIREÇÂO da BS_s, em
               ddif1 = min(abs(anguloRelUE360 - mtAngsStering(s,indice)));             % ~ [0º, 360º]
               ddif2 = min(abs(anguloRelUE180 - wrapTo180(mtAngsStering(s,indice))));  % ~ [-180º, +180º]
               mtDifAngsHor_gr(s, jj) = min(ddif1, ddif2);
           end
       end
end

