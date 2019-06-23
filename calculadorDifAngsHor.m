function [ mtDifAngsHor_gr ] = calculadorDifAngsHor( S, numUE, vtUePos, vtBsSetor, mtUeSector, mtAngsStering )
    % 
    % calcula a diferen�a entre o �ng. AZIMUTAL de cada UE e o ang. de DIRE��O da BS de cada SETOR
    % 
    % IN:
    %   S - n�mero de SETORES
    %   numUE - n�mero de USU�RIOS
    %   vtUePos - vetor com a s posi��es dos UE's na horizontal
    %   mtUeSector - matriz com os indices dos UE's por setores (linhas da mt.)
    %   mtAngsStering - matriz com os angs. de DIRE��O de cada setor (linha da mt.)
    %
    % OUT:
    %     mtDifAngsHor_gr - mt. com as diferen�as dos �ngulos, onde:
    %     [linha_X, coluna_Y] = [Setor_X, UE_Y]
    
    % la�o percorrendo cada UE
    for jj = 1:numUE, 
       
       % [setor, inst]  = [ n� do setor que UE_jj pertence, instante de tempo q UE_jj est� ATIVO]    
       [setor, inst] = find(mtUeSector == jj);
       
       % posi��o do UE_jj em rela��o a sua BS_setor
       pos_x = real(vtUePos(jj)) - real(vtBsSetor(setor));
       pos_y = imag(vtUePos(jj)) - imag(vtBsSetor(setor));
       z_aux = pos_x + 1j*pos_y;
    
       % ang. AZIMUTAL do UE_jj em rela��o a sua BS_setor, em:
       anguloUE180 = (180/pi).*angle(z_aux);  % ~ [-180�, +180�]
       anguloUE360 = wrapTo360(anguloUE180);  % ~ [0�, +360�]
       
       % angs. de DIRE��O p/ BS_setor, em:
       angFiGrupo_360 = mtAngsStering(setor,:);     % ~ [0�, 360] 
       angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % ~ [-180�, +180�]
    
       % calculando a dif. entre o ang. AZIMUTAL do UE_jj e ang. de DIRE��O da BS_setor
       dif1 = min(abs(anguloUE360 - angFiGrupo_360));
       dif2 = min(abs(anguloUE180 - angFiGrupo_180));
       difAngulo = min(dif1, dif2);
       mtDifAngsHor_gr(setor, jj) = difAngulo;
    
       % la�o percorrendo os SETORES 
       for s = 1:S,
            
           % se o setor_s � diferente do 'setor' do UE_jj (la�o externo), ent�o
           if s ~= setor,
               
               % posi��o do UE_jj em rela��o a c�lula de BS_s
               pos_xRel = real(vtUePos(jj)) - real(vtBsSetor(s));
               pos_yRel = imag(vtUePos(jj)) - imag(vtBsSetor(s));
               zRel = pos_xRel + 1j*pos_yRel;
               
               % ang. AZIMUTAL do UE_jj em rela��o a BS_s, em:
               anguloRelUE180 = (180/pi)*angle(zRel);       % ~ [-180�, +180�]
               anguloRelUE360 = wrapTo360(anguloRelUE180);  % ~ [0�, 360�]
               
               % posi��o do UE ATIVO no setor_s no INSTANTE de TEMPO 'inst'
               xBS_rel = real(vtUePos(mtUeSector(s, inst))) - real(vtBsSetor(s));
               yBS_rel = imag(vtUePos(mtUeSector(s, inst))) - imag(vtBsSetor(s));
               zBS_rel = xBS_rel + 1j*yBS_rel;
               
               % ang. AZIMUTAL do zBS_rel em rela��o a BS_s
               angAzimUE_ATIVO180 = (180/pi).*angle(zBS_rel); 
               angAzimUE_ATIVO360 = wrapTo360(angAzimUE_ATIVO180);
               
               % [ ... , indice] = [ ... , indice do ang. de DIRE��O da BS_s] 
               [diferenca, indice] = min(abs(angAzimUE_ATIVO360 - mtAngsStering(s,:)));
            
               % calculando a dif. entre o ang. AZIMUTAL do UE_jj e ang. de DIRE��O da BS_s, em
               ddif1 = min(abs(anguloRelUE360 - mtAngsStering(s,indice)));             % ~ [0�, 360�]
               ddif2 = min(abs(anguloRelUE180 - wrapTo180(mtAngsStering(s,indice))));  % ~ [-180�, +180�]
               mtDifAngsHor_gr(s, jj) = min(ddif1, ddif2);
           end
       end
end

