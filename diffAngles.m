function mtDif = diffAngles(ang1, ang2)
    %{  
        Calculo da diferença minima entre os angulos dos vetores de entrada

        INPUT:
             ang1 - vetor de angulos [radianos]
             ang2 - vetor de angulos [radianos]
        
        OUTPUT:
            mtDif - matriz da diferença minima entre angulos [graus]
    %}
    
    [X, Y] = meshgrid(ang1, ang2);  % cria 2 matrizes com linhas e colunas repetidas de ang1 e ang2 
    
    %res = [X(:) Y(:) wrapToPi(Y(:))]; % cada linha - [elem_X elem_Y elem_Y_180]
    % res = [X(:) Y(:)];                 % cada linha - [elem_X elem_Y]
     
    mtDif = (180/pi).*abs(X(:)- Y(:));
    
    mtDif = reshape(mtDif, length(ang2), length(ang1)); 
    mtDif = mtDif';
    
    % percorrendo a linha da matriz res
%     for ii = 1:size(res,1)
%         mtDif(ii,:) = res(ii,1) - res(ii, 2:3);
%     end
%     
%     mtDif = abs(mtDif);
%     mtDif = (180/pi).*min(mtDif, [], 2);
%     
%     mtDif = reshape(mtDif, length(ang1), length(ang2));
end