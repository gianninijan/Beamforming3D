function [ mtAngsStering, mtAngsTild_gr ] = geracaoGrupos(Bh, Bv, FatorSetor, M)
    
    
    angHorInic = 0;   
    angHorFinal = 360;
    passo = (angHorFinal - angHorInic)/(Bh*FatorSetor);
    angsFiSt_gr = linspace(angHorInic + (passo/2), angHorFinal - (passo/2), Bh*FatorSetor);

    % calculando a matriz de angulos de STEERING p/ cada célula
    mtAngsFiSt_gr = reshape(angsFiSt_gr, Bh, FatorSetor);    
    mtAngsFiSt_gr = mtAngsFiSt_gr';                          
    mtAngsStering = repmat(mtAngsFiSt_gr, M, 1);     % [linha, coluna] = [setor, grupos de angs. de STEERING para o setor (linha)]

    % angDownTild_gr = linspace(-90, 90, Bv);
    angDownTild_gr = zeros(1, Bv);
    angTILD_Inicio = 4;
    passo = 2;
    %angTILD_Final = 10;
    %angDownTild_gr = [angTILD_Inicio:2:angTILD_Final];
    for jj = 1:Bv, 
        angDownTild_gr(jj) = angTILD_Inicio + (jj - 1)*passo;
    end

    % matriz de angulos de INCLINAÇÃO p/ uma celula;
    mtAngsTild_gr = repmat(angDownTild_gr, FatorSetor*M, 1); % [linha, coluna] = [setor, grupos do ang. Inclinação (TILD) por setor] 
end

