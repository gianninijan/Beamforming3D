% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture
clear all;
clc;
close all;

%% SETUP SIMULATION  

M = 7;                                              % numero de celulas (1 anel)
%M = 19;                                              % numero de celulas(2 anel)
FatorSetor = 3;                                      % Fator de setorização, i.e, setores/celulas
S = M*FatorSetor;                                    % número de setores. S = {1, 2, 3, ..., }      
UEcadaSetor = 100;                                   % numero de UE's por (micro)setor
numUE = UEcadaSetor*S;                               % numero de UE's total
R = 250;                                             % raio da pequena celula
xBS = 0;                                             % Posição do eixo x da BS
yBS = 0;                                             % Posição do eixo y da BS
bandWidth = 10e+6;                                   % largura de banda 
fc = 3.5;                                            % frequencia de portadora [Ghz] na small-cell
Am = 25;                                             % Atenuação Maxima [dB]
SLA = 20;                                            % Limite de nível de lobulo-lateral [dB]
G_BS = 5;                                            % Ganha da antena BS de pequena celula [dBi]
P_BS = 36;                                           % Potencia de TX da BS por setor [dBm]
sigma = 4;                                           % desvio padrão do vetor de sombreamento.[dB]
No = -174;                                           % densidade espectral de potencia [ dBm/Hz ]
FN = 5;                                              % figura de ruido [5 dB]
H_BS = 10;                                           % altura da antena da Estação-Base [metros] 
H_UE = 1.5;                                          % altura da antena da estação-móvel [metros]

%% GERANDO A POSIÇÃO DA BS DE CADA CÉLULA %%

% vetor da POSIÇÃO das BSs
vtBS = [];

% CÉLULA CENTRAL
vtBS(1) = 0*exp(1j*0);

% 1º ANEL
vtBS(2) = 2*R*exp(1j*0);
vtBS(3) = 2*R*exp(1j*pi/3);
vtBS(4) = 2*R*exp(1j*2*pi/3);
vtBS(5) = 2*R*exp(-1j*pi);
vtBS(6) = 2*R*exp(-1j*2*pi/3);
vtBS(7) = 2*R*exp(-1j*pi/3);

% 2º ANEL
% vtBS(8) = 4*R*exp(1j*0); 
% vtBS(9) = (real(vtBS(3)) + 2*R) + 1j*imag(vtBS(3));
% vtBS(10) = 4*R*exp(1j*pi/3);
% vtBS(11) = 2*sqrt(3)*R*exp(1j*pi/2);
% vtBS(12) = 4*R*exp(1j*2*pi/3);
% vtBS(13) = (real(vtBS(4)) - 2*R) + 1j*imag(vtBS(4));
% vtBS(14) = 4*R*exp(-1j*pi);
% vtBS(15) = real(vtBS(13)) - 1j*imag(vtBS(13));
% vtBS(16) = real(vtBS(12)) - 1j*imag(vtBS(12));
% vtBS(17) = 2*sqrt(3)*R*exp(-1j*pi/2);
% vtBS(18) = real(vtBS(10)) - 1j*imag((vtBS(10)));
% vtBS(19) = real(vtBS(9)) - 1j*imag((vtBS(9)));

% sobrepor gráficos
hold on

% laço percorrendo cada BS p/ plotar suas posições e sua ÁREA DE COBERTURA
for ii = 1:length(vtBS)
    circle(real(vtBS(ii)),imag(vtBS(ii)),R)                     % imprimindo a ÁREA DE COBERTURA da BS_ii
    plot(real(vtBS(ii)),imag(vtBS(ii)),'*r','MarkerSize',16)    % imprimindo a POSIÇÃO da BS_ii
end

% vetor com as POSIÇÕES da BS de CADA SETOR
vtBsSetor = repelem(vtBS, FatorSetor);          % = [BS_s1, BS_s2, BS_s3, BS_s4, BS_s5, ..., BS_sS], onde S é o número de setores do sistema


%% GERANDO A POSIÇÃO DE CADA UE %%

% vetor com as POSIÇÕES dos UE's em RELAÇÃO a ORIGEM (0, 0)
vtUePos = []; 

% vetor com o ANG. AZIMUTAL de INICIO (º) de CADA SETOR 
vtAngIncSetor = repmat([0, 120, 240], 1, M);

% laço percorre o numero de UE com intuito de gerar a posição de cada UE p/ q tenhamos um UE ATIVO em cada setor por SLOT de TEMPO
for ii = 1:numUE,
    
    % LAÇO INFINITO
    while true,
        
        % gerando a posição do UE_ii p/ 
        vtUePos(ii) = (6*R*rand - 3*R) + 1j*(6*R*rand - 3*R);      % 1 ANEL
        %vtUePos(ii) = (10*R*rand - 5*R) + 1j*(10*R*rand - 5*R);   % 2 ANÉIS 
        
        % calculando o setor do UEii
         ind = mod(ii, S);      % ind = indice do setor na qual UE_ii faz parte 
         if ind == 0, 
             ind = S;
         end
         
        % distancia do UE_ii p/ BS_ind
        distUE = norm(vtUePos(ii) - vtBsSetor(ind));
        
        % posição do UE_ii em relação a BS_ind
        pos_x = real(vtUePos(ii)) - real(vtBsSetor(ind));
        pos_y = imag(vtUePos(ii)) - imag(vtBsSetor(ind));
        z_aux = pos_x + 1j*pos_y; 
        
        % ang. AZIMUTAL do UE_ii em relação a BS_ind
        anguloUE = wrapTo360((180/pi)*angle(z_aux));    % ~ [0º, 360º]
        
        % se as CONDIÇÕES abaixo em RELAÇÃO a UE_ii e BS_ind forem VERDADEIRAS, então:
        if (distUE < R) && (distUE >= 10) && (anguloUE >= vtAngIncSetor(ind)) && (anguloUE < (vtAngIncSetor(ind) + (360/FatorSetor))),
            break;      % posição VALIDA p/ UE_ii  e QUEBRA o LAÇO INFINITO
        end
    end
end

% plotando as POSIÇÕES dos UE's    
plot(vtUePos,'*b','MarkerSize',12);       % UE = +AZUL
grid on;                                  
legend('Celula', 'BS', 'UEs')             
title('Cenário de Multi-células densas')
% hold off;


%% CONSTRUINDO A MATRIZ INDICE DO UE'S POR SETORES (linhas) %%

% cada elemento representa o indice do UE no vtUePos q pertence ao Setor da Linha e está ativo no instante de tempo 1.  
mtUeSector = zeros(S, UEcadaSetor);  % [linhas, colunas] = [setores, SLOT de TEMPO]

% matriz que encontra os indices na borda da celula
mtBordaCelula = zeros(S, UEcadaSetor); % [linhas, colunas] = [setores, indices dos UE's]

% vetor com os indices dos UE's do vt 'vtUePos' que estão nas bordas 
vtIndBorda = [];

% laço percorrendo cada UE 
for ii = 1:numUE,
    
    % calcula a celula do usuário ii através da distância dele para cada BS
    [dif, indCel] = min(abs(vtUePos(ii) - vtBS)); % = [ distancia do UE_ii p/ BS_indCell,  numero da CELULA que UE_ii pertence ]
    
    % posição do UE_ii em relação BS_indCel
    pos_x = real(vtUePos(ii)) - real(vtBS(indCel));
    pos_y = imag(vtUePos(ii)) - imag(vtBS(indCel));
    z_aux = pos_x + 1j*pos_y;
    
    % ang. AZIMUTAL de UE_ii em relação a célula_indCel
    angUE = wrapTo2Pi(angle(z_aux));                            % ~ [0, 2*pi] radianos
    
    % calculando qual setor do UEii em relação a sua BS
    indSetor = 0;           % indice do SETOR na qual UE_ii pertence
    
    if (angUE >= 0) && (angUE < (2*pi/3)),
        indSetor = 1;
        
    elseif (angUE >= (2*pi/3)) && (angUE < (4*pi/3)),
         indSetor = 2;
         
    else
        indSetor = 3;
    end
    
    linha = FatorSetor*(indCel-1) + indSetor;
    coluna = find( mtUeSector(linha, :) == 0, 1, 'first');
    mtUeSector(linha, coluna) = ii;
    
    if (abs(z_aux) >= 0.9*R) && (abs(z_aux) <= R),
        colBorda = find( mtBordaCelula(linha, :) == 0, 1, 'first');
        mtBordaCelula(linha, colBorda) = ii; 
    end
    
    % se UE_ii satisfazer as condições abaixo, então o mesmo se encontra na BORDA
    if (abs(z_aux) >= 0.9*R) && (abs(z_aux) <= R),
        indice = length(vtIndBorda);
        vtIndBorda(indice + 1) = ii;
    end
    
end

% TESTANDO A MATRIZ 'mtBordaCelula'
% find(mtBordaCelula(i,:));                                         % retorna os elementos não nulos da linha_i
% plot(vtUePos(mtBordaCelula(ii ,find(mtBordaCelula(ii,:)))),'ro'); % imprimir os valores que estão na borda

% TESTANDO A MATRIZ 'mtUeSector' 
% resultado = vtUePos(mtUeSector(jj,:));  % captura as distancia dos UE's que estão dentro do setor_jj
% plot(resultado,'ro');                   % plota o valores do 'resultados' com um circulo azul  
 
hold off


%% CALCULAR A PERDA DE CAMINHO 

% MATRIZ de DISTÂNCIA de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                            

% laço percorrendo cada SETOR
for ii = 1:S
    
    % laço percorrendo cada UE
    for jj = 1:numUE
        
        % calculando a DISTÂNCIA de UE_jj para BS_ii
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));    % [Linha_x, Coluna_y] = [Setor_x , UE_y]   
    
    end    
end

% matriz de PERDA DE PERCURSO de cada UE para cada BS 
mtPL = PathLoss(mtDist, fc);                                 % [Linha_x, Coluna_y] = [Setor_x, UE_y]

PL = [];

% laço percorrendo cada UE's p/ calcular o PL do UE_ii em relação a sua BS_setor 
for ii = 1:numUE,
  
    [setor, inst] = find(mtUeSector == ii); % = [ nº do setor que UE_ii pertence ,  INTERVALO de TEMPO que o UE_ii está ATIVO]
    
    % Perda de caminho p UE_ii
    PL(ii) = mtPL(setor, ii);

end

figure;                                      
plot(sort(min(mtDist)), sort(PL))
grid on;
xlabel('d (M)')
ylabel('PL (DB)')
title('Perda de Percurso ')


%% DADOS COMUNS PARA OS BEAMFORMINGS %%

% matriz de SHADOWING normal em dB
mtNormal = sigma.*randn(S,numUE);               % [LINHA_X, COLUNA_Y] = [SETOR_X, UE_Y]

% matriz com os ANGs. de ELEVAÇÃO [em, GRAUS (º)]  de cada UE p/ BSs de cada setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));     % [LINHA_X, COLUNA_Y] = [SETOR_X, UE_Y]

% POTÊNCIA do RUÍDO em
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);       % escala LINEAR 

% POTÊNCIA da ANTENA TRANSMISSORA em
Pot = dbm2lin(P_BS);                            % escala LINEAR      

% vetor de DESVANECIMENTO RÁPIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];


%%  FORMATAÇÃO BIDIMENSIONAL (2D) DE FEIXES   %%

% angs de DIREÇÃO FIXO p/ BS de cada setor (i.e, angulo AZIMUTAL de ganho máximo da antena) 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de direção p/ cada CELULA [em, Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo direção p/ cada SETOR [em, Radianos]

% MATRIZ de diferença entre o ang. AZIMUAL de cada UE p/ o ang. de DIREÇÃO da BS de cada setor (em, GRAUS º)
mtdifAngsHor_2D = zeros(S, numUE);      % [linha_X, coluna_Y] = [setor_X, UE_Y]

% laço percorrendo cada UE p/ calcular os valores da matriz 'mtdifAngsHor_2D'
for ii = 1:numUE,
    
    [setor, inst] = find(mtUeSector == ii); % = [ nº do setor que UE_ii pertence ,  INTERVALO de TEMPO que o UE_ii está ATIVO]
    
    % posição do UEii em relação a BS_setor
    pos_x = real(vtUePos(ii)) - real(vtBsSetor(setor));
    pos_y = imag(vtUePos(ii)) - imag(vtBsSetor(setor));
    z_aux = pos_x + 1j*pos_y;  
    
    % ang. AZIMUTAL do UEii em relação a posição da BS_setor, em
    anguloUE180 = (180/pi)*angle(z_aux);         % ~ [-180º, +180º]
    anguloUE360 = wrapTo360(anguloUE180);        % ~ [0º, 360º]
    
    % ang. de DIREÇÃO da BS_setor
    angStrTo360 = (180/pi)*vtAngST(setor);       % ~ [0º, 360º]
    angStrTo180 = wrapTo180(angStrTo360);        % ~ [-180º, +180º]
    
    % calculando a diferença entre o ang. de AZIMUTAL da UEii com o ang. de DIREÇÃO da BS_setor
    mtdifAngsHor_2D(setor, ii) = min(abs(anguloUE360 - angStrTo360), abs(anguloUE180 - angStrTo180));
    
    % laço percorrendo cada SETOR
    for s = 1:S,
        
        % se o SETOR_s é diferente do SETOR_setor na qual UEii pertence (laço externo), então
        if s ~= setor,
            
            % Posição do UEii em relação a BS,s 
            pos_xRel = real(vtUePos(ii)) - real(vtBsSetor(s));
            pos_yRel = imag(vtUePos(ii)) - imag(vtBsSetor(s));
            zRel = pos_xRel + 1j*pos_yRel;
            
            % angulo AZIMUTAL de UE_ii em relação a posição da BS_s
            anguloRelUE180 = (180/pi)*angle(zRel);             % ~ [-180º, +180º]
            anguloRelUE360 = wrapTo360(anguloRelUE180);        % ~ [0º, 360º]
    
            % angulo de DIREÇÃO da BS do setor_s
            angStr360 = (180/pi)*vtAngST(s);                   % ~ [0º, 360º]
            angStr180 = wrapTo180(angStr360);                  % ~ [-180º, +180º]
            
            % calculando a diferença entre o ang. de AZIMUTAL da UEii com o angulo de STEERING da BS_s
            mtdifAngsHor_2D(s, ii) = min(abs(anguloRelUE360 - angStr360), abs(anguloRelUE180 - angStr180));
        end
    end
end

% MATRIZ de diferença entre o angulo ELEVAÇÃO de cada UEs para o angulo de INCLINAÇÃO (TILD) da BS de cada setor (em, GRAUS º)
mtdifAngsVer_2D = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% ang. INCLINAÇÃO FIXO p/ BS de cada setor, em
angDownTild_2d = 8;                % [GRAUS, º] -> Pág. 4836, Simulation Setup

% calculando os valores da matriz 'mtdifAngsVer_2D'
mtdifAngsVer_2D = abs(mtThetaUE - angDownTild_2d);

% LARGURAS DE FEIXE em 3 dB, na
fi3dB_2D = 70;                     % HORIZONTAL [GRAUS] -> Pág. 4837, Simulation Setup
theta3dB_2D = 10;                  % VERTICAL   [GRAUS] -> Pág. 4837, Simulation Setup

% caluclando o padrão de RADIAÇÃO 
Ah_2D = -min(12.*((mtdifAngsHor_2D./fi3dB_2D).^2), Am);      % HORIZONTAL
Av_2D = -min(12.*((mtdifAngsVer_2D./theta3dB_2D).^2), SLA);  % VERTICAL   
A_2D = -min(-(Ah_2D + Av_2D), Am);                           % TOTAL

% COEFICIENTES DO CANAL AO QUADRADO (SEM FAST-FADING), em
H_2d = G_BS + A_2D - mtPL + mtNormal;   % [dB]
h_2d = db2lin(H_2d);                    % ESCALA LINEAR

% vetor de SINR p/ 2DBF
Y2D = [];   

% laço percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
    
    [setor, inst] = find(mtUeSector == ii);  % = [ nº do setor que UE_ii pertence ,  INTERVALO de TEMPO que o UE_ii está ATIVO]
    
    % calcula os coeficientes de canais INTERFERENTES p/ UE_ii
    aux = sum(h_2d(:, ii)) - h_2d(setor, ii);
    
    % calculo da SNIR do UE_ii
    Y2D(ii) = (Pot*h_2d(setor, ii))/(Pot*aux + PN);
    
end

% SINR em dB p/ 2DBF
Y2D_dB = 10*log10(Y2D);

% EFICIÊNCIA ESPECTRAL 
R_2d = log2(1 + Y2D);        % de cada UE
RTotal_2D = sum(R_2d, 2);    % TOTAL 
Rmedia_2D = RTotal_2D/numUE; % MÉDIA

% CAPACIDADE
C_2d = bandWidth.*log2(1 + Y2D);                            % de cada UE 
CTotal_2D = sum(C_2d, 2);                                   % TOTAL
Cmedia_2D = CTotal_2D/numUE;                                % MÉDIA
Cborda_2d = bandWidth.*log2(1 + Y2D(vtIndBorda));           % de cada UE na BORDA
CMediaborda_2d = (sum(Cborda_2d, 2))/(length(vtIndBorda));  % média dos UE's na Borda

figure;             % gera uma nova figura 
cdfplot(Y2D_dB)     % plota a CDF do SINR p/ 2DBF em dB 
hold on;


%% FORMATAÇÃO DE FEIXES TRIDIMENSIONAIS (3D) POR USUÁRIOS (3DBF UE-SPECIFIC) 

% MATRIZ de DIFERENÇA entre o ang. AZIMUTAL de cada UEs para o ang. de DIREÇÃO da BS de cada setor (em, GRAUS º)
mtDifAngHor_Esp = zeros(S, numUE);  % [linha_X, coluna_Y] = [setor_X, UE_Y]

% laço percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngHor_Esp'
for jj = 1:numUE,
    
    [setor, inst] = find(mtUeSector == jj); % = [ nº do setor que UE_jj pertence, INTERVALO DE TEMPO que o UE_jj está ativo]
     
    % laço percorrendo cada SETOR
    for s = 1:S,
        
        % se o SETOR_s é diferente do SETOR_setor na qual UEjj pertence (laço externo), então
        if s ~= setor,
            
            % Posição do UEjj em relação a BS_s 
            pos_xRel = real(vtUePos(jj)) - real(vtBsSetor(s));
            pos_yRel = imag(vtUePos(jj)) - imag(vtBsSetor(s));
            zRel = pos_xRel + 1j*pos_yRel;
            
            % angulo AZIMUTAL de UE_jj em relação a posição da BS_s, em
            anguloRelUE180 = (180/pi)*angle(zRel);             % ~ [-180º, +180º]
            anguloRelUE360 = wrapTo360(anguloRelUE180);        % ~ [0º, 360º]
            
            % Posição do UE_ATIVO no setor_s no instante 'inst'
            xBS_rel = real(vtUePos(mtUeSector(s, inst))) - real(vtBsSetor(s));
            yBS_rel = imag(vtUePos(mtUeSector(s, inst))) - imag(vtBsSetor(s));
            zBS_rel = xBS_rel + 1j*yBS_rel;
            
            % angulo de STEERING da BS (= ang. AZIMUTAL do UE ATIVO) do setor_s no SLOT DE TEMPO 'inst'
            angsStrBS180 = (180/pi).*angle(zBS_rel);       % [-180º, +180º]
            angsStrBS360 = wrapTo360(angsStrBS180);        % [0º, 360º]
            
            % calculando a diferença entre o ang. de AZIMUTAL da UEjj com o ang. de STEERING da BS_s  
            mtDifAngHor_Esp(s, jj) = min(abs(anguloRelUE360 - angsStrBS360), abs(anguloRelUE180 - angsStrBS180));    
        end 
    end 
end


% matriz de DIFERENÇA entre o ang. ELEVAÇÃO de cada UE p/ o ang. de INCLINAÇÃO da BS de cada setor (em, GRAUS º)
mtDifAngsVer_Esp = zeros(S, numUE);  % [linha_X, coluna_Y] = [setor_X, UE_Y]

% laço percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngsVer_Esp'
for jj = 1:numUE,
    
    [setor, inst] = find(mtUeSector == jj);  % [ nº do setor que UEjj pertence , INTERVALO de TEMPO que o UEjj está ativo ]
        
    % laço percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' é diferente do 'setor' UEjj (laço externo), então
        if s ~= setor,
            
            % angulo de ELEVAÇÂO do UEjj em relação BS_s
            angElevRelUE_Esp = mtThetaUE(s, jj);                   % ~ [º,º]
            
            % ang. de INCLINAÇÃO (TILD) da BS_s será igual o ang. de Elevação do UE ativo no SLOT DE TEMPO 'inst' nesse setor_s
            angVertUeAtivo = mtThetaUE(s, mtUeSector(s, inst));  % = ang. de INCLINAÇÃO da BS_s
            
            % calculando a diferença entre o ang. de ELEVAÇÃO do UEjj com BS_s E o ang. de INCLINAÇÂO (TILD) da BS_s 
            mtDifAngsVer_Esp(s, jj) = abs(angElevRelUE_Esp - angVertUeAtivo); 
        end
    end
end

% TENSORES para calcular o padrão RADIAÇÃO
Av_Esp = []; % VERTICAL  p/ cada ang. theta_3dB
Ah_Esp = []; % HORIZONTAL p/ cada ang. \phi_3dB
A_ESP = [];  % TOTAL  

% vetor de LARGURAS DE FEIXE em 3 dB, na
% vtFi3dB_Esp = [20 10 5];                  % HORIZONTAL [GRAUS] --> (Fig. 3, p. 4836)
% vtTheta3dB_Esp = [10 10 5];               % VERTICAL   [GRAUS] --> (Fig. 3, p. 4836)

% LARGURAS DE FEIXE em 3 dB, na
fi3dB_3DBFbyUE = 10;                 % HORIZONTAL [GRAUS]
theta3dB_3DBFbyUE = 10;              % VERTICAL   [GRAUS]

% calculando o PADRÃO de RADIAÇÃO
Ah_Esp = -min(12.*((mtDifAngHor_Esp./fi3dB_3DBFbyUE).^2), Am);       % VERTICAL
Av_Esp = -min(12.*((mtDifAngsVer_Esp./theta3dB_3DBFbyUE).^2), SLA);  % HORIZONTAL
A_Esp = -min(-(Ah_Esp + Av_Esp), Am);                                % TOTAL 

% COEFICIENTE do CANAL ao QUADRADO (sem fast-fading) em,
H_ESP = G_BS + A_Esp - mtPL + mtNormal;   % dB
h_esp = db2lin(H_ESP);    % ESCALA LINEAR, onde: [linha_X, coluna_Y] = [Setor_X, UE_Y]

% vetr SINR para 3DBF POR UE
Y_ESP = [];  

% laço percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
    
    % 'setor' é o setor que UE_ii pertence e o instante 'inst' que o UEjj está ativo
    [setor, inst] = find(mtUeSector == ii);      
        
    % soma o h_esp interferentes do UE_ii
    aux = sum(h_esp(:, ii)) - h_esp(setor, ii);
        
    % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
    Y_ESP(ii) = (Pot*h_esp(setor, ii))/(Pot*aux + PN);
    
end

% SINR p/ 3DBF por UE, em dB
YESP_dB = 10.*log10(Y_ESP);

% EFICIÊNCIA ESPECTRAL 
R_esp = log2(1 + Y_ESP);            % DE CADA UE
RTotal_esp = sum(R_esp, 2);         % TOTAL
Rmedia_esp = RTotal_esp./numUE;     % MÉDIA

% CAPACIDADE DE CADA USUÁRIO
C_esp = bandWidth.*log2(1 + Y_ESP);                                   % DE CADA UE 
CTotal_esp = sum(C_esp, 2);                                           % TOTAL 
CMedia_esp = CTotal_esp./numUE;                                       % MÉDIA
CBorda_Esp = bandWidth.*log2(1 + Y_ESP(vtIndBorda));                  % NA BORDA
CMediaBorda_Esp = (sum(CBorda_Esp, 2))./(length(vtIndBorda));         % MÉDIA DOS USUÁRIOS NA BORDA

% PLOTANDO A CDF DA SINR_dB
hold on
cdfplot(YESP_dB)


%% FORMATAÇÃO DE FEIXES TRIDIMENSIONAIS (3D) POR GRUPOS DE UE (3DBF UE-GROUP) 

% numero de feixes
Bh = 8;                           % HORIZONTAIS p/ cada setor
Bv = 2;                           % VERTICAIS p/ cada setor
B = Bh*Bv;                        % TOTAIS (ou, GRUPOS)

% vetor de PARTIÇÕES
vtBh = [8, 16, 16, 32, 32, 64, 128];  %16, 16, 32, 64, 64];    % HORIZONTAIS
vtBv = [2, 2, 4, 4, 8, 8, 8];     %4, 4, 4, 8];       % VERTICAIS
vtB = vtBh.*vtBv;       % TOTAL (ou, GRUPOS)

% Tensor para os coeficientes de canal para cada par de Grupos
h_gr = [];

% laço percorrendo cada grupo
for jj = 1:length(vtB), 
    
    % matriz com o grupo de ang. de DIREÇÃO p/ cada SETOR
    mtAngsStering = geradorAngsDIRECAO(vtBh(jj), M, FatorSetor);  % [linha_X, coluna_Y] = [setor_X, ang_Direcao_Y do setor_X]
    
    % matriz DIFERENÇA do ang. AZIMUTAL de cada UE para o ang. DIREÇÃO da BS de cada SETOR [em, GRAUS (º)]
    mtDifAngsHor_gr = calculadorDifAngsHor( S, numUE, vtUePos, vtBsSetor, mtUeSector, mtAngsStering ); % [linha_X, coluna_Y] = [setor_X, UE_Y]

    % matriz de angulos de INCLINAÇÃO p/ cada SETOR;
    mtAngsTild_gr = geradorAngsINCLINACAO(vtBv(jj), S, 0, max(max(mtThetaUE)));       % [linha_X, coluna_Y] = [setor_X, ang_Inclinacao_Y do setor_X]

    % matriz DIFERENÇA do ang de ELEVAÇÃO de cada UE p/ o ang. de INCLINAÇÃO da BS de cada SETOR (em, GRAUS º)
    mtdifAngsVer_gr = calculadorDifAngsVert( S, numUE, mtThetaUE, mtUeSector, mtAngsTild_gr);
    
    % LARGURA DE FEIXES DE 3 dB, na
    fi3dB_3DBFbyGr = 10;     % HORIZONTAL [GRAUS]
    Theta3dB_3DBFbtGr = 10;  % VERTICAL [GRAUS]
    
    % calculando o PADRÃO de RADIAÇÃO
    Ah_gr = -min(12.*((mtDifAngsHor_gr./fi3dB_3DBFbyGr).^2), Am);       % HORIZONTAL
    Av_gr = -min(12.*((mtdifAngsVer_gr./Theta3dB_3DBFbtGr).^2), SLA);   % VERTICAL
    A_GR = -min(-(Av_gr + Ah_gr), Am);                                  % TOTAL
    
    % COEFICIENTES do CANAL ao QUADRADO (SEM FAST-FADING)
    H_GR = G_BS + A_GR - mtPL + mtNormal;   % [dB]
    h_gr(:,:,jj) = db2lin(H_GR);                    % ESCALA LINEAR
end

% vetor SINR para Beamforming GRUPO
Y_GR = [];   % SINR

% laço percorrendo cada par (Bh, Bv)
for grupo = 1:length(vtB), 
    
    % laço percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,

        [setor, inst] = find(mtUeSector == ii);  % = [nº do setor q UE_ii pertence, INTERVALO de TEMPO que o UE_ii está ativo]

        % soma h_gr interferentes para cada UE 
        aux = sum(h_gr(:, ii, grupo)) - h_gr(setor, ii, grupo);

        % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
        Y_GR(grupo, ii) = (Pot*h_gr(setor, ii, grupo))/(Pot*aux + PN);     

    end
end

% SINR em dB
YGR_dB = 10.*log10(Y_GR);
hold on;
cdfplot(YGR_dB(1,:))
cdfplot(YGR_dB(2,:))
cdfplot(YGR_dB(3,:))
cdfplot(YGR_dB(4,:))
cdfplot(YGR_dB(5,:))
cdfplot(YGR_dB(6,:))
cdfplot(YGR_dB(7,:))
legend('FFC', 'FFU', 'FFG - 16 Grupos', 'FFG - 32 Grupos', 'FFG - 64 Grupos', 'FFG - 128 Grupos')
xlabel('SINR (dB)')
ylabel('CDF')
title('')
hold off;

%% GRÁFICOS DA CAPACIDADE

