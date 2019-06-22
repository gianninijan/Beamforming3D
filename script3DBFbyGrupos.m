%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture
clear all;
clc;
close all;


%% SETUP SIMULATION  

%M = 7;                                               % numero de celulas (1 anel)
M = 19;                                            % numero de celulas (2 anel)
FatorSetor = 3;                                      % Fator de setoriza��o, i.e, setores/celulas
S = M*FatorSetor;                                    % n�mero de setores. S = {1, 2, 3, ..., }      
UEcadaSetor = 50;                                   % numero de UE's por (micro)setor
numUE = UEcadaSetor*S;                               % numero de UE's total
R = 250;                                             % raio da pequena celula
xBS = 0;                                             % Posi��o do eixo x da BS
yBS = 0;                                             % Posi��o do eixo y da BS
bandWidth = 10e+6;                                   % largura de banda 
fc = 3.5;                                            % frequencia de portadora [Ghz] na small-cell
Am = 25;                                             % Atenua��o Maxima [dB]
SLA = 20;                                            % Limite de n�vel de lobulo-lateral [dB]
G_BS = 5;                                            % Ganha da antena BS de pequena celula [dBi]
P_BS = 36;                                           % Potencia de TX da BS por setor [dBm]
sigma = 4;                                           % desvio padr�o do vetor de sombreamento.[dB]
No = -174;                                           % densidade espectral de potencia [ dBm/Hz ]
FN = 5;                                              % figura de ruido [5 dB]
H_BS = 10;                                           % altura da antena da Esta��o-Base [metros] 
H_UE = 1.5;                                          % altura da antena da esta��o-m�vel [metros]


%% GERANDO A POSI��O DA BS DE CADA C�LULA %%

% VALORES DE TESTE
vtBS = [];

% celula central
vtBS(1) = 0*exp(1j*0);

% primeiro anel
vtBS(2:7) = [2*R*exp(1j*0), 2*R*exp(1j*pi/3), 2*R*exp(1j*2*pi/3), 2*R*exp(-1j*pi), 2*R*exp(-1j*2*pi/3), 2*R*exp(-1j*pi/3)];

% segundo anel
vtBS(8) = 4*R*exp(1j*0); 
vtBS(9) = (real(vtBS(3)) + 2*R) + 1j*imag(vtBS(3));
vtBS(10) = 4*R*exp(1j*pi/3);
vtBS(11) = 2*sqrt(3)*R*exp(1j*pi/2);
vtBS(12) = 4*R*exp(1j*2*pi/3);
vtBS(13) = (real(vtBS(4)) - 2*R) + 1j*imag(vtBS(4));
vtBS(14) = 4*R*exp(-1j*pi);
vtBS(15) = real(vtBS(13)) - 1j*imag(vtBS(13));
vtBS(16) = real(vtBS(12)) - 1j*imag(vtBS(12));
vtBS(17) = 2*sqrt(3)*R*exp(-1j*pi/2);
vtBS(18) = real(vtBS(10)) - 1j*imag((vtBS(10)));
vtBS(19) = real(vtBS(9)) - 1j*imag((vtBS(9)));

% sobrepor gr�ficos
hold on

% la�o percorrendo cada BS p/ plotar suas posi��es e sua �REA DE COBERTURA
for ii = 1:length(vtBS)
    circle(real(vtBS(ii)),imag(vtBS(ii)),R)
    plot(real(vtBS(ii)),imag(vtBS(ii)),'xr')
end

vtBsSetor = repelem(vtBS, FatorSetor);             % vetor posi��o da BS's p/ cada setor 


%% GERANDO A POSI��O DE CADA UE ~ U(10,R) %%

% vetor p/ as posi��es aleat�rios dos UE's em RELA��O a ORIGEM (0, 0) do sistema
vtUePos = []; 

% vetor com o angulo Horizontal de inicio (em, graus �) de cada setor 
vtAngIncSetor = repmat([0, 120, 240], 1, M);

% la�o percorre o numero de UE's p/ gerar a posi��o de cada UE de forma que teremos um UE ATIVO em cada setor por SLOT de TEMPO
for ii = 1:numUE,
    
    while true,
        
        % calculando a posi��o do UEii 
        % vtUePos(ii) = (6*R*rand - 3*R) + 1j*(6*R*rand - 3*R);      % p/ 1 anel
        vtUePos(ii) = (10*R*rand - 5*R) + 1j*(10*R*rand - 5*R);   % p/ 2 an�is 
        
        % calculando o indice na qual representa o setor do UEii
         ind = mod(ii, S);
         if ind == 0, 
             ind = S;
         end
         
        % distancia do UE p/ BS do setor
        distUE = norm(vtUePos(ii) - vtBsSetor(ind));
        
        % posi��o do UEii em rela��o a BS do setor
        pos_x = real(vtUePos(ii)) - real(vtBsSetor(ind));
        pos_y = imag(vtUePos(ii)) - imag(vtBsSetor(ind));
        z_aux = pos_x + 1j*pos_y; 
        
        % angulo do UE em rela��o a posi��o da BS
        anguloUE = wrapTo360((180/pi)*angle(z_aux)); % ~[0�, 360�]
        
        % se distancia do UE p/ BS do setor for menor que o Raio e maior que 10, ent�o:
        if (distUE < R) && (distUE >= 10) && (anguloUE >= vtAngIncSetor(ind)) && (anguloUE < (vtAngIncSetor(ind) + (360/FatorSetor))),
            % plot(vtUePos(ii),'*b','MarkerSize',12);       % plotando os UE's (+azul)
            break;
        end
    end
end

% plotando as POSI��ES dos UE's    
plot(vtUePos,'*b','MarkerSize',12);       % UE = +AZUL
grid on;
legend('Celula', 'BS', 'UEs')
title('Cen�rio de Multi-c�lulas densas')
% hold off;


%% CONSTRUINDO A MATRIZ INDICE DO UE'S POR SETORES (linhas) %%

% cada elemento representa o indice do UE no vtUePos q pertence ao Setor da Linha e est� ativo no instante de tempo 1.  
mtUeSector = zeros(S, UEcadaSetor);  % [linhas, colunas] = [setores, SLOT de TEMPO]

% matriz que encontra os indices na borda da celula
mtBordaCelula = zeros(S, UEcadaSetor); % [linhas, colunas] = [setores, indices dos UE's]
vtIndBorda = [];

% la�o percorrendo cada UE 
for ii = 1:numUE,
    
    % calcula a celula do usu�rio ii atrav�s da dist�ncia dele para cada BS
    % [dif, indCel] = [diferen�a do UE_ii p/ BS_setor, indice da celula que UE_ii pertence]
    [dif, indCel] = min(abs(vtUePos(ii) - vtBS)); 
    
    % Encontrando o angulo do UEii em rela��o a sua BS
    pos_x = real(vtUePos(ii)) - real(vtBS(indCel));
    pos_y = imag(vtUePos(ii)) - imag(vtBS(indCel));
    z_aux = pos_x + 1j*pos_y;
    angUE = wrapTo2Pi(angle(z_aux));                            % angUE ~ [0, 2*pi] radianos
    
    % calculando qual setor do UEii em rela��o a sua BS
    indSetor = 0;
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
    
    if (abs(z_aux) >= 0.9*R) && (abs(z_aux) <= R),
        indice = length(vtIndBorda);
        vtIndBorda(indice + 1) = ii;
    end
end


% find(mtBordaCelula(i,:)); % retorna os elementos n�o nulos da linha_i

% TESTANDO O CODIGO DA MATRIZ 'mtUeSector' 
%resultado = vtUePos(mtUeSector(1,:));   % captura as distancia dos UE's que est�o dentro do setor X
%plot(resultado,'ro');                   % plota o valores do 'resultados' com um circulo azul  
 hold off


%% CALCULAR A PERDA DE CAMINHO (en: Path Loss) DO UE'S P/ BS'S %%

% MATRIZ de DIST�NCIA de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                            

% la�o percorrendo cada setor (coluna)
for ii = 1:S
    
    % la�o percorrendo os UE's (linha)
    for jj = 1:numUE
        
        % calculando a DIST�NCIA de UEjj para BSii
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % [Linha, Coluna] = [BS do Setor, UE]   
    
    end    
end

% matriz de PATH LOSS de cada UE para cada BS 
mtPL = PathLoss(mtDist, fc);                   % [linhas, colunas] = [SETORES, UE's]

PL = [];

% la�o percorrendo os UE's p/ calcular o PL 
for ii = 1:numUE,

    % 'setor' � o n� do setor que UEjj pertence; 'inst' � o SLOT de TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == ii);
    
    % Perda de caminho de cada UE
    PL(ii) = mtPL(setor, ii);
end

% PL = min(mtPL);                                % valor minimo
% vtDistUEtoBS = min(mtDist);

%[vtDistUEtoBS, POS] = min(mtDist, [], 1);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(min(mtDist)), sort(PL))
grid on;
xlabel('d (M)')
ylabel('PL (DB)')
title('Perda de Percurso ')


%% DADOS COMUNS PARA OS BEAMFORMINGS %%

% vetor de DESVANECIMENTO R�PIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];

% matriz de SHADOWING normal em dB
mtNormal = sigma.*randn(S,numUE);    % [Linhas, Colunas] = [Setores, UE's]

% matriz de Angulos de ELEVA��O [em, GRAUS (�)]  de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [LINHA, COLUNA] = [SETOR, UE]

% Pot�ncia do Ru�do Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Pot�ncia da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% angulos de STEERING FIXO p/ BS de cada setor, i.e, angulo azimutal na qual teremos o ganho m�ximo da antena 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de boresight p/ cada celula [Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS] -> P�g. 4837, Simulation Setup
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS] -> P�g. 4837, Simulation Setup
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS] -> P�g. 4836, Simulation Setup

% MATRIZ de diferen�a entre o angulo AZIMUAL de cada UE's para o angulo de BORESIGHT da BS de cada setor (em, GRAUS �)
mtdifAngsHor_2D = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtdifAngsHor_2D'
for ii = 1:numUE,
    
    % 'setor' � o n� do setor que UEii pertence; 'inst' � o SLOT de TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == ii);
    
    % posi��o do UEii em rela��o a BS_setor
    pos_x = real(vtUePos(ii)) - real(vtBsSetor(setor));
    pos_y = imag(vtUePos(ii)) - imag(vtBsSetor(setor));
    z_aux = pos_x + 1j*pos_y;  
    
    % angulo AZIMUTAL do UEii em rela��o a posi��o da BS_setor, em
    anguloUE180 = (180/pi)*angle(z_aux);         % ~ [-180�, +180�]
    anguloUE360 = wrapTo360(anguloUE180);        % ~ [0�, 360�]
    
    % angulo de STEERING da BS do 'setor'
    angStrTo360 = (180/pi)*vtAngST(setor);       % ~ [0�, 360�]
    angStrTo180 = wrapTo180(angStrTo360);        % ~ [-180�, +180�]
    
    % calculando a diferen�a entre o ang. de AZIMUTAL da UEii com o angulo de STEERING da BS_setor
    mtdifAngsHor_2D(setor, ii) = min(abs(anguloUE360 - angStrTo360), abs(anguloUE180 - angStrTo180));
    
    % la�o percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' � diferente do 'setor' UEii (la�o externo), ent�o
        if s ~= setor,
            
            % Posi��o do UEii em rela��o a BS,s 
            pos_xRel = real(vtUePos(ii)) - real(vtBsSetor(s));
            pos_yRel = imag(vtUePos(ii)) - imag(vtBsSetor(s));
            zRel = pos_xRel + 1j*pos_yRel;
            
            % angulo AZIMUTAL de UEii em rela��o a posi��o da BS,s
            anguloRelUE180 = (180/pi)*angle(zRel);             % ~ [-180�, +180�]
            anguloRelUE360 = wrapTo360(anguloRelUE180);        % ~ [0�, 360�]
    
            % angulo de STEERING da BS do setor_s
            angStr360 = (180/pi)*vtAngST(s);                   % ~ [0�, 360�]
            angStr180 = wrapTo180(angStr360);                  % ~ [-180�, +180�]
            
            % calculando a diferen�a entre o ang. de AZIMUTAL da UEii com o angulo de STEERING da BS_s
            mtdifAngsHor_2D(s, ii) = min(abs(anguloRelUE360 - angStr360), abs(anguloRelUE180 - angStr180));
        end
    end
end

% CALCULANDO O  PADR�O DE RADIA��O HORIZONTAL
Ah_2D = -min(12.*((mtdifAngsHor_2D./fi3dB_2D).^2), Am);

% MATRIZ de diferen�a entre o angulo ELEVA��O de cada UEs para o angulo de INCLINA��O (TILD) da BS de cada setor (em, GRAUS �)
mtdifAngsVer_2D = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% calculando os valores da matriz 'mtdifAngsVer_2D'
mtdifAngsVer_2D = abs(mtThetaUE - angDownTild_2d);

% CALCULANDO O PADR�O DE RADIA��O VERTICAL
Av_2D = -min(12.*((mtdifAngsVer_2D./theta3dB_2D).^2), SLA);    % [linhas, colunas] = [setores, UE's]

% CALCULANDO O PADR�O DE RADIA��O TOTAL
A_2D = -min(-(Ah_2D + Av_2D), Am);

% COEFICIENTES DO CANAL AO QUADRADO EM dB (SEM FAST-FADING)
H_2d = G_BS + A_2D - mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em ESCALA LINEAR (sem fast-fading)
h_2d = db2lin(H_2d);

% SINR
Y2D = [];   % SINR

% la�o percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
    
    % 'setor' que o UEjj pertence e o SLOT DE TEMPO 'inst' que o UEii est� ativo
    [setor, inst] = find(mtUeSector == ii);
    
    % calcula os coeficientes de canais INTERFERENTES p/ UEii
    aux = sum(h_2d(:, ii)) - h_2d(setor, ii);
    
    % calculo da SNIR
    Y2D(ii) = (Pot*h_2d(setor, ii))/(Pot*aux + PN);
    
end

% SINR em dB
Y2D_dB = 10*log10(Y2D);

% EFICI�NCIA ESPECTRAL DE CADA USU�RIO
R_2d = log2(1 + Y2D);

% EFICI�NCIA ESPECTRAL TOTAL
RTotal_2D = sum(R_2d, 2);

figure;             % gera uma nova figura
cdfplot(Y2D_dB)     % plota a CDF do SINR p/ 2DBF em dB 
hold on;


%% BEAMFORMING UE ESPECIFICO (3DBF UE-SPECIFIC) %%

% valores tirados do artigo
vtFi3dB_Esp = [30 20 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS] --> (Fig. 3, p. 4836)
vtTheta3dB_Esp = [10 10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS] --> (Fig. 3, p. 4836)
%vtFi3dB_Esp = [30 30 30 30 30 30 30 30 30];
%vtTheta3dB_Esp = [8:2:24];

% MATRIZ de DIFEREN�A entre o ang. ELEVA��O de cada UE para o ang. de INCLINA��O (TILD) da BS de cada setor (em, GRAUS �)
mtDifAngsVer_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngsVer_Esp'
for jj = 1:numUE,
    
    % 'setor' � o n� do setor que UEjj pertence; 'inst' � o SLOT de TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == jj);  
        
    % la�o percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' � diferente do 'setor' UEjj (la�o externo), ent�o
        if s ~= setor,
            
            % angulo de ELEVA��O do UEjj em rela��o BS_s
            angElevRelUE_Esp = mtThetaUE(s, jj);                   % ~ [�,�]
            
            % ang. de INCLINA��O (TILD) da BS_s ser� igual o ang. de Eleva��o do UE ativo no SLOT DE TEMPO 'inst' nesse setor_s
            angVertUeAtivo = mtThetaUE(s, mtUeSector(s, inst));  % = ang. de INCLINA��O da BS_s
            
            % calculando a diferen�a entre o ang. de ELEVA��O do UEjj com BS_s E o ang. de INCLINA��O (TILD) da BS_s 
            mtDifAngsVer_Esp(s, jj) = abs(angElevRelUE_Esp - angVertUeAtivo); 
        end
    end
end    


% MATRIZ de DIFEREN�A entre o angulo AZIMUTAL de cada UEs para o angulo de STEERING da BS de cada setor (em, GRAUS �)
mtDifAngHor_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngHor_Esp'
for jj = 1:numUE,
    
    % 'setor' � o n� do setor que UEjj pertence; 'inst' � o SLOT DE TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == jj);
     
    % la�o percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' � diferente do setor UEjj (la�o externo), ent�o
        if s ~= setor,
            
            % Posi��o do UEjj em rela��o a BS_s 
            pos_xRel = real(vtUePos(jj)) - real(vtBsSetor(s));
            pos_yRel = imag(vtUePos(jj)) - imag(vtBsSetor(s));
            zRel = pos_xRel + 1j*pos_yRel;
            
            % angulo AZIMUTAL de UEjj em rela��o a posi��o da BS_s, em
            anguloRelUE180 = (180/pi)*angle(zRel);             % ~ [-180�, +180�]
            anguloRelUE360 = wrapTo360(anguloRelUE180);        % ~ [0�, 360�]
            
            % Posi��o do UE_ATIVO no setor_s no instante 'inst'
            xBS_rel = real(vtUePos(mtUeSector(s, inst))) - real(vtBsSetor(s));
            yBS_rel = imag(vtUePos(mtUeSector(s, inst))) - imag(vtBsSetor(s));
            zBS_rel = xBS_rel + 1j*yBS_rel;
            
            % angulo de STEERING da BS (= ang. AZIMUTAL do UE ATIVO) do setor_s no SLOT DE TEMPO 'inst'
            angsStrBS180 = (180/pi).*angle(zBS_rel);       % [-180�, +180�]
            angsStrBS360 = wrapTo360(angsStrBS180);        % [0�, 360�]
            
            % calculando a diferen�a entre o ang. de AZIMUTAL da UEjj com o ang. de STEERING da BS_s  
            mtDifAngHor_Esp(s, jj) = min(abs(anguloRelUE360 - angsStrBS360), abs(anguloRelUE180 - angsStrBS180));    
        end 
    end 
end

% TENSORES para calcular o
Av_Esp = []; % padr�o RADIA��O VERTICAL da antena para cada ANGULO \theta_3dB
Ah_Esp = []; % padr�o RADIA��O HORIZONTAL da antena para cada angulo \phi_3dB
A_ESP = [];  % padr�o RADIA��O TOTAL p/ Beamforming ESPECIFICO   
H_ESP = [];  % COEFICIENTE do CANAL ao quadrado em dB (sem fast-fading)

% la�o percorrendo cada angulo \theta_3dB, \phi_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    % padr�o de radia��o VERTICAL
    Av_Esp(:,:,ii) = -min(12.*((mtDifAngsVer_Esp./vtTheta3dB_Esp(ii)).^2), SLA);  
    
    % padr�o de radia��o HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtDifAngHor_Esp./vtFi3dB_Esp(ii)).^2), Am);
    
    % Padr�o de RADIA��O TOTAL
    A_ESP(:,:,ii) = -min(-(Ah_Esp(:,:,ii) + Av_Esp(:,:,ii)), Am);

    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_ESP(:,:,ii) = G_BS + A_ESP(:,:,ii) - mtPL + mtNormal;   % [dB]
    
end

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR para BEAMFORMING ESPECIFICO
Y_ESP = [];   % SINR

% la�o percorrendo as dimens�es do tensor h_esp
for jj = 1:length(vtFi3dB_Esp),
    
    % la�o percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj est� ativo
        [setor, inst] = find(mtUeSector == ii);      
        
        % soma o h_esp interferentes de cada UE
        aux = sum(h_esp(:, ii, jj)) - h_esp(setor, ii, jj);
        
        % calculo da SNIR, onde cada linha ser� Y_ESP p/ Fi3dB e theta3dB
        Y_ESP(jj, ii) = (Pot*h_esp(setor, ii, jj))/(Pot*aux + PN);
    
    end
end

% SINR em dB
YESP_dB = 10*log10(Y_ESP);

% EFICI�NCIA ESPECTRAL DE CADA USU�RIO
R_esp = log2(1 + Y_ESP);

% EFICI�NCIA ESPECTRAL TOTAL
RTotal_esp = sum(R_esp, 2);

hold on;
% figure;
cdfplot(YESP_dB(1, :))
%cdfplot(YESP_dB(2, :))
%cdfplot(YESP_dB(3, :))
%legend('Conventional', 'UE especifica - (\theta_{3dB}, \phi_{3dB) = (70�, 10�)}');


%% BEAMFORMING GRUPO-ESPECIFICO (3DBF UE-GROUP) %%

% valores tirados do artigo 
%vtFi3dB_gr = [30 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS] ---> (Fig. 3, p. 4836)
%vtTheta3dB_gr = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS] ---> (Fig. 3, p. 4836)
Fi3dB_gr = 30;
Theta3dB_gr = 10;
% B = 32;                            % numero de padr�es de feixes (ou, GRUPOS) = Bh*Bv
% Bh = 8;                            % numero de feixes HORIZONTAIS p/ cada setor
% Bv = B/Bh;                         % numero de feixes VERTICAIS p/ cada setor

%vtBh = [8 16 32 32 32 64 64 64 64 128 128 128 128 256 256 256 256 512 512 512 512 1024 1024 1024];
%vtBv = [2  2  2  4  8  2  4  8 16   2   4   8  16   2   4   8  16   2   4   8 16     2    4    8];
vtBh = [8 16 32 64 64 128 128];
vtBv = [2  2  2  2  4 2  4];

particoes = length(vtBh);

% matriz DIFEREN�A do ang. AZIMUTAL de cada UE para o ang. STEERING da BS de cada setor
mtDifAngsHor_gr = zeros(S, numUE, particoes);  % [linha, coluna] = [Setor, UE] 

% diferen�a do ang de ELEVA��O de cada UE p/ o ang. de INCLINA��O da BS de cada setor (em, GRAUS �)
mtdifAngsVer_gr = zeros(S, numUE, particoes);

for ii = 1:particoes, 
    
    [mtAngsStering, mtAngsTild_gr] = geracaoGrupos(vtBh(ii), vtBv(ii), FatorSetor, M);

    % la�o percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsHor_gr'
    for jj = 1:numUE,

        % 'setor' que o UEjj pertence e o slot_de_tempo 'inst' que o UEjj est� ativo 
        [setor, inst] = find(mtUeSector == jj);    

        % posi��o do UEjj em rela��o a BS_setor na qual ele pertence
        pos_x = real(vtUePos(jj)) - real(vtBsSetor(setor));
        pos_y = imag(vtUePos(jj)) - imag(vtBsSetor(setor));
        z_aux = pos_x + 1j*pos_y;

        % angulo AZIMUTAL do UEjj em rela��o a posi��o da BS_setor na qual ele pertence, em
        anguloUE180 = (180/pi).*angle(z_aux);  % ~ [-180�, +180�]
        anguloUE360 = wrapTo360(anguloUE180);  % ~ [0�, +360�]

        % GRUPO de angs de STEERING p/ BS_setor na qual UEjj faz parte, em:
        angFiGrupo_360 = mtAngsStering(setor,:);     % ~ [0�, 360] 
        angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % ~ [-180�, +180�]

        % calculando a diferen�a entre o ang. AZIMUTAL do UEjj (em rela��o a BS_setor na qual UEjj pertence) e ang. STEERING da BS_setor 
        dif1 = min(abs(anguloUE360 - angFiGrupo_360));
        dif2 = min(abs(anguloUE180 - angFiGrupo_180));
        difAngulo = min(dif1, dif2);
        mtDifAngsHor_gr(setor, jj, ii) = difAngulo;

        % la�o percorrendo os SETORES 
        for s = 1:S,

            % se o setor 's' � diferente do setor do UEjj (la�o externo), ent�o 
            if s ~= setor,

                % Posi��o do UEjj em rela��o a BS_s
                pos_xRel = real(vtUePos(jj)) - real(vtBsSetor(s));
                pos_yRel = imag(vtUePos(jj)) - imag(vtBsSetor(s));
                zRel = pos_xRel + 1j*pos_yRel;

                % angulo AZIMUTAL de UEjj em rela��o a BS_s, em:
                anguloRelUE180 = (180/pi)*angle(zRel);             % ~ [-180�, +180�]
                anguloRelUE360 = wrapTo360(anguloRelUE180);        % ~ [0�, 360�]

                % Posi��o do UE ATIVO no setor_s no instante 'inst'
                xBS_rel = real(vtUePos(mtUeSector(s, inst))) - real(vtBsSetor(s));
                yBS_rel = imag(vtUePos(mtUeSector(s, inst))) - imag(vtBsSetor(s));
                zBS_rel = xBS_rel + 1j*yBS_rel;

                % ang. AZIMUTAL do UE_ATIVO no slot_de_tempo 'inst' no setor_s, em:
                angAzimUE_ATIVO180 = (180/pi).*angle(zBS_rel); 
                angAzimUE_ATIVO360 = wrapTo360(angAzimUE_ATIVO180);

                % encontrando qual ang. de STEEREING da BS_s dentro do angs. de STEERING dispon�veis p/ setor_s 
                [diferenca, indice] = min(abs(angAzimUE_ATIVO360 - mtAngsStering(s,:)));

                % calculando a diferen�a entre o ang. AZIMUTAL do UEjj (em rela��o a BS_s) e o ang. de STEERING da BS_s, em
                ddif1 = min(abs(anguloRelUE360 - mtAngsStering(s,indice)));             % ~ [0�, 360�]
                ddif2 = min(abs(anguloRelUE180 - wrapTo180(mtAngsStering(s,indice))));  % ~ [-180�, +180�]
                mtDifAngsHor_gr(s, jj, ii) = min(ddif1, ddif2);
            end
        end
    end
    
    % la�o percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsVer_gr'
    for jj= 1:numUE,

        % 'setor' � o n� do setor que o UEjj pertence; 'inst' � o slot_de_tempo na qual UEjj est� ativo
        [setor, inst] = find(mtUeSector == jj);       

        % angulo de ELEVA��O do UEjj em rela��o BS_setor, na qual ele pertence
        angElevUE = mtThetaUE(setor, jj); % em, graus (�)

        % calculando a DIFEREN�A entre o ang. de ELEVA��O do UEjj e o Angs de INCLIN�O (TILD) dispon�veis p/ BS_setor na qual ele pertence 
        mtdifAngsVer_gr(setor, jj, ii) = min(abs(angElevUE - mtAngsTild_gr(setor,:)));

        % la�o percorrendo todos os  SETORES
        for s = 1:S,

            % se o setor 's' � diferente do 'setor' que cont�m o UEjj (la�o externo), ent�o
            if s ~= setor,

                % ang. de ELEVA��O do UEjj em RELA��O a BS_s
                angElevUeRel = mtThetaUE(s,jj);

                % ang. de ELEVA��O do UE_ATIVO no setor 's' no slot_de_tempo 'inst'
                angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));

                % encontrando o ang. de ELEVA��O (TILD) da BS_s dentro do GRUPO de angs. de INCLINA��O d�sponiveis para o setor_s
                [diferenca, indice] = min(abs(angElevUeAtivo - mtAngsTild_gr(s,:)));

                % calculando a diferen�a entre o ang. ELEVA��O do UEjj (em rela��o a BS_s) e o ang. de INCLINA��O da BS_s, em
                mtdifAngsVer_gr(s, jj, ii) = min(abs(angElevUeRel - mtAngsTild_gr(s, indice)));    
            end
        end    
    end

end

% TENSORES para calcular o padr�o de RADIA��O 
Ah_gr = [];  % HORIZONTAL da antena para cada angulo \phi_3dB 
Av_gr = [];  % VERTICAL da antena para cada angulo \theta_3dB
A_GR = [];   % TOTAL p/ Beamforming GRUPO

% la�o percorrendo cada angulo \theta_3dB, \phi_3dB p/ calcular cada dimens�o do setor
for ii = 1:particoes,
    
    % Ah para cada \fi_3dB
    Ah_gr(:,:,ii) = -min(12.*((mtDifAngsHor_gr(:,:,ii)./Fi3dB_gr).^2), Am);

    % Av para cada \theta_3dB
    Av_gr(:,:,ii) = -min(12.*((mtdifAngsVer_gr(:,:,ii)./Theta3dB_gr).^2), SLA);

    % Padr�o de radia��o TOTAL
    A_GR(:,:,ii) = -min(-(Av_gr(:,:,ii) + Ah_gr(:,:,ii)), Am);

    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_GR(:,:,ii) = G_BS + A_GR(:,:,ii) - mtPL + mtNormal;   % [dB]
end

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_gr = db2lin(H_GR);

% SINR para Beamforming GRUPO
Y_GR = [];   % SINR

% la�o percorrendo cada dimensao do tensor p/ calcular cada dimens�o do setor
for jj = 1:particoes,
    
    % la�o percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj est� ativo
        [setor, inst] = find(mtUeSector == ii);
        
        % soma h_gr interferentes para cada UE 
        aux = sum(h_gr(:, ii, jj)) - h_gr(setor, ii, jj);
        
        % calculo da SNIR, onde cada linha ser� Y_ESP p/ Fi3dB e theta3dB
        Y_GR(jj, ii) = (Pot*h_gr(setor, ii, jj))/(Pot*aux + PN);     % linha: SNIR p/ cada \theta_3dB e \fi_3dB
        
    end
end

% SINR em dB
YGR_dB = 10.*log10(Y_GR);
cdfplot(YGR_dB(1,:))
cdfplot(YGR_dB(2,:))
cdfplot(YGR_dB(3,:))
cdfplot(YGR_dB(4,:))
cdfplot(YGR_dB(5,:))
cdfplot(YGR_dB(6,:))
cdfplot(YGR_dB(7,:))
legend('2DBF', '3DBF usu�rio', '3DBF 16 grupo', '3DBF 32 grupo', '3DBF 64 grupo', '3DBF 128 grupo', '3DBF 256 grupo');
xlabel('SINR (dB)')
ylabel('CDF')
title('')
