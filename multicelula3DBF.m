%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture
clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 7;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setoriza��o, i.e, setores/celulas
S = M*FatorSetor;                                    % n�mero de setores. S = {1, 2, 3, ..., }      
UEcadaSetor = 5;                                     % numero de UE's por (micro)setor
numUE = UEcadaSetor*S;                               % numero de UE's total
R = 80;                                              % raio da pequena celula
xBS = 0;                                             % Posi��o do eixo x da BS
yBS = 0;                                             % Posi��o do eixo y da BS
% vtSector = [ R*exp( 1j*[0 2*pi/3 4*pi/3] ) ];        % vetor marca��o dos pontos de sectoriza��o
bandWidth = 10e+6;                                   % largura de banda 
fc = 3.5;                                            % frequencia de portadora [Ghz] na small-cell
Am = 25;                                             % Atenua��o Maxima [dB]
SLA = 20;                                            % Limite de n�vel de lobulo-lateral [dB]
G_BS = 5;                                            % Ganha da antena BS de pequena celula [dBi]
P_BS = 24;                                           % Potencia de TX da BS por setor [dBm]
sigma = 4;                                           % desvio padr�o do vetor de sombreamento.[dB]
No = -174;                                           % densidade espectral de potencia [ dBm/Hz ]
FN = 5;                                              % figura de ruido [5 dB]
H_BS = 10;                                           % altura da antena da Esta��o-Base [metros] 
H_UE = 1.5;                                          % altura da antena da esta��o-m�vel [metros]


%% GERANDO A POSI��O DA BS DE CADA C�LULA %%

% VALORES DE TESTE
vtBS = [0*exp(1j*0), 2*R*exp(1j*0), 2*R*exp(1j*pi/3), 2*R*exp(1j*2*pi/3), 2*R*exp(-j*pi), 2*R*exp(-1j*2*pi/3) 2*R*exp(-1j*pi/3)];

% sobrepor gr�ficos
hold on

% la�o percorrendo cada BS p/ plotar suas posi��es e sua �REA DE COBERTURA
for ii = 1:length(vtBS)
    circle(real(vtBS(ii)),imag(vtBS(ii)),R)
    plot(real(vtBS(ii)),imag(vtBS(ii)),'xr')
end

% vtBS = [];                                          % vetor com a posi��o da BS (esta��o base) de cada celula
% raio_bs = 0;                                        % raio das BS's - gerar aleatorio
% ang_bs = 0;                                         % angulos de posi��es das BS's - gerar aleatorio
% vtBS(1) = raio_bs.*exp(1j*ang_bs);                  % primeira BS ser� na origem do sistema
% vtBS = [0, 2*R*exp(1j*0)]

% la�o que percorre o n�mero de c�lulas p/ gerar a posi��o das BS 
% for ii = 2:M,
%     
%     while true, 
%         
%         % calculando a posi��o da BSii 
%         vtBS(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
%         
%         if ,
%         end
%     end
%     
% end

vtBsSetor = repelem(vtBS, FatorSetor);             % vetor posi��o da BS's p/ cada setor 


%% GERANDO A POSI��O DE CADA UE ~ U(10,R) %%

% vetor p/ as posi��es aleat�rios dos UE's 
vtUePos = []; 

% vetor com o angulo Horizontal de inicio (em, graus �) de cada setor 
vtAngIncSetor = repmat([0, 120, 240], 1, M);

% la�o percorre o numero de UE's p/ gerar a posi��o de cada UE de forma que teremos um UE ATIVO em cada setor por SLOT de TEMPO
for ii = 1:numUE,
    
    while true,
        
        % calculando a posi��o do UEii 
        vtUePos(ii) = (6*R*rand - 3*R) + 1j*(6*R*rand - 3*R);
        
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
% hold off;


%% CONSTRUINDO A MATRIZ INDICE DO UE'S POR SETORES (linhas) %%

mtUeSector = zeros(S, UEcadaSetor);  % [linhas, colunas] = [setores, indice do UE de cada setor ativo no SLOT de TEMPO]

% la�o percorrendo cada UE 
for ii = 1:numUE,
    
    % dif = valor da menor diferen�a entre o ang. Azimutal do UE e posi��o da BS's de cada celula
    % indCel = indice da BS da celula na qual UEii pertence
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
end


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
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')


%% DADOS COMUNS PARA OS BEAMFORMINGS %%

% vetor de DESVANECIMENTO R�PIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];

% matriz de SHADOWING normal em dB
mtNormal = sigma.*randn(S,numUE);    % [Linhas, Colunas] = [Setores, UE's]

% matriz de Angulos de ELEVA��O [em, GRAUS (�)] de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [LINHA, COLUNA] = [SETOR, UE]

% Pot�ncia do Ru�do Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Pot�ncia da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% angulos de BORESIGHT FIXO p/ BS de cada setor, i.e, angulo azimutal na qual teremos o ganho m�ximo da antena 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de boresight p/ cada celula [Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS] -> P�g. 4837, Simulation Setup
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS] -> P�g. 4837, Simulation Setup
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS] -> P�g. 4836, Simulation Setup

% MATRIZ de diferen�a entre o angulo ELEVA��O de cada UEs para o angulo de INCLINA��O (TILD) da BS de cada setor (em, GRAUS �)
mtdifAngsHor_2D = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtdifAngsHor_2D'
for ii = 1:numUE,
    
    % 'setor' � o n� do setor que UEii pertence; 'inst' � o SLOT de TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == ii);
    
    % posi��o do UEii em rela��o a BS do setor
    pos_x = real(vtUePos(ii)) - real(vtBsSetor(setor));
    pos_y = imag(vtUePos(ii)) - imag(vtBsSetor(setor));
    z_aux = pos_x + 1j*pos_y;  
    
    % angulo do UEii em rela��o a posi��o da BS,setor
    anguloUE180 = (180/pi)*angle(z_aux);         % ~ [-180�, +180�]
    anguloUE360 = wrapTo360(anguloUE180);        % ~ [0�, 360�]
    
    % angulo de STEERING da BS do 'setor'
    angStrTo360 = (180/pi)*vtAngST(setor);       % ~ [0�, 360�]
    angStrTo180 = wrapTo180(angStrTo360);        % ~ [-180�, +180�]
    
    % calculando a diferen�a entre o ang. de AZIMUTAL da UEii com o angulo de STEERING da BS_setor
    mtdifAngsHor_2D(setor, ii) = min(abs(anguloUE360 - angStrTo360), abs(anguloUE180 - angStrTo180));
    
    % la�o percorrendo cada angulo de STEERING da c�lula
    for s = 1:S,
        
        % se o setor 's' � diferente do 'setor' UEii (la�o externo), ent�o
        if s ~= setor,
            
            % Posi��o do UEii em rela��o a BS,s 
            pos_xRel = real(vtUePos(ii)) - real(vtBsSetor(s));
            pos_yRel = imag(vtUePos(ii)) - imag(vtBsSetor(s));
            zRel = pos_xRel + 1j*pos_yRel;
            
            % angulo do UEii em rela��o a posi��o da BS,s
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



% MATRIZ de diferen�a entre o angulo AZIMUTAL de cada UEs para o angulo de BORESIGHT da BS de cada setor (em, GRAUS �)
mtdifAngsVer_2D = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% calculando os valores da matriz 'mtdifAngsVer_2D'
mtdifAngsVer_2D = abs(mtThetaUE - angDownTild_2d);

% CALCULANDO O PADR�O DE RADIA��O VERTICAL
Av_2D = -min(12.*((mtdifAngsVer_2D./theta3dB_2D).^2), SLA);    % [linhas, colunas] = [setores, UE's]




% CALCULANDO O PADR�O DE RADIA��O TOTAL
A_2D = -min(-(Ah_2D + Av_2D), Am);

% COEFICIENTES DO CANAL AO QUADRADO EM dB (SEM FAST-FADING)
H_2d = G_BS + A_2D + mtPL + mtNormal;   % [dB]

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
    
figure;             % gera uma nova figura
cdfplot(Y2D_dB)     % plota a CDF do SINR p/ 2DBF em dB 
hold on;


%% BEAMFORMING UE ESPECIFICO %%

% valores tirados do artigo
vtFi3dB_Esp = [30 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS] --> (Fig. 3, p. 4836)
vtTheta3dB_Esp = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS] --> (Fig. 3, p. 4836)

% matriz de diferen�a entre o angulo ELEVA��O de cada UEs para o angulo de INCLINA��O (TILD) da BS de cada setor (em, GRAUS �)
mtDifAngsVer_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngsVer_Esp'
for jj = 1:numUE,
    
    % 'setor' � o n� do setor que UEjj pertence; 'inst' � o SLOT de TEMPO que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == jj);  
    
    % se temos algum UE INATIVO no slot de tempo 'inst', ent�o pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % angulo de ELEVA��O do UEjj ser� igual o angulo de INCLINA��O do BS_setor:
    angElevUE_Esp = mtThetaUE(setor, jj); % ~ [�,�]
        
    % la�o percorrendo cada SETOR
    for s = 1:S,
        % se o setor 's' � diferente do 'setor' UEjj (la�o externo), ent�o
        if s ~= setor,
            angVertUeAtivo = mtThetaUE(s, mtUeSector(s, inst)); % ang. de INCLINA��O da BS_s sera igual o ang. de Eleva��o do UE ativo nesse setor_s 
            mtDifAngsVer_Esp(s, jj) = abs(angElevUE_Esp - angVertUeAtivo); 
        end
    end
end    


% matriz de diferen�a entre o angulo AZIMUTAL de cada UEs para o angulo de STEERING da BS de cada setor (em, GRAUS �)
mtdifAngHor_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% la�o percorrendo cada UE's p/ calcular os valores da matriz 'mtdifAngHor_Esp'
for jj = 1:numUE,
    
    % angulo azimutal do UE, em
    angUE180 = (180/pi)*angle(vtUePos(jj));       % ~ [-180�, 180�]
    angUE360 = wrapTo360(angUE180);               % ~ [0�, 360�] 
    
    % 'setor' � o n� do setor que UEjj pertence; 'inst' � o instante que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == jj);       
    
    % se temos algum UE INATIVO no slot de tempo 'inst', ent�o pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % la�o percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' � diferente do setor UEjj (la�o externo), ent�o
        if s ~= setor,
            % angulo horizontal do UE ativo no instante 'inst' no setor 'setor'
            angUeAtivo = wrapTo360((180/pi).*angle(vtUePos(mtUeSector(s, inst)))); % [0�, 360�]
            
            % calcular a diferen�a entre o angulos horizontais em
            dif1 = min(abs(angUE360-angUeAtivo));  % ~[0�, 360�]
            dif2 = min(abs(angUE180-wrapTo180(angUeAtivo))); % ~[-180�, +180�]
            
            mtdifAngHor_Esp(s, jj) = min(dif1, dif2);    
        end 
    end 
end

% TENSORES para calcular 
Av_Esp = []; % padr�o RADIA��O VERTICAL da antena para cada ANGULO \theta_3dB
Ah_Esp = []; % padr�o RADIA��O HORIZONTAL da antena para cada angulo \phi_3dB
A_ESP = [];  % padr�o RADIA��O TOTAL p/ Beamforming ESPECIFICO   
H_ESP = [];  % COEFICIENTE do CANAL ao quadrado em dB (sem fast-fading)

% la�o percorrendo cada angulo \theta_3dB, \phi_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    % padr�o de radia��o VERTICAL
    Av_Esp(:,:,ii) = -min(12.*((mtDifAngsVer_Esp./vtTheta3dB_Esp(ii)).^2), SLA);  
    
    % padr�o de radia��o HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtdifAngHor_Esp./vtFi3dB_Esp(ii)).^2), Am);
    
    % Padr�o de RADIA��O TOTAL
    A_ESP(:,:,ii) = -min(-(Ah_Esp(:,:,ii) + Av_Esp(:,:,ii)), Am);

    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_ESP(:,:,ii) = G_BS + A_ESP(:,:,ii) + mtPL + mtNormal;   % [dB]
    
end

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR para Beamforming ESPECIFICO
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

hold on;
% figure;
cdfplot(YESP_dB(1, :))
%cdfplot(YESP_dB(2, :))
%cdfplot(YESP_dB(3, :))
%legend('Conventional', 'UE especifica - (\theta_{3dB}, \phi_{3dB) = (70�, 10�)}');


%% BEAMFORMING GRUPO-ESPECIFICO

% valores tirados do artigo 
vtFi3dB_gr = [70 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_gr = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]
B = 16;                            % numero de padr�es de feixes (ou, GRUPOS)
Bh = 8;                            % numero de feixes horizontais p/ cada setor
Bv = B/Bh;                         % numero de feixes vertiais p/ cada setor

% mtDirFeixes = zeros(Bh, Bv);     % matriz de dire��o de feixes, onde cada elemento: (Fi_st, Theta_tild)
% mtDirFeixes = cell(Bh, Bv);      % celula de dire��o de feixes

angsFiSt_gr = zeros(1, Bh);
andDownTild_gr = zeros(1, Bv);

angHorInic = 0;
angHorFinal = 360;
passo = (angHorFinal - angHorInic)/(Bh*FatorSetor);

angsFiSt_gr = linspace(angHorInic + (passo/2), angHorFinal - (passo/2), Bh*FatorSetor);

% matriz de angulos \fi_st para uma celula
mtAngsFiSt_gr = reshape(angsFiSt_gr, Bh , FatorSetor*M);   
mtAngsFiSt_gr = mtAngsFiSt_gr';                  % [linha, coluna] = [setor, angulo Thi_st p/ cada grupo do setor]

% angDownTild_gr = linspace(-90, 90, Bv);
angDownTild_gr = [4 8];

% matriz de angulos \theta_thild p/ uma celula;
mtAngsThetaTild_gr = repmat(angDownTild_gr, FatorSetor*M, 1); % linha: setor, coluna: \theta_tild 's p/ cada setor

% matriz DIFEREN�A do ang. AZIMUTAL de cada UE para o ang. STEERING (\phi_st) da BS de cada setor
mtdifAngsHor_gr = zeros(S, numUE);  % [linha, coluna] = [Setor, UE] 

% la�o percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsHor_gr'
for jj = 1:numUE,
    
    % 'setor' que o UEjj pertence e o slot_de_tempo 'inst' que o UEjj est� ativo 
    [setor, inst] = find(mtUeSector == jj);    
    
    % se temos algum UE INATIVO no slot de tempo 'inst', ent�o pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % calcula a diferen�a entre angulos
    angFiGrupo_360 = mtAngsFiSt_gr(setor,:);     % GRUPOS de angulos de steering [0�, 360] disponiveis para o setor_'setor' onde UEjj pertence
    angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % GRUPOS de angulos de steering [-180�, +180�] disponiveis para o setor_'setor' onde UEjj pertence
    
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  % angulo azimutal do UE atual em [-180, 180]
    anguloUE360 = wrapTo360(anguloUE180);        % angulo azimutal do UE atual em [0, 360]
   
    dif1 = min(abs(anguloUE360 - angFiGrupo_360));
    dif2 = min(abs(anguloUE180 - angFiGrupo_180));
    difAngulo = min(dif1, dif2);
    
    mtdifAngsHor_gr(setor, jj) = difAngulo;
    
    % la�o percorrendo os SETORES 
    for s = 1:S, 
        % se o setor 's' � diferente do setor do UEjj (la�o externo), ent�o 
        if s ~= setor,
            angUeAtivo = wrapTo360((180/pi).*angle(vtUePos(mtUeSector(s, inst)))); % UE ativo no instante 'inst' no setor 'setor'
            [M, II] = min(abs(angUeAtivo - mtAngsFiSt_gr(s,:)));
            ddif1 = min(abs(anguloUE360 - mtAngsFiSt_gr(s,II)));
            ddif2 = min(abs(anguloUE180 - wrapTo180(mtAngsFiSt_gr(s,II))));
            mtdifAngsHor_gr(s, jj) = min(ddif1, ddif2);
        end
    end
end


% diferen�a de angulos de eleva��o de cada UE de grupo
mtdifAngsVer_gr = zeros(S, numUE);

% la�o percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsVer_gr'
for jj= 1:numUE,
    
    % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj est� ativo
    [setor, inst] = find(mtUeSector == jj);       
    
    % se temos algum UE INATIVO no slot de tempo 'inst', ent�o pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % calculo o angulo de eleva��o THETA do UE
    angElevUE180 = mtThetaUE(setor, jj);
    
    mtdifAngsVer_gr(setor, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(setor,:)));
    
    % la�o percorrendo todos os  SETORES
    for s = 1:S,
        % Percorrendo os OUTROS setores, i.e, setores diferente do setor que cont�m o UEjj (la�o externo)
        if s ~= setor,
            angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));    % ang. de ELEVA��O do UE ativo no setor 's' no instante 'inst'
            [M, II] = min(abs(angElevUeAtivo - mtAngsThetaTild_gr(s,:)));
            mtdifAngsVer_gr(s, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(s, II)));    
        end
    end    
end

% TENSORES para calcular o padr�o de RADIA��O 
Ah_gr = [];  % HORIZONTAL da antena para cada angulo \phi_3dB 
Av_gr = [];  % VERTICAL da antena para cada angulo \theta_3dB
A_GR = [];   % TOTAL p/ Beamforming GRUPO

% la�o percorrendo cada angulo \theta_3dB, \phi_3dB p/ calcular cada dimens�o do setor
for ii = 1:length(vtTheta3dB_Esp),

    % Ah para cada \fi_3dB
    Ah_gr(:,:,ii) = -min(12.*((mtdifAngsHor_gr./vtFi3dB_gr(ii)).^2), Am);
    
    % Av para cada \theta_3dB
    Av_gr(:,:,ii) = -min(12.*((mtdifAngsVer_gr./vtTheta3dB_gr(ii)).^2), SLA);
    
    % Padr�o de radia��o TOTAL
    A_GR(:,:,ii) = -min(-(Av_gr(:,:,ii) + Ah_gr(:,:,ii)), Am);
    
    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_GR(:,:,ii) = G_BS + A_GR(:,:,ii) + mtPL + mtNormal;   % [dB]
end


% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_gr = db2lin(H_GR);

% SINR para Beamforming GRUPO
Y_GR = [];   % SINR

% la�o percorrendo cada dimensao do tensor p/ calcular cada dimens�o do setor
for jj = 1:length(vtFi3dB_Esp),
  
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
YGR_dB = 10*log10(Y_GR);
cdfplot(YGR_dB(1,:))
legend('Conventional', 'UE especifica', '16 grupo');


