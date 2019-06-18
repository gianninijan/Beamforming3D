%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture
clear all;
clc;
close all;

%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setorização, i.e, setores/celulas
S = M*FatorSetor;                                    % número de setores. S = {1, 2, 3, ..., }      
UEcadaSetor = 300;                                   % numero de UE's por (micro)setor
numUE = UEcadaSetor*S;                               % numero de UE's total
R = 250;                                             % raio da pequena celula
xBS = 0;                                             % Posição do eixo x da BS
yBS = 0;                                             % Posição do eixo y da BS
vtSector = [ R*exp( 1j*[0 2*pi/3 4*pi/3] ) ];        % vetor marcação dos pontos de sectorização
bandWidth = 10e+6;                                   % largura de banda 
fc = 3.5;                                            % frequencia de portadora [Ghz] na small-cell
Am = 25;                                             % Atenuação Maxima [dB]
SLA = 20;                                            % Limite de nível de lobulo-lateral [dB]
G_BS = 5;                                            % Ganha da antena BS de pequena celula [dBi]
P_BS = 24;                                           % Potencia de TX da BS por setor [dBm]
sigma = 4;                                           % desvio padrão do vetor de sombreamento.[dB]
No = -174;                                           % densidade espectral de potencia [ dBm/Hz ]
FN = 5;                                              % figura de ruido [5 dB]
H_BS = 10;                                           % altura da antena da Estação-Base [metros] 
H_UE = 1.5;                                          % altura da antena da estação-móvel [metros]


%% POSIÇÕES DA BS E UE'S %%

% GERANDO AS POSIÇÕES DO BS's
raio_bs = 0;                                        % raio das BS's - gerar aleatorio
ang_bs = 0;                                         % angulos de posições das BS's - gerar aleatorio
posBS = raio_bs.*exp(1j*ang_bs);                    % forma exponencial das BS's
vtBsSetor = repelem(posBS, FatorSetor);             % vetor posição da BS's p/ cada setor

% GERANDO AS POSIÇÕES DO UE's ~ U(10, R) DE FORMA QUE SEJA UE ATIVO POR SETOR EM CADA SLOT DE TEMPO
vtUePos = [];   % vetor p/ as posições aleatórios dos UE's 

% vetor com o angulo Horizontal (em, graus º) de cada setor
vtAngIncSetor = repmat([0, 120, 240], 1, UEcadaSetor);

% laço percorre o numero de UE's p/ gerar sua posição
for ii = 1:numUE,
    
    while true,
        % calculando a posição do UEii 
        vtUePos(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
        
        % angulo do UE em graus
        anguloUE = wrapTo360((180/pi)*angle(vtUePos(ii))); % ~ [0º, 360º]
        
        % indice
        ind = mod(ii, FatorSetor);
        if ind == 0, 
            ind = FatorSetor;
        end
        
        % distancia do UE p/ BS do setor
        distUE = abs(vtUePos(ii) - vtBsSetor(ind));
        
        % angulo do UE em relação a posição da BS

        % se distancia do UE p/ BS do setor for menor que o Raio e maior que 10, então:
        if (distUE < R) && (distUE > 10) &&  (anguloUE >= vtAngIncSetor(ind))  && (anguloUE < (vtAngIncSetor(ind) + (360/FatorSetor))),
            break;
        end
    end
end

% plotando as POSIÇÕES da BS e UE's    
hold on;
circle(xBS,yBS,R);    
plot(xBS,yBS,'*r','MarkerSize',16);       % plotando a BS  (*vermelho)
plot(vtUePos,'*b','MarkerSize',12);       % plotando os UE's (+azul)
grid on;

% gerando as LINHAS dos SETORES. 
x = zeros(1,2*length(vtSector));
y = zeros(1,2*length(vtSector));
x([2:2:2*length(vtSector)]) = real(vtSector);
y([2:2:2*length(vtSector)]) = imag(vtSector);
plot(x,y,'y');
legend('Celula', 'Estação base', 'Usuários','Setores')

% Encontrando os indices dos UE's de cada sector.
angUE = wrapTo2Pi(angle(vtUePos));                            % angUE ~ [0, 2*pi] radianos
sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que estão no setor 1
sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que estão no setor 2
sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que estão no setor 3

% criando a MATRIZ de SETORES (linhas) e indice do UE's de cada setor 
maximo = max([length(sector1), length(sector2), length(sector3)]);
mtUeSector = max(S, maximo);                % [linhas, colunas] = [setores, slot_de_tempo]
mtUeSector(1, 1:length(sector1)) = sector1;
mtUeSector(2, 1:length(sector2)) = sector2;
mtUeSector(3, 1:length(sector3)) = sector3;

% Testando o codigo acima dos setores: 
% resultado = vtUePos(sector2);          % captura as distancia dos UE's que estão dentro do setor 1
% plot(resultado,'ro');                  % plota o valores do 'resultados' com um circulo azul  
% 
% hold off


%% CALCULAR O PATH-LOSS DOS UE'S

% MATRIZ de DISTÂNCIA de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                            

% laço percorrendo cada setor (coluna)
for ii = 1:S
    
    % laço percorrendo os UE's (linha)
    for jj = 1:numUE
        
        % calculando a DISTÂNCIA de UEjj para cada BSii
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % DIM( mtDist ) = [Linhas, Colunas] = [Setores, UE's]
        
    end
    
end


% matriz de PATH LOSS de cada UE para cada BS 
mtPL = PathLoss(mtDist, fc);                   % [linha, coluna] = [SETORES, UE's]
PL = min(mtPL);                                % valor minimo

vtDistUEtoBS = min(mtDist);
%[vtDistUEtoBS, POS] = min(mtDist, [], 1);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Perda de Caminho ')


%% DADOS COMUNS PARA OS BEAMFORMINGS

% vetor de DESVANECIMENTO RÁPIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];

% matriz de SHADOWING normal em dB
mtNormal = sigma.*randn(S,numUE);    % dim( mtNormal ) = [Linhas, Colunas] = [Setores, Usuários]

% matriz de Angulos de ELEVAÇÃO [em, GRAUS (º)] de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [LINHA, COLUNA] = [SETOR, UE]

% Potência do Ruído Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Potência da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% angulos de BORESIGHT FIXO p/ a antena de cada setor, i.e, angulo azimutal na qual teremos o ganho máximo da antena 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de boresight p/ cada celula [Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS] -> Pág. 4837, Simulation Setup
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS] -> Pág. 4837, Simulation Setup
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS] -> Pág. 4836, Simulation Setup

% calculando o padrão de radiação VERTICAL da antena
% matriz DIFERENÇA entre o ang. AZIMUTAL de cada UE p/ o ang. de STEERING da BS de cada setor
mtdifAngsVer_2D = abs(mtThetaUE - angDownTild_2d);

Av_2D = -min(12.*((mtdifAngsVer_2D./theta3dB_2D).^2), SLA);    % [linhas, colunas] = [setores, UE's]
 

% calculando o padrão de radiação HORIZONTAL da antena
% matriz DIFERENÇA entre o ang. AZIMUTAL de cada UE p/ o ang. de STEERING da BS de cada setor
mtdifAngsHor_2D= [];   
 
% laço percorrendo cada UE p/ calcular os valores da matriz 'mtdifAngsHor_2D'
for jj = 1:numUE,
    
    % angulo de cada UE ~ [-180º, 180º]
    anguloUE180 = (180/pi).*angle(vtUePos(jj)); % em, graus (º)   
   
    % angulo de cada UE ~ [0º, 360º]
    anguloUE360 = wrapTo360(anguloUE180); % em, graus (º)
    
    % angulo de steering ~ [0º, 360º]
    angStrTo360 = (180/pi).*vtAngST;    % em, graus (º)
    
    % angulo de steering ~ [-180º, 180º]
    angStrTo180 = wrapTo180(angStrTo360); % em, graus (º)

    % laço percorrendo cada setor
    for ii = 1:S,
        
        % calcula a diferença entre angulos
        dif1 = abs(anguloUE360 - angStrTo360(ii));
        dif2 = abs(anguloUE180 - angStrTo180(ii));
        difAngulo = min(dif1, dif2);
        
        mtdifAngsHor_2D(ii, jj) = difAngulo;
    end
end

% calculando padrão de radiação HORIZONTAL
Ah_2D = -min(12.*((mtdifAngsHor_2D./fi3dB_2D).^2), Am);

% Padrão de radiação TOTAL
A_2D = -min(-(Ah_2D + Av_2D), Am);

% COEFICIENTES do CANAL ao quadrado em dB (sem fast-fading)
H_2d = G_BS + A_2D - mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_2d = db2lin(H_2d);

% SINR
Y2D = [];   % SINR

% capturando o setor de cada UE atraves do INDICES
[X,I] = min(sqrt(mtdifAngsHor_2D.^2), [], 1); % X = menores elemento de cada coluna; I = posição do menor elemento na coluna, i.e, localização na linha do menor elemento de cada coluna

% laço percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
    
    % calcula os coeficientes de canais INTERFERENTES p/ UEii
    aux = sum(h_2d(:, ii)) - h_2d(I(ii), ii);
    
    % calculo da SNIR
    Y2D(ii) = (Pot*h_2d(I(ii), ii))/(Pot*aux + PN);
    
end

% SINR em dB
Y2D_dB = 10*log10(Y2D);


%% BEAMFORMING UE ESPECIFICO %%

% valores tirados do artigo
vtFi3dB_Esp = [30 20 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_Esp = [10 10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]

% matriz de diferença entre o angulo ELEVAÇÃO de cada UEs para o angulo de INCLINAÇÃO (TILD) da BS de cada setor (em, GRAUS º)
mtDifAngsVer_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% laço percorrendo cada UE's p/ calcular os valores da matriz 'mtDifAngsVer_Esp'
for jj = 1:numUE,
    
    % 'setor' é o nº do setor que UEjj pertence; 'inst' é o SLOT de TEMPO que o UEjj está ativo
    [setor, inst] = find(mtUeSector == jj);  
    
    % se temos algum UE INATIVO no slot de tempo 'inst', então pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
        
    % laço percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' é diferente do 'setor' UEjj (laço externo), então
        if s ~= setor,
            
            % angulo de ELEVAÇÂO do UEjj em relação a BS_s:
            angElevUE_Esp = mtThetaUE(s, jj); % ~ [º,º]
            
            % ang. de INCLINAÇÃO (TILD) da BS_s será igual o ang. de Elevação do UE ativo no SLOT DE TEMPO 'inst' nesse setor_s
            angVertUeAtivo = mtThetaUE(s, mtUeSector(s, inst));
            
            % calculando a diferença entre os ANGULOS VERTICAIS ACIMA
            mtDifAngsVer_Esp(s, jj) = abs(angElevUE_Esp - angVertUeAtivo); 
        end
    end
end    


% matriz de diferença entre o angulo AZIMUTAL de cada UEs para o angulo de STEERING da BS de cada setor (em, GRAUS º)
mtdifAngHor_Esp = zeros(S, numUE);  % [linha, coluna] = [setor, UE]

% laço percorrendo cada UE's p/ calcular os valores da matriz 'mtdifAngHor_Esp'
for jj = 1:numUE,
    
    % angulo azimutal do UE, em
    angUE180 = (180/pi)*angle(vtUePos(jj));       % ~ [-180º, 180º]
    angUE360 = wrapTo360(angUE180);               % ~ [0º, 360º] 
    
    % 'setor' é o nº do setor que UEjj pertence; 'inst' é o instante que o UEjj está ativo
    [setor, inst] = find(mtUeSector == jj);       
    
    % se temos algum UE INATIVO no slot de tempo 'inst', então pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % laço percorrendo cada SETOR
    for s = 1:S,
        
        % se o setor 's' é diferente do setor UEjj (laço externo), então
        if s ~= setor,
            % angulo horizontal do UE ativo no instante 'inst' no setor 'setor'
            angUeAtivo = wrapTo360((180/pi).*angle(vtUePos(mtUeSector(s, inst)))); % [0º, 360º]
            
            % calcular a diferença entre o angulos horizontais em
            dif1 = min(abs(angUE360-angUeAtivo));  % ~[0º, 360º]
            dif2 = min(abs(angUE180-wrapTo180(angUeAtivo))); % ~[-180º, +180º]
            
            mtdifAngHor_Esp(s, jj) = min(dif1, dif2);    
        end 
    end 
end

% TENSORES para calcular 
Av_Esp = []; % padrão RADIAÇÃO VERTICAL da antena para cada ANGULO \theta_3dB
Ah_Esp = []; % padrão RADIAÇÃO HORIZONTAL da antena para cada angulo \phi_3dB
A_ESP = [];  % padrão RADIAÇÃO TOTAL p/ Beamforming ESPECIFICO   
H_ESP = [];  % COEFICIENTE do CANAL ao quadrado em dB (sem fast-fading)

% laço percorrendo cada angulo \theta_3dB, \phi_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    % padrão de radiação VERTICAL
    Av_Esp(:,:,ii) = -min(12.*((mtDifAngsVer_Esp./vtTheta3dB_Esp(ii)).^2), SLA);  
    
    % padrão de radiação HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtdifAngHor_Esp./vtFi3dB_Esp(ii)).^2), Am);
    
    % Padrão de RADIAÇÃO TOTAL
    A_ESP(:,:,ii) = -min(-(Ah_Esp(:,:,ii) + Av_Esp(:,:,ii)), Am);

    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_ESP(:,:,ii) = G_BS + A_ESP(:,:,ii) - mtPL + mtNormal;   % [dB]
    
end

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR para Beamforming ESPECIFICO
Y_ESP = [];   % SINR

% laço percorrendo as dimensões do tensor h_esp
for jj = 1:length(vtFi3dB_Esp),
    
    % laço percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj está ativo
        [setor, inst] = find(mtUeSector == ii);      
        
        % soma o h_esp interferentes de cada UE
        aux = sum(h_esp(:, ii, jj)) - h_esp(setor, ii, jj);
        
        % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
        Y_ESP(jj, ii) = (Pot*h_esp(setor, ii, jj))/(Pot*aux + PN);
    
    end
end

% SINR em dB
YESP_dB = 10*log10(Y_ESP);


%% BEAMFORMING GRUPO-ESPECIFICO

% valores tirados do artigo 
vtFi3dB_gr = [30 20 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_gr = [10 10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]
B = 16;                            % numero de padrões de feixes (ou, GRUPOS)
Bh = 8;                            % numero de feixes horizontais p/ cada setor
Bv = B/Bh;                         % numero de feixes vertiais p/ cada setor

% mtDirFeixes = zeros(Bh, Bv);     % matriz de direção de feixes, onde cada elemento: (Fi_st, Theta_tild)
% mtDirFeixes = cell(Bh, Bv);      % celula de direção de feixes

angsFiSt_gr = zeros(1, Bh);
andDownTild_gr = zeros(1, Bv);

angHorInic = 0;
angHorFinal = 360;
passo = (angHorFinal - angHorInic)/(Bh*FatorSetor);

angsFiSt_gr = linspace(angHorInic + (passo/2), angHorFinal - (passo/2), Bh*FatorSetor);

% matriz de angulos \fi_st para uma celula
mtAngsFiSt_gr = reshape(angsFiSt_gr, Bh , FatorSetor);   
mtAngsFiSt_gr = mtAngsFiSt_gr';                  % [linha, coluna] = [setor, angulo Thi_st p/ cada grupo do setor]

% angDownTild_gr = linspace(-90, 90, Bv);
angDownTild_gr = [8 10];

% matriz de angulos \theta_thild p/ uma celula;
mtAngsThetaTild_gr = repmat(angDownTild_gr, FatorSetor*M, 1); % linha: setor, coluna: \theta_tild 's p/ cada setor

% matriz DIFERENÇA do ang. AZIMUTAL de cada UE para o ang. STEERING (\phi_st) da BS de cada setor
mtdifAngsHor_gr = zeros(S, numUE);  % [linha, coluna] = [Setor, UE] 

% laço percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsHor_gr'
for jj = 1:numUE,
    
    % 'setor' que o UEjj pertence e o slot_de_tempo 'inst' que o UEjj está ativo 
    [setor, inst] = find(mtUeSector == jj);    
    
    % se temos algum UE INATIVO no slot de tempo 'inst', então pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % calcula a diferença entre angulos
    angFiGrupo_360 = mtAngsFiSt_gr(setor,:);     % GRUPOS de angulos de steering [0º, 360] disponiveis para o setor_'setor' onde UEjj pertence
    angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % GRUPOS de angulos de steering [-180º, +180º] disponiveis para o setor_'setor' onde UEjj pertence
    
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  % angulo azimutal do UE atual em [-180, 180]
    anguloUE360 = wrapTo360(anguloUE180);        % angulo azimutal do UE atual em [0, 360]
   
    dif1 = min(abs(anguloUE360 - angFiGrupo_360));
    dif2 = min(abs(anguloUE180 - angFiGrupo_180));
    difAngulo = min(dif1, dif2);
    
    mtdifAngsHor_gr(setor, jj) = difAngulo;
    
    % laço percorrendo os SETORES 
    for s = 1:S, 
        % se o setor 's' é diferente do setor do UEjj (laço externo), então 
        if s ~= setor,
            angUeAtivo = wrapTo360((180/pi).*angle(vtUePos(mtUeSector(s, inst)))); % UE ativo no instante 'inst' no setor 'setor'
            [M, II] = min(abs(angUeAtivo - mtAngsFiSt_gr(s,:)));
            ddif1 = min(abs(anguloUE360 - mtAngsFiSt_gr(s,II)));
            ddif2 = min(abs(anguloUE180 - wrapTo180(mtAngsFiSt_gr(s,II))));
            mtdifAngsHor_gr(s, jj) = min(ddif1, ddif2);
        end
    end
end


% diferença de angulos de elevação de cada UE de grupo
mtdifAngsVer_gr = zeros(S, numUE);

% laço percorrendo todos UE's p/ calcular os valores da matriz 'mtdifAngsVer_gr'
for jj= 1:numUE,
    
    % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj está ativo
    [setor, inst] = find(mtUeSector == jj);       
    
    % se temos algum UE INATIVO no slot de tempo 'inst', então pulamos os calculos abaixo para UEjj
    if find(mtUeSector(:, inst) == 0), 
        continue;
    end
    
    % calculo o angulo de elevação THETA do UE
    angElevUE180 = mtThetaUE(setor, jj);
    
    mtdifAngsVer_gr(setor, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(setor,:)));
    
    % laço percorrendo todos os  SETORES
    for s = 1:S,
        % Percorrendo os OUTROS setores, i.e, setores diferente do setor que contém o UEjj (laço externo)
        if s ~= setor,
            angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));    % ang. de ELEVAÇÂO do UE ativo no setor 's' no instante 'inst'
            [M, II] = min(abs(angElevUeAtivo - mtAngsThetaTild_gr(s,:)));
            mtdifAngsVer_gr(s, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(s, II)));    
        end
    end    
end

% TENSORES para calcular o padrão de RADIAÇÃO 
Ah_gr = [];  % HORIZONTAL da antena para cada angulo \phi_3dB 
Av_gr = [];  % VERTICAL da antena para cada angulo \theta_3dB
A_GR = [];   % TOTAL p/ Beamforming GRUPO

% laço percorrendo cada angulo \theta_3dB, \phi_3dB p/ calcular cada dimensão do setor
for ii = 1:length(vtTheta3dB_Esp),

    % Ah para cada \fi_3dB
    Ah_gr(:,:,ii) = -min(12.*((mtdifAngsHor_gr./vtFi3dB_gr(ii)).^2), Am);
    
    % Av para cada \theta_3dB
    Av_gr(:,:,ii) = -min(12.*((mtdifAngsVer_gr./vtTheta3dB_gr(ii)).^2), SLA);
    
    % Padrão de radiação TOTAL
    A_GR(:,:,ii) = -min(-(Av_gr(:,:,ii) + Ah_gr(:,:,ii)), Am);
    
    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_GR(:,:,ii) = G_BS + A_GR(:,:,ii) - mtPL + mtNormal;   % [dB]
end


% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_gr = db2lin(H_GR);

% SINR para Beamforming GRUPO
Y_GR = [];   % SINR

% laço percorrendo cada dimensao do tensor p/ calcular cada dimensão do setor
for jj = 1:length(vtFi3dB_Esp),
  
    % laço percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj está ativo
        [setor, inst] = find(mtUeSector == ii);
        
        % soma h_gr interferentes para cada UE 
        aux = sum(h_gr(:, ii, jj)) - h_gr(setor, ii, jj);
        
        % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
        Y_GR(jj, ii) = (Pot*h_gr(setor, ii, jj))/(Pot*aux + PN);     % linha: SNIR p/ cada \theta_3dB e \fi_3dB
        
    end
end

% SINR em dB
YGR_dB = 10*log10(Y_GR);


%% PLOTANDO OS GRAFICOS


figure;
cdfplot(Y2D_dB)          % plotar SINR p/ 2DBF ou Beamforming Convencional
hold on;

% figure;
cdfplot(YESP_dB(1, :))   % plotar SINR p/ 3DBF UE-Especifico
cdfplot(YESP_dB(2, :))
cdfplot(YESP_dB(3, :))
%legend('Conventional', 'UE especifica - (\theta_{3dB}, \phi_{3dB) = (70º, 10º)}');

cdfplot(YGR_dB(1,:))     % plotar SINR p/ 3DBF Grupo-UE    
legend('2DBF', '3DBF-Usuários', '3DBF-Grupos (16 grupo)');
xlabel('SINR (dB)')
ylabel('CDF')
title('')
