%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

bandWidth = 10e+6;                                   % largura de banda 
fc = 3.5;                                            % frequencia de portadora [Ghz] na small-cell
R = 40;                                              % raio da pequena celula
numUE = 20;                                          % numero de UE's por macro-setores
xBS = 0;                                             % Posição do eixo x da BS
yBS = 0;                                             % Posição do eixo y da BS
vtSector = [ R*exp( 1j*[0 2*pi/3 4*pi/3] ) ];        % vetor marcação dos pontos de sectorização
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

% gerando os vetores para UE's (equipamento do usuário) aleatorios e iid.                       
angUE = 2*pi*rand(1,numUE);                            % angulo azimutal aleatorio dos UE's iid entre [0,2*PI] (RADIANOS)
raioUE = (R-10).*rand(1,numUE)+10;                     % raio aleatorio dos UE's iid entre (10,R)
vtUePos = raioUE.*exp(1j*angUE);                       % vetor de posições de UE's

% plotando as POSIÇÕES da BS e UE's    
hold on;
circle(xBS,yBS,R);    
plot(xBS,yBS,'*r','MarkerSize',16);                    % plotando a BS
plot(vtUePos,'*b','MarkerSize',12);                    % plotando os UE's, representados por sinal de + da cor verde
grid on;

% gerando as LINHAS dos SETORES. 
x = zeros(1,2*length(vtSector));
y = zeros(1,2*length(vtSector));
x([2:2:2*length(vtSector)]) = real(vtSector);
y([2:2:2*length(vtSector)]) = imag(vtSector);
plot(x,y,'y');
legend('Celula', 'BS', 'UEs','Setores')

% Encontrando os indices dos UE's de cada sector.
sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que estão no setor 1
sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que estão no setor 2
sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que estão no setor 3

% Testando o codigo acima dos setores: 
% resultado = vtUePos(sector1);          % captura as distancia dos UE's que estão dentro do setor 1
% plot(resultado,'bo');                  % plota o valores do 'resultados' com um circulo azul  

hold off


%% CALCULAR O PATH-LOSS DOS UE'S
vtDistUEtoBS = abs(vtUePos);                                     % Calculando as distancias entre UE's e BS:

% PL = 36.7*log10(sort(vtDistUEtoBS)) + 22.7 + 26*log10(fc);     % Path Loss em [db]          
PL = 36.7*log10((vtDistUEtoBS)) + 22.7 + 26*log10(fc);
figure;                                                          % gera uma nova figura para plotar os graficos
%plot(sort(vtDistUEtoBS), PL)
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (KM)')
ylabel('PL (DB)')
title('Path Loss ')


%%  BEAMFORMING 2D %%
%vtShadowing = sigma*randn(1, numUE);                          % Vetor de sombreamento para cada UE's [dB]
vtLogNormal = lognrnd(0,db2lin(4),1, numUE);                   % Vetor de shadowing Log Normal (?)
vtFastFad = (1/sqrt(2))*[randn(1,numUE) + 1j*randn(1,numUE)];  % vetor de Desvanecimento Rapido para cada UE's  
fi3dB_2D = 70;                                                 % largura de feixe de 3 dB na horizontal [º (GRAUS)]
theta3dB_2D = 10;                                              % largura de feixe de 3 dB na vertical   [º (GRAUS)]
angDownTild_2d = 8;                                            % angulo de down-tild (FIXO) [º (GRAUS)]
ang_st = [60, 180, 300];                                       % vetor de angulos de sterring para cada setor
angUEd = (180/pi).*angUE;                                      % angulo azimutal dos UE's em [º (GRAUS)]

AH_2d = zeros(1, numUE);                                       % vetor de padrão de radiação horizontal
    
% Padrão de radiação na HORIZONTAL 
AH_2d(sector1) = -min((12.*(((angUEd(sector1) - ang_st(1))/fi3dB_2D).^2)), Am);
AH_2d(sector2) = -min((12.*(((angUEd(sector2) - ang_st(2))/fi3dB_2D).^2)), Am);
AH_2d(sector3) = -min((12.*(((angUEd(sector3) - ang_st(3))/fi3dB_2D).^2)), Am);

% Padrão de radiação na VERTICAL
thetaUEs = atand((H_BS - H_UE)./(vtDistUEtoBS));                                   % [º em GRAUS]
AV_2d = -min(12.*(((thetaUEs-angDownTild_2d)./theta3dB_2D).^2), SLA);            

% Padrão de radiação TOTAL
A = -min(-(AH_2d + AV_2d), Am);

% Passando os parametros para escala linear
G_L = db2lin(G_BS);                                       % Ganho LINEAR da antena BS 
A_L = db2lin(A);                                          % Padrão da Antena Total em escala Linear 
PL_L = db2lin(PL);

h_2D = sqrt(numUE.*A_L.*PL_L.*vtLogNormal).*vtFastFad;    % coeficientes do canal para o 2D
Pot = dbm2lin(P_BS);                                      % Potência da antena transmissora
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);                 % Potência do Ruído

% coeficientes do canal interferentes
hI1_2D = sqrt(numUE.*A_L.*PL_L.*vtLogNormal).*vtFastFad;
hI2_2D = sqrt(numUE.*A_L.*PL_L.*vtLogNormal).*vtFastFad; 

SINR_2D = zeros(1, numUE);
% SINR_2D(sector1) = Pot*abs(h_2D(sector1))./(Pot.*() + PN);

