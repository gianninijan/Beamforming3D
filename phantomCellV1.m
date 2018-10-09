%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setorização, i.e, setores/celulas
S = M*FatorSetor;                                    % número de setores. S = {1, 2, 3, ..., }      
numUE = 20;                                          % numero de UE's por macro-setores
R = 80;                                              % raio da pequena celula
R_BS = 1e-6;                                         % raio da posição BS
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

% vetor das posições dos UE's (equipamento do usuário) aleatorios e iid.                       
angUE = 2*pi*rand(1,numUE);                         % angulo azimutal aleatorio dos UE's iid entre [0,2*PI] (RADIANOS)
raioUE = (R-10).*rand(1,numUE)+10;                  % raio aleatorio dos UE's iid entre (10,R)
vtUePos = raioUE.*exp(1j*angUE);                    % vetor de posições de UE's

% vetor das posições da BS de cada setor
ang_st = [pi/3, pi, 5*pi/3];                              % angulos de sterring p/ cada celula
raio_bs = 0;                                              % raio das BS's
ang_bs = 0;                                               % angulos de posições das BS's
posBS = raio_bs*exp(1j*ang_bs);                           % forma exponencial das BS's
mtBSPos = [ repelem(posBS, S)' repmat(ang_st(:), M, 1)];  % cada linha [posBS_setor ang_stl]

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

mtDist = zeros(S,numUE);                         % matriz de distancia de cada UE para BS de cada setor               

% laço percorrendo cada setor
for ii = 1:S
    
    % laço percorrendo todos os usuários
    for jj = 1:numUE
        mtDist(ii,jj) = norm(vtUePos(jj) - mtBSPos(ii,1));  % elementos da matriz de distância
    end
end

vtDistUEtoBS = abs(vtUePos);                     % Calculando as distancias entre UE's e BS:
PL = PathLoss(vtDistUEtoBS, fc);                 % Path Loss em [db]
figure;                                          % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')

mtPL = PathLoss(mtDist, fc);                    % matriz de PATH LOSS de cada setor para cada usuário
% PL = min(mtPL);                               % valor minimo


%%  BEAMFORMING 2D %%

%vtShadowing = sigma*randn(1, numUE);                          % Vetor de sombreamento para cada UE's [dB]
vtLogNormal = lognrnd(0,db2lin(4),1, numUE);                   % Vetor de shadowing Log Normal (?)
vtFastFad = (1/sqrt(2))*[randn(1,numUE) + 1j*randn(1,numUE)];  % vetor de Desvanecimento Rapido para cada UE's  

mtLogNormal = lognrnd(0,db2lin(4),S, numUE);
mtFastFad = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];

fi3dB_2D = 70;                                                 % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                                              % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                                            % angulo de down-tild (FIXO) [GRAUS]
angUEd = (180/pi).*angUE;                                      % angulo azimutal dos UE's em [GRAUS]

% Padrão de radiação na HORIZONTAL 
AH_2d = zeros(1, numUE);                                       % vetor de padrão de radiação horizontal
AH_2d(sector1) = padrao_Horizontal(angUE(sector1), ang_st(1), fi3dB_2D, Am);
AH_2d(sector2) = padrao_Horizontal(angUE(sector2), ang_st(2), fi3dB_2D, Am);
AH_2d(sector3) = padrao_Horizontal(angUE(sector3), ang_st(3), fi3dB_2D, Am);

% Padrão de radiação na VERTICAL
thetaUEs = atand((H_BS - H_UE)./(vtDistUEtoBS));                           % [GRAUS]
AV_2d = padrao_Vertical(thetaUEs, angDownTild_2d, theta3dB_2D, SLA);       % [dB]
           
% Padrão de radiação TOTAL
A = -min(-(AH_2d + AV_2d), Am);

% Passando os parametros para escala linear
G_L = db2lin(G_BS);                                       % Ganho LINEAR da antena BS 
A_L = db2lin(A);                                          % Padrão da Antena Total em escala Linear 
PL_L = db2lin(PL);

h_2D = canal(G_L, A_L, PL_L, vtLogNormal, vtFastFad);     % coeficientes do canal para o 2D
     
Pot = dbm2lin(P_BS);                                      % Potência da antena transmissora
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);                 % Potência do Ruído

% coeficientes do canal interferentes
hi21_2D = sqrt(G_L.*A_L.*PL_L.*vtLogNormal).*vtFastFad;   % sinal interferente do setor 2 p/ os usuarios do setor 1
hi31_2D = sqrt(G_L.*A_L.*PL_L.*vtLogNormal).*vtFastFad;   % sinal interferente do setor 3 p/ os usuarios do setor 1

SINR_2D = zeros(1, numUE);
% SINR_2D(sector1) = Pot*abs(h_2D(sector1))./(Pot.*() + PN);


% angulos_teste = [pi/3 3*pi/4 5*pi/4 7*pi/4 ];      % = [60º 135º 225º 315º]
% lambdaWrapped = (180/pi).*wrapToPi(angulos_teste); % = [60º 135º -135º -45º] % Se angulo<0: angulo + 360 
