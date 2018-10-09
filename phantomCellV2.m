%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setoriza��o, i.e, setores/celulas
S = M*FatorSetor;                                    % n�mero de setores. S = {1, 2, 3, ..., }      
numUE = 20;                                          % numero de UE's por macro-setores
R = 80;                                              % raio da pequena celula
xBS = 0;                                             % Posi��o do eixo x da BS
yBS = 0;                                             % Posi��o do eixo y da BS
vtSector = [ R*exp( 1j*[0 2*pi/3 4*pi/3] ) ];        % vetor marca��o dos pontos de sectoriza��o
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

%% POSI��ES DA BS E UE'S %%

% vetor das posi��es dos UE's (equipamento do usu�rio) aleatorios e iid.                       
% angUE = 2*pi*rand(1,numUE);             % angulo azimutal aleatorio dos UE's iid entre [0,2*PI] (RADIANOS)             
% raioUE = (R-10).*rand(1,numUE)+10;      % raio aleatorio dos UE's iid entre (10,R)
% vtUePos = raioUE.*exp(1j*angUE);        % vetor de posi��es de UE's
% angUE = [pi/3, pi, 5*pi/3];             % angulos de teste
% raioUE = (R/2).*ones(1,numUE);          % raios de teste
% vtUePos = raioUE.*exp(1j*angUE);        % vetor de posi��o das UE's de teste

vtUePos = [];                             % vetor de posi��o de usu�rios

% gerando os UE's uniformemente entre [10, R]
for ii = 1:numUE,
    while true,
         vtUePos(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
         
         if (abs(vtUePos(ii)) < R) && (abs(vtUePos(ii)) > 10),
             break;
         end
    end
end

% vetor das posi��es da BS de cada setor
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de sterring p/ cada celula [Radianos]
raio_bs = 0;                                        % raio das BS's - gerar aleatorio
ang_bs = 0;                                         % angulos de posi��es das BS's - gerar aleatorio
posBS = raio_bs.*exp(1j*ang_bs);                    % forma exponencial das BS's
vtBsSetor = repelem(posBS, FatorSetor);             % vetor de posi��o da BS de cada setor
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% plotando as POSI��ES da BS e UE's    
hold on;
circle(xBS,yBS,R);    
plot(xBS,yBS,'*r','MarkerSize',16);       % plotando a BS
plot(vtUePos,'*b','MarkerSize',12);       % plotando os UE's, representados por sinal de + da cor verde
grid on;

% gerando as LINHAS dos SETORES. 
x = zeros(1,2*length(vtSector));
y = zeros(1,2*length(vtSector));
x([2:2:2*length(vtSector)]) = real(vtSector);
y([2:2:2*length(vtSector)]) = imag(vtSector);
plot(x,y,'y');
legend('Celula', 'BS', 'UEs','Setores')

% Encontrando os indices dos UE's de cada sector.
% sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que est�o no setor 1
% sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que est�o no setor 2
% sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que est�o no setor 3

% Testando o codigo acima dos setores: 
% resultado = vtUePos(sector1);          % captura as distancia dos UE's que est�o dentro do setor 1
% plot(resultado,'bo');                  % plota o valores do 'resultados' com um circulo azul  

hold off


%% CALCULAR O PATH-LOSS DOS UE'S

mtDist = zeros(S,numUE);              % matriz de distancia de cada UE (coluna) para BS (linha) de cada setor               

% la�o percorrendo cada setor
for ii = 1:S
    
    % la�o percorrendo todos os usu�rios
    for jj = 1:numUE
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % elementos da matriz de dist�ncia
    end
    
end

% para calcular o setor de cada UE devemos levar em conta a distancia para 
% cada BS do setor e angulo
mtPL = PathLoss(mtDist, fc);                  % matriz de PATH LOSS de cada setor para cada usu�rio
PL = min(mtPL);                               % valor minimo
vtDistUEtoBS = min(mtDist);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')


%%  BEAMFORMING 2D %%
% mtLogNormal_2D = lognrnd(0,db2lin(4),S, numUE);                    % matriz de shadowing Log Normal (?)
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)]; % vetor de Desvanecimento Rapido para cada UE's
% load('desvanecimento.mat');
mtNormal_2D = sigma.*randn(S,numUE);    % Linhas: Setores; Colunas: Usu�rios

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS]

% Angulo \theta de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [GRAUS]

% calculando o padr�o vertical da antena
Av_2D = -min(12.*(((mtThetaUE - angDownTild_2d)./theta3dB_2D).^2), SLA);    % linhas: setores, colunas: UE's

% calculando o padr�o horizontal da antena
Ah_2D = [];

% la�o percorrendo cada setor (linhas de Ah_2D) 
for ii = 1:S,
    
    % la�o percorrendo cada UE (colunas de Ah_2D)
    for jj = 1:numUE,
        
    end
    
end

% calculando a diferen�a entre angulos para as formulas da A_V e A_H.
% mtDifFhis_2D = diffAngles(vtAngST, angUE);    % linha: cada UE, coluna: phiST_setor [GRAUS]

% calculando os padr�es verticais e horizontais
% Ah_2D = -min(12.*((mtDifFhis_2D./fi3dB_2D).^2), Am);


% Padr�o de radia��o TOTAL
A = -min(-(Ah_2D + Av_2D), Am);

% Passando os parametros para escala linear
% G_L = db2lin(G_BS);                        % Ganho LINEAR da antena BS 
% mtA_L = db2lin(A);                         % matriz Padr�o da Antena Total em escala Linear 
% mtPL_L = db2lin(mtPL);                     % matriz de path loss de cada UE para cada setor
% Pot = dbm2lin(P_BS);                       % Pot�ncia da antena transmissora
% PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  % Pot�ncia do Ru�do

% calculando os coeficientes do canal 
h_2D = canal(G_L, mtA_L, mtPL_L, mtLogNormal_2D, mtFastFad_2D);

% capturando o setor de cada UE atraves do INDICES
[X,I] = min(mtDifFhis_2D, [], 1);

% modulo ao quadrado da matriz h_2D
modQuadH_2d = abs(h_2D).^2;

Y = [];   % SINR

% percorrendo cada coluna de 
for ii = 1:size(mtDifFhis_2D,2)
    
    % somando todos os elementos do setor, menos o 
    aux = sum(modQuadH_2d(:, ii)) - modQuadH_2d(I(ii), ii);
    
    % calculo da SNIR
    Y(ii) = (Pot*modQuadH_2d(I(ii), ii))/(Pot*aux + PN);
end 

figure;
% cdfplot(Y)

% quantSector1 = length(find(I==1));
% quantSector2 = length(find(I==2));
% quantSector3 = length(find(I==3));

% angulos azimutais dos UE's em [GRAUS]
% angUEd = (180/pi).*angUE;          



