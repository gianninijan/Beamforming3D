%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setorização, i.e, setores/celulas
S = M*FatorSetor;                                    % número de setores. S = {1, 2, 3, ..., }      
numUE = 1000;                                           % numero de UE's por macro-setores
R = 80;                                              % raio da pequena celula
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

% angUE = [pi/3, pi, 5*pi/3];             % angulos de teste
% raioUE = (R/2).*ones(1,numUE);          % raios de teste
% vtUePos = raioUE.*exp(1j*angUE);        % vetor de posição das UE's de teste

vtUePos = [];                             % vetor de posição de usuários

% gerando os UE's uniformemente entre [10, R]
for ii = 1:numUE,
    while true,
         vtUePos(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
         
         if (abs(vtUePos(ii)) < R) && (abs(vtUePos(ii)) > 10),
             break;
         end
    end
end

% vetor das posições da BS de cada setor
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de sterring p/ cada celula [Radianos]
raio_bs = 0;                                        % raio das BS's - gerar aleatorio
ang_bs = 0;                                         % angulos de posições das BS's - gerar aleatorio
posBS = raio_bs.*exp(1j*ang_bs);                    % forma exponencial das BS's
vtBsSetor = repelem(posBS, FatorSetor);             % vetor de posição da BS de cada setor
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% plotando as POSIÇÕES da BS e UE's    
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
angUE = wrapTo2Pi(angle(vtUePos));
sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que estão no setor 1
sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que estão no setor 2
sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que estão no setor 3

% Testando o codigo acima dos setores: 
resultado = vtUePos(sector2);          % captura as distancia dos UE's que estão dentro do setor 1
plot(resultado,'ro');                  % plota o valores do 'resultados' com um circulo azul  

hold off


%% CALCULAR O PATH-LOSS DOS UE'S

% matriz de distancia de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                             

% laço percorrendo cada setor
for ii = 1:S
    
    % laço percorrendo todos os usuários
    for jj = 1:numUE
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % elementos da matriz de distância
    end
    
end

% P/ calcular o setor de cada UE devemos levar em conta a distancia p/ cada BS do setor e angulo

mtPL = PathLoss(mtDist, fc);                  % matriz de PATH LOSS de cada setor para cada usuário
PL = min(mtPL);                               % valor minimo
vtDistUEtoBS = min(mtDist);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')


%%  BEAMFORMING 2D %%
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)]; % vetor de Desvanecimento Rapido para cada UE's
% load('desvanecimento.mat');

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS]

% matriz de shadowing normal
mtNormal_2D = sigma.*randn(S,numUE);    % Linhas: Setores; Colunas: Usuários

% Angulo \theta de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [GRAUS]

% calculando o padrão vertical da antena
Av_2D = -min(12.*(((mtThetaUE - angDownTild_2d)./theta3dB_2D).^2), SLA);    % linhas: setores, colunas: UE's

% calculando o padrão horizontal da antena
Ah_2D = [];

% posições dos UE's 
posicoes_d = (180/pi).*angle(vtUePos); % [-180, 180]

% diferença do angulo azimutal de cada UE p/ cada angulo de steering de cada setor
mtdifAngulos_d = [];

% laço percorrendo cada UE
for jj = 1:numUE,
    
    % angulo de cada usuario em [-180, 180]
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  
   
    % angulo de cada usuario em [0, 360]
    anguloUE360 = wrapTo360(anguloUE180);
    
    % calculo o angulo de steering para cada setor entre [0, 360]
    angStrTo360 = (180/pi).*vtAngST;
    
    % calculo o angulo de steering para cada setor entre [-180, 180]
    angStrTo180 = wrapTo180(angStrTo360);

    % laço percorrendo cada setor
    for ii = 1:S,
        
        % calcula a diferença entre angulos
        dif1 = abs(anguloUE360 - angStrTo360(ii));
        dif2 = abs(anguloUE180 - angStrTo180(ii));
        difAngulo = min(dif1, dif2);
        
        mtdifAngulos_d(ii, jj) = difAngulo;
    end
end

% padrão de radiação HORIZONTAL
Ah_2D = -min(12.*((mtdifAngulos_d./fi3dB_2D).^2), Am);

% Padrão de radiação TOTAL
A = -min(-(Ah_2D + Av_2D), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_2d = G_BS + A + mtPL + mtNormal_2D;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_2d = db2lin(H_2d);

% Potência do Ruído Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Potência da antena transmissora
Pot = dbm2lin(P_BS);                       

% SINR
Y = [];   % SINR

% capturando o setor de cada UE atraves do INDICES
[X,I] = min(sqrt(mtdifAngulos_d.^2), [], 1);

% laço percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
     
    aux = sum(h_2d(:, ii)) - h_2d(I(ii), ii);
    
    % calculo da SNIR
    Y(ii) = (Pot*h_2d(I(ii), ii))/(Pot*aux + PN);
    
end

% SINR em dB
Y_dB = 10*log10(Y);

figure;
cdfplot(Y_dB)

% angST = [60 180 300];
% angUE = [30 60 90 120 150 180 210 240 270 300 330 360];
% angST180 = wrapTo180(angST);
% angUE180 = wrapTo180(angUE);
% mt = [];
% mt2 = [];
% 
% for ii = 1:length(angUE),
%     mt(:, ii) = angUE(ii) - angST';
% end
% 
% for ii = 1:length(angUE),
%     mt2(:, ii) = angUE180(ii) - angST180';
% end