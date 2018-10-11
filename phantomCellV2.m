%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setorização, i.e, setores/celulas
S = M*FatorSetor;                                    % número de setores. S = {1, 2, 3, ..., }      
numUE = 20;                                         % numero de UE's por macro-setores
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

%% Dados Comuns para os Beamforming

% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)]; % vetor de Desvanecimento Rapido para cada UE's
% load('desvanecimento.mat');

% matriz de shadowing normal
mtNormal = sigma.*randn(S,numUE);    % Linhas: Setores; Colunas: Usuários

% Angulo \theta de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [GRAUS]

% Potência do Ruído Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Potência da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS]

% calculando o padrão vertical da antena
Av_2D = -min(12.*(((mtThetaUE - angDownTild_2d)./theta3dB_2D).^2), SLA);    % linhas: setores, colunas: UE's

% calculando o padrão horizontal da antena
Ah_2D = [];

% posições dos UE's 
% posicoes_d = (180/pi).*angle(vtUePos); % [-180, 180]

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
A_2D = -min(-(Ah_2D + Av_2D), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_2d = G_BS + A_2D + mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_2d = db2lin(H_2d);

% SINR
Y2D = [];   % SINR

% capturando o setor de cada UE atraves do INDICES
[X,I] = min(sqrt(mtdifAngulos_d.^2), [], 1);

% laço percorrendo cada UE para calcular a SINR
for ii = 1:numUE,
     
    aux = sum(h_2d(:, ii)) - h_2d(I(ii), ii);
    
    % calculo da SNIR
    Y2D(ii) = (Pot*h_2d(I(ii), ii))/(Pot*aux + PN);
    
end

% SINR em dB
Y2D_dB = 10*log10(Y2D);

figure;
cdfplot(Y2D_dB)
hold on;

%% BEAMFORMING UE ESPECIFICO

% valores tirados do artigo
vtFi3dB_Esp = [30 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_Esp = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]

% angulo de donwtild p/ Beamforming UE especifico será igual a \theta de cada usuário
angDownTild_Esp = atand((H_BS - H_UE)./(mtDist));

% tensor para calcular o padrão de radiação vertical da antena para cada angulo \theta_3dB
Av_Esp = [];

% laço percorrendo cada angulo theta_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    Av_Esp(:,:,ii) = -min(12.*(((mtThetaUE - angDownTild_Esp)./vtTheta3dB_Esp(ii)).^2), SLA);  

end

% tensor para calcular o padrão de radiação horizontal da antena para cada angulo \theta_3dB
Ah_Esp = [];

% diferença do angulo azimutal de cada UE p/ cada angulo de steering de cada setor
mtdifAnguloEsp_d = mtdifAngulos_d;

% laço percorrendo cada UE's
for jj = 1:numUE,
    
    % mudando apenas a diferça de angulos na horizontal p/ o setor de cada UE
    mtdifAnguloEsp_d( I(jj), jj) = 0;
    
end


% laço percorrendo cada angulo fi_3dB
for ii = 1:length(vtFi3dB_Esp),

    % padrão de radiação HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtdifAnguloEsp_d./vtFi3dB_Esp(ii)).^2), Am);
    
end

% Padrão de radiação TOTAL
A_ESP = -min(-(Ah_Esp + Av_Esp), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_ESP = G_BS + A_ESP + mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR
Y_ESP = [];   % SINR

% laço percorrendo cada dimensao do tensor
for jj = 1:length(vtFi3dB_Esp),
    
    % laço percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        aux = sum(h_esp(:, ii, jj)) - h_esp(I(ii), ii, jj);
        
        % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
        Y_ESP(jj, ii) = (Pot*h_esp(I(ii), ii, jj))/(Pot*aux + PN);
        
    end
    
end

% SINR em dB
YESP_dB = 10*log10(Y_ESP);

hold off;
figure;
cdfplot(YESP_dB(1, :))
cdfplot(YESP_dB(2, :))
cdfplot(YESP_dB(3, :))
