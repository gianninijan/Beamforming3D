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

% VALORES DE TESTE P/ numUE = 8;
% angUE = [0*pi, 5*pi/6, 11*pi/6, pi/3, pi/2, pi, 7*pi/6, 35*pi/18];    % angulos de teste
% raioUE = [R/2, R/4, R/8, R/4, R/8, R/2, R/8, R/4];                    % raios de teste
% vtUePos = raioUE.*exp(1j*angUE);        % vetor de posição das UE's de teste


% GERANDO AS POSIÇÕES DO UE's ~ U(10, R)

vtUePos = [];                             % vetor de posição de usuários aleatórios 

for ii = 1:numUE,
    while true,
         vtUePos(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
         
         if (abs(vtUePos(ii)) < R) && (abs(vtUePos(ii)) > 10),
             break;
         end
    end
end


% vetor das posições da BS de cada setor
raio_bs = 0;                                        % raio das BS's - gerar aleatorio
ang_bs = 0;                                         % angulos de posições das BS's - gerar aleatorio
posBS = raio_bs.*exp(1j*ang_bs);                    % forma exponencial das BS's
vtBsSetor = repelem(posBS, FatorSetor);             % vetor posição da BS's p/ cada setor


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
legend('Celula', 'BS', 'UEs','Setores')

% Encontrando os indices dos UE's de cada sector.
angUE = wrapTo2Pi(angle(vtUePos));                            % angUE ~ [0, 2*pi] radianos
sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que estão no setor 1
sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que estão no setor 2
sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que estão no setor 3

% criando a matriz de setores (linhas) e UE's de cada setor 
maximo = max([length(sector1), length(sector2), length(sector3)]);
mtUeSector = max(S, maximo);
mtUeSector(1, 1:length(sector1)) = sector1;
mtUeSector(2, 1:length(sector2)) = sector2;
mtUeSector(3, 1:length(sector3)) = sector3;

% Testando o codigo acima dos setores: 
% resultado = vtUePos(sector2);          % captura as distancia dos UE's que estão dentro do setor 1
% plot(resultado,'ro');                  % plota o valores do 'resultados' com um circulo azul  
% 
% hold off


%% CALCULAR O PATH-LOSS DOS UE'S

% matriz de DISTÂNCIA de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                             

% laço percorrendo cada setor (coluna)
for ii = 1:S
    
    % laço percorrendo todos os usuários (linha)
    for jj = 1:numUE
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % DIM( mtDist ) = [Linhas, Colunas] = [Setores, UE's]
    end
    
end

% P/ calcular o setor de cada UE devemos levar em conta a distancia p/ cada BS do setor e angulo
mtPL = PathLoss(mtDist, fc);                  % matriz de PATH LOSS de cada setor para cada usuário
PL = min(mtPL);                               % valor minimo

vtDistUEtoBS = min(mtDist);
%[vtDistUEtoBS, POS] = min(mtDist, [], 1);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')


%% DADOS COMUNS PARA OS BEAMFORMINGS

% vetor de DESVANECIMENTO RÁPIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];

% matriz de shadowing normal em dB
mtNormal = sigma.*randn(S,numUE);    % dim( mtNormal ) = [Linhas, Colunas] = [Setores, Usuários]

% Angulo \theta de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [GRAUS]

% Potência do Ruído Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Potência da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% angulos de BORESIGHT FIXO p/ a antena de cada setor, i.e, angulo azimutal na qual teremos o ganho máximo da antena 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de boresight p/ cada celula [Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS]

% calculando o padrão vertical da antena
Av_2D = -min(12.*(((mtThetaUE - angDownTild_2d)./theta3dB_2D).^2), SLA);    % [linhas, colunas] = [setores, UE's]

% calculando o padrão horizontal da antena
Ah_2D = [];

% posições dos UE's 
% posicoes_d = (180/pi).*angle(vtUePos); % [-180, 180]

% diferença do angulo azimutal de cada UE p/ cada angulo de steering de cada setor
mtdifAngulos_d = [];

% laço percorrendo cada UE
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
[X,I] = min(sqrt(mtdifAngulos_d.^2), [], 1); % X = menores elemento de cada coluna; I = posição do menor elemento na coluna, i.e, localização na linha do menor elemento de cada coluna

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

%% BEAMFORMING UE ESPECIFICO %%

% valores tirados do artigo
vtFi3dB_Esp = [70 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_Esp = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]

% angulo de donwtild p/ Beamforming UE especifico será igual a \theta de cada usuário
angDownTild_Esp = [];

% Calculando o angulo de downtild 
% laço percorrendo os usuários  
for jj = 1:numUE,
    
    % se UEjj pertence ao setor 1, então: 
    if find(sector1 == jj),
        
        [i, inst] = find(sector1 == jj);                     % inst (instante) é a posição do elemento jj no vetor sector1
        
        angDownTild_Esp(1,jj) = mtThetaUE(1, jj);            % atand((H_BS - H_UE)./(norm(vtUePos(jj) - vtBsSetor(1))));
        
        % se no instante inst temos UE ativo no setor2, então
        if inst <= length(sector2), 
            angDownTild_Esp(2,jj) = mtThetaUE(2, sector2(inst));   
        else
            angDownTild_Esp(2,jj) = angDownTild_2d;
        end
        
        % se no instante inst temos UE ativo no setor3, então
        if inst <= length(sector3),
            angDownTild_Esp(3,jj) = mtThetaUE(3, sector3(inst)); 
        else
            angDownTild_Esp(3,jj) = angDownTild_2d;
        end
    
    % se UEjj pertence ao setor 2, então:     
    elseif find(sector2 == jj),
        [i, inst] = find(sector2 == jj);            % inst (instante) é a posição do elemento jj no vetor sector2
        
        angDownTild_Esp(2,jj) = mtThetaUE(2, jj);
        
        % se no instante inst temos UE ativo no setor1, então
        if inst <= length(sector1),
            angDownTild_Esp(1,jj) = mtThetaUE(1, sector1(inst)); 
        else
            angDownTild_Esp(1,jj) = angDownTild_2d;
        end
        
        % se no instante inst temos UE ativo no setor3, então
        if inst <= length(sector3),
            angDownTild_Esp(3,jj) = mtThetaUE(3, sector3(inst)); 
        else
            angDownTild_Esp(3,jj) = angDownTild_2d;
        end

    % se UEj pertence ao setor 3, então:     
    elseif find(sector3 == jj),
        [i, inst] = find(sector3 == jj);  % inst (instante) é a posição do elemento jj no vetor sector3
        
        angDownTild_Esp(3,jj) = mtThetaUE(3, jj);  
        
        % se no instante inst temos UE ativo no setor1, então
        if inst <= length(sector1), 
            angDownTild_Esp(1,jj) = mtThetaUE(1, sector1(inst));   
        else
           angDownTild_Esp(1,jj) = angDownTild_2d;      
        end
        
        % se no instante inst temos UE ativo no setor2, então
        if inst <= length(sector2), 
            angDownTild_Esp(2,jj) = mtThetaUE(2, sector2(inst));   
        else
            angDownTild_Esp(2,jj) = angDownTild_2d;
        end
        
   end
    
end

% matriz que calcula a diferença entre os angulos no numerador Ah (padrão de radiação horizontal)
mtdifAnguloEsp_d = zeros(S, numUE);

% calculando a difença entre o angulo_fi de cada usuario com o angulo_fi dos setores 
% laço percorrendo cada UE's
for jj = 1:numUE,
    
    angUE180 = (180/pi)*angle(vtUePos(jj));       % ~ [-180º, 180º]
    angUE360 = wrapTo360(angUE180);               % ~ [0º, 360º]
    
    % se UEjj pertence ao setor 1, então: 
    if find(sector1 == jj),
        
        [i, inst] = find(sector1 == jj);                     % inst (instante) é a posição do elemento jj no vetor sector1
        
        % se no instante inst temos UE ativo no setor2, então
        if inst <= length(sector2),
            angUE180_S2 = (180/pi)*angle(vtUePos(sector2(inst)));          % ~ [-180º, 180º]
            angUE360_S2 = wrapTo360(angUE180_S2);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S2);
            dif2 = abs(angUE180 - angUE180_S2);
            mtdifAnguloEsp_d(2, jj) = min(dif1, dif2);
        else
           angUE360_S2 = (180/pi)*vtAngST(2);
           angUE180_S2 = wrapTo180(angUE360_S2);
           dif1 = abs(angUE360 - angUE360_S2);
           dif2 = abs(angUE180 - angUE180_S2);
           mtdifAnguloEsp_d(2, jj) = min(dif1, dif2); 
        end
        
        % se no instante inst temos UE ativo no setor3, então
        if inst <= length(sector3),
            angUE180_S3 = (180/pi)*angle(vtUePos(sector3(inst)));          % ~ [-180º, 180º]
            angUE360_S3 = wrapTo360(angUE180_S3);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S3);
            dif2 = abs(angUE180 - angUE180_S3);
            mtdifAnguloEsp_d(3, jj) = min(dif1, dif2);
        else
           angUE360_S3 = (180/pi)*vtAngST(3);                              % ~ [0º, 360º]
           angUE180_S3 = wrapTo180(angUE360_S3);                           % ~ [-180º, 180º]
           dif1 = abs(angUE360 - angUE360_S3);
           dif2 = abs(angUE180 - angUE180_S3);
           mtdifAnguloEsp_d(3, jj) = min(dif1, dif2); 
        end
    
    % se UEjj pertence ao setor 2, então: 
    elseif find(sector2 == jj),
        
        [i, inst] = find(sector2 == jj);                     % inst (instante) é a posição do elemento jj no vetor sector2
        
        % se no instante inst temos UE ativo no setor1, então
        if inst <= length(sector1),
            angUE180_S1 = (180/pi)*angle(vtUePos(sector1(inst)));          % ~ [-180º, 180º]
            angUE360_S1 = wrapTo360(angUE180_S1);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S1);
            dif2 = abs(angUE180 - angUE180_S1);
            mtdifAnguloEsp_d(1, jj) = min(dif1, dif2);
        else
           angUE360_S1 = (180/pi)*vtAngST(1);
           angUE180_S1 = wrapTo180(angUE360_S1);
           dif1 = abs(angUE360 - angUE360_S1);
           dif2 = abs(angUE180 - angUE180_S1);
           mtdifAnguloEsp_d(1, jj) = min(dif1, dif2); 
        end
        
        % se no instante inst temos UE ativo no setor3, então
        if inst <= length(sector3),
            angUE180_S3 = (180/pi)*angle(vtUePos(sector3(inst)));          % ~ [-180º, 180º]
            angUE360_S3 = wrapTo360(angUE180_S3);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S3);
            dif2 = abs(angUE180 - angUE180_S3);
            mtdifAnguloEsp_d(3, jj) = min(dif1, dif2);
        else
           angUE360_S3 = (180/pi)*vtAngST(3);
           angUE180_S3 = wrapTo180(angUE360_S3);
           dif1 = abs(angUE360 - angUE360_S3);
           dif2 = abs(angUE180 - angUE180_S3);
           mtdifAnguloEsp_d(3, jj) = min(dif1, dif2); 
        end
    
    % se UEjj pertence ao setor 3, então:
    elseif find(sector3 == jj),
        
        [i, inst] = find(sector3 == jj);                     % inst (instante) é a posição do elemento jj no vetor sector3
        
        % se no instante inst temos UE ativo no setor1, então
        if inst <= length(sector1),
            angUE180_S1 = (180/pi)*angle(vtUePos(sector1(inst)));          % ~ [-180º, 180º]
            angUE360_S1 = wrapTo360(angUE180_S1);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S1);
            dif2 = abs(angUE180 - angUE180_S1);
            mtdifAnguloEsp_d(1, jj) = min(dif1, dif2);
        else
           angUE360_S1 = (180/pi)*vtAngST(1);
           angUE180_S1 = wrapTo180(angUE360_S1);
           dif1 = abs(angUE360 - angUE360_S1);
           dif2 = abs(angUE180 - angUE180_S1);
           mtdifAnguloEsp_d(1, jj) = min(dif1, dif2); 
        end
        
        % se no instante inst temos UE ativo no setor2, então
        if inst <= length(sector2),
            angUE180_S2 = (180/pi)*angle(vtUePos(sector2(inst)));          % ~ [-180º, 180º]
            angUE360_S2 = wrapTo360(angUE180_S2);                          % ~ [0º, 360º]
            dif1 = abs(angUE360 - angUE360_S2);
            dif2 = abs(angUE180 - angUE180_S2);
            mtdifAnguloEsp_d(2, jj) = min(dif1, dif2);
        else
           angUE360_S2 = (180/pi)*vtAngST(2);
           angUE180_S2 = wrapTo360(angUE360_S2);
           dif1 = abs(angUE360 - angUE360_S2);
           dif2 = abs(angUE180 - angUE180_S2);
           mtdifAnguloEsp_d(2, jj) = min(dif1, dif2); 
        end
    end  
end

% TENSORES para calcular o padrão de 
Av_Esp = []; % RADIAÇÃO VERTICAL da antena para cada ANGULO \theta_3dB
Ah_Esp = []; % RADIAÇÃO HORIZONTAL da antena para cada angulo \phi_3dB
A_ESP = [];  % RADIAÇÃO TOTAL p/ Beamforming ESPECIFICO   

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_ESP = [];

% laço percorrendo cada angulo \theta_3dB, \phi_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    % padrão de radiação VERTICAL
    Av_Esp(:,:,ii) = -min(12.*(((mtThetaUE - angDownTild_Esp)./vtTheta3dB_Esp(ii)).^2), SLA);  
    
    % padrão de radiação HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtdifAnguloEsp_d./vtFi3dB_Esp(ii)).^2), Am);
    
    % Padrão de RADIAÇÃO TOTAL
    A_ESP(:,:,ii) = -min(-(Ah_Esp(:,:,ii) + Av_Esp(:,:,ii)), Am);

    % coeficientes do canal ao quadrado em dB (sem fast-fading)
    H_ESP(:,:,ii) = G_BS + A_ESP(:,:,ii) + mtPL + mtNormal;   % [dB]
end

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR para Beamforming ESPECIFICO
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

hold on;
% figure;
cdfplot(YESP_dB(1, :))
%cdfplot(YESP_dB(2, :))
%cdfplot(YESP_dB(3, :))
%legend('Conventional', 'UE especifica - (\theta_{3dB}, \phi_{3dB) = (70º, 10º)}');


%% BEAMFORMING GRUPO-ESPECIFICO

% valores tirados do artigo 
vtFi3dB_gr = [70 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_gr = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]
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
mtAngsFiSt_gr = reshape(angsFiSt_gr, Bh , FatorSetor*M);   
mtAngsFiSt_gr = mtAngsFiSt_gr';                  % [linha, coluna] = [setor, angulo Thi_st p/ cada grupo do setor]

% angDownTild_gr = linspace(-90, 90, Bv);
angDownTild_gr = [4 8];

% matriz de angulos \theta_thild p/ uma celula;
mtAngsThetaTild_gr = repmat(angDownTild_gr, FatorSetor*M, 1); % linha: setor, coluna: \theta_tild 's p/ cada setor

% diferença do angulo azimutal de cada UE de grupo
mtdifAngsHor_gr = zeros(S, numUE);

% laço percorrendo todos UE's
for jj = 1:numUE,
    
    [setor, inst] = find(mtUeSector == jj);      % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj está ativo 
    
    % calcula a diferença entre angulos
    angFiGrupo_360 = mtAngsFiSt_gr(setor,:);     % grupo de angulos de steering [0º, 360] disponiveis para o setor 'setor'
    angFiGrupo_180 = wrapTo180(angFiGrupo_360);  % grupo de angulos de steering [-180º, +180º] disponiveis para o setor 'setor'
    
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  % angulo azimutal do UE atual em [-180, 180]
    anguloUE360 = wrapTo360(anguloUE180);        % angulo azimutal do UE atual em [0, 360]
    
    %dif1 = min(sqrt(((anguloUE360 - angFiGrupo_360).^2)));
    dif1 = min(abs(anguloUE360 - angFiGrupo_360));
    %dif2 = min(sqrt(((anguloUE180 - angFiGrupo_180).^2)));
    dif2 = min(abs(anguloUE180 - angFiGrupo_180));
    difAngulo = min(dif1, dif2);
    
    mtdifAngsHor_gr(setor, jj) = difAngulo;
    
    % laço percorrendo os SETORES 
    for s = 1:S,

        % se o setor 's' é diferente do setor do UEjj (laço externo), então 
        if s ~= setor,
            
            % se não temos usuário ativo no setor 's' no instante 'inst', então
            if mtUeSector(s, inst) == 0,
                ddif1 = min(abs(anguloUE360 - ((180/pi).*vtAngST(s))));
                ddif2 = min(abs(anguloUE180 - wrapTo180(((180/pi).*vtAngST(s)))));
                mtdifAngsHor_gr(s, jj) = min(ddif1, ddif2);
            else
                angUeAtivo = wrapTo360((180/pi).*angle(vtUePos(mtUeSector(s, inst)))); % UE ativo no instante 'inst' no setor 'setor'
                [M, II] = min(abs(angUeAtivo - mtAngsFiSt_gr(s,:)));
                ddif1 = min(abs(anguloUE360 - mtAngsFiSt_gr(s,II)));
                ddif2 = min(abs(anguloUE180 - wrapTo180(mtAngsFiSt_gr(s,II))));
                mtdifAngsHor_gr(s, jj) = min(ddif1, ddif2);
            end
        end
    end
end


% diferença de angulos de elevação de cada UE de grupo
mtdifAngsVer_gr = zeros(S, numUE);

% laço percorrendo todos UE's
for jj= 1:numUE,
    
    % 'setor' que o UEjj pertence e o instante 'inst' que o UEjj está ativo
    [setor, inst] = find(mtUeSector == jj);       
    
    % calculo o angulo de elevação THETA do UE
    angElevUE180 = mtThetaUE(setor, jj);
    
    mtdifAngsVer_gr(setor, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(setor,:)));
    
    % laço percorrendo todos os  SETORES
    for s = 1:S,
        
        % Percorrendo os OUTROS setores, i.e, setores diferente do setor que contém o UEjj (laço externo)
        if s ~= setor, 
            
            % se NÃO temos UE ATIVO no setor 's' no instante 'inst', então
            if mtUeSector(s, inst) == 0,
                
                % diferença entre o Angulo de Elevação do UEjj e o angulo de Downtild do setor será
                mtdifAngsVer_gr(s, jj) = min(abs(angElevUE180 - angDownTild_2d));
            
            % se temos UE ATIVO no setor 's' no instante 'inst', então    
            else
                angElevUeAtivo = mtThetaUE(s, mtUeSector(s, inst));    % ang. de ELEVAÇÂO do UE ativo no setor 's' no instante 'inst'
                [M, II] = min(abs(angElevUeAtivo - mtAngsThetaTild_gr(s,:)));
                mtdifAngsVer_gr(s, jj) = min(abs(angElevUE180 - mtAngsThetaTild_gr(s, II)));
            end
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
    H_GR(:,:,ii) = G_BS + A_GR(:,:,ii) + mtPL + mtNormal;   % [dB]
end


% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_gr = db2lin(H_GR);

% SINR para Beamforming GRUPO
Y_GR = [];   % SINR

% laço percorrendo cada dimensao do tensor p/ calcular cada dimensão do setor
for jj = 1:length(vtFi3dB_Esp),
  
    % laço percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % soma todo os coeficientes de canal até o UE menos o coef. do canal do setor que ele se encontra
        aux = sum(h_gr(:, ii, jj)) - h_gr(I(ii), ii, jj);
        
        % calculo da SNIR, onde cada linha será Y_ESP p/ Fi3dB e theta3dB
        Y_GR(jj, ii) = (Pot*h_gr(I(ii), ii, jj))/(Pot*aux + PN);     % linha: SNIR p/ cada \theta_3dB e \fi_3dB
        
    end
end

% SINR em dB
YGR_dB = 10*log10(Y_GR);
cdfplot(YGR_dB(1,:))
legend('Conventional', 'UE especifica', '16 grupo');


