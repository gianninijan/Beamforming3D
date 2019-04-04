%% (ARTIGO1) %%
% 3D Beamforming Capacity Improvement in Macrocell-Assisted Small Cell Architeture

clear all;
clc;
close all;


%% SETUP SIMULATION  

M = 1;                                               % numero de celulas
FatorSetor = 3;                                      % Fator de setoriza��o, i.e, setores/celulas
S = M*FatorSetor;                                    % n�mero de setores. S = {1, 2, 3, ..., }      
numUE = 10;                                        % numero de UE's por macro-setores
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

% angUE = [pi/3, pi, 5*pi/3];             % angulos de teste
% raioUE = (R/2).*ones(1,numUE);          % raios de teste
% vtUePos = raioUE.*exp(1j*angUE);        % vetor de posi��o das UE's de teste

vtUePos = [];                             % vetor de posi��o de usu�rios

% GERANDO AS POSI��ES DO UE's ~ U(10, R)
for ii = 1:numUE,
    while true,
         vtUePos(ii) = (2*R*rand(1,1) - R) + 1j.*(2*R*rand(1,1) - R);
         
         if (abs(vtUePos(ii)) < R) && (abs(vtUePos(ii)) > 10),
             break;
         end
    end
end

% vetor de posi��o de testes
% vtUePos = [-71.1058-28.4379i  35.3306+37.7640i  44.0459-63.8599i  43.3912+44.9478i -74.6801-15.4539i];

% vetor das posi��es da BS de cada setor
raio_bs = 0;                                        % raio das BS's - gerar aleatorio
ang_bs = 0;                                         % angulos de posi��es das BS's - gerar aleatorio
posBS = raio_bs.*exp(1j*ang_bs);                    % forma exponencial das BS's
vtBsSetor = repelem(posBS, FatorSetor);             % vetor posi��o da BS's p/ cada setor


% plotando as POSI��ES da BS e UE's    
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
sector1 = find( (angUE >= 0) & (angUE < (2*pi/3)) );          % indices dos UE's que est�o no setor 1
sector2 = find( (angUE >= (2*pi/3)) & (angUE < (4*pi/3)) );   % indices dos UE's que est�o no setor 2
sector3 = find( (angUE >= (4*pi/3)) & (angUE < (2*pi)) );     % indices dos UE's que est�o no setor 3

% Testando o codigo acima dos setores: 
% resultado = vtUePos(sector2);          % captura as distancia dos UE's que est�o dentro do setor 1
% plot(resultado,'ro');                  % plota o valores do 'resultados' com um circulo azul  
% 
% hold off


%% CALCULAR O PATH-LOSS DOS UE'S

% matriz de DIST�NCIA de cada UE (coluna) para BS (linha) de cada setor
mtDist = zeros(S,numUE);                             

% la�o percorrendo cada setor (coluna)
for ii = 1:S
    
    % la�o percorrendo todos os usu�rios (linha)
    for jj = 1:numUE
        mtDist(ii,jj) = norm(vtUePos(jj) - vtBsSetor(ii));      % DIM( mtDist ) = [Linhas, Colunas] = [Setores, UE's]
    end
    
end

% P/ calcular o setor de cada UE devemos levar em conta a distancia p/ cada BS do setor e angulo
mtPL = PathLoss(mtDist, fc);                  % matriz de PATH LOSS de cada setor para cada usu�rio
PL = min(mtPL);                               % valor minimo

vtDistUEtoBS = min(mtDist);
%[vtDistUEtoBS, POS] = min(mtDist, [], 1);
figure;                                       % gera uma nova figura para plotar os graficos
plot(sort(vtDistUEtoBS), sort(PL))
xlabel('d (M)')
ylabel('PL (DB)')
title('Path Loss ')


%% DADOS COMUNS PARA OS BEAMFORMINGS

% vetor de DESVANECIMENTO R�PIDO para cada UE's
% mtFastFad_2D = (1/sqrt(2))*[randn(S, numUE) + 1j*randn(S, numUE)];
% load('desvanecimento.mat');

% matriz de shadowing normal
mtNormal = sigma.*randn(S,numUE);    % dim( mtNormal ) = [Linhas, Colunas] = [Setores, Usu�rios]

% Angulo \theta de cada UE p/ cada BS do setor 
mtThetaUE = atand((H_BS - H_UE)./(mtDist));   % [GRAUS]

% Pot�ncia do Ru�do Linear
PN = dbm2lin(No+ 10*log10(bandWidth)+FN);  

% Pot�ncia da antena transmissora
Pot = dbm2lin(P_BS);                       


%%  BEAMFORMING 2D %%

% angulos de BORESIGHT FIXO p/ a antena de cada setor, i.e, angulo azimutal na qual teremos o ganho m�ximo da antena 
ang_st = [pi/3, pi, 5*pi/3];                        % angulos de boresight p/ cada celula [Radianos]
vtAngST = repmat(ang_st, 1, M);                     % angulo steering de cada setor [Radianos]

% valores tirados do artigo
fi3dB_2D = 70;                     % largura de feixe de 3 dB na horizontal [GRAUS]
theta3dB_2D = 10;                  % largura de feixe de 3 dB na vertical   [GRAUS]
angDownTild_2d = 8;                % angulo de down-tild (FIXO) [GRAUS]

% calculando o padr�o vertical da antena
Av_2D = -min(12.*(((mtThetaUE - angDownTild_2d)./theta3dB_2D).^2), SLA);    % [linhas, colunas] = [setores, UE's]

% calculando o padr�o horizontal da antena
Ah_2D = [];

% posi��es dos UE's 
% posicoes_d = (180/pi).*angle(vtUePos); % [-180, 180]

% diferen�a do angulo azimutal de cada UE p/ cada angulo de steering de cada setor
mtdifAngulos_d = [];

% la�o percorrendo cada UE
for jj = 1:numUE,
    
    % angulo de cada UE em [-180�, 180�]
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  
   
    % angulo de cada UE em [0�, 360�]
    anguloUE360 = wrapTo360(anguloUE180);
    
    % calculo o angulo de steering para cada setor entre [0�, 360�]
    angStrTo360 = (180/pi).*vtAngST;    % em, graus (�)
    
    % calculo o angulo de steering para cada setor entre [-180�, 180�]
    angStrTo180 = wrapTo180(angStrTo360);

    % la�o percorrendo cada setor
    for ii = 1:S,
        
        % calcula a diferen�a entre angulos
        dif1 = abs(anguloUE360 - angStrTo360(ii));
        dif2 = abs(anguloUE180 - angStrTo180(ii));
        difAngulo = min(dif1, dif2);
        
        mtdifAngulos_d(ii, jj) = difAngulo;
    end
end

% padr�o de radia��o HORIZONTAL
Ah_2D = -min(12.*((mtdifAngulos_d./fi3dB_2D).^2), Am);

% Padr�o de radia��o TOTAL
A_2D = -min(-(Ah_2D + Av_2D), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_2d = G_BS + A_2D + mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_2d = db2lin(H_2d);

% SINR
Y2D = [];   % SINR

% capturando o setor de cada UE atraves do INDICES
[X,I] = min(sqrt(mtdifAngulos_d.^2), [], 1);

% la�o percorrendo cada UE para calcular a SINR
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
vtFi3dB_Esp = [70 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_Esp = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]

% angulo de donwtild p/ Beamforming UE especifico ser� igual a \theta de cada usu�rio
angDownTild_Esp = [];

minimo = min([length(sector1), length(sector2), length(sector3)]);

% Calculando o angulo de downtild
% la�o percorrendo os usu�rios at� o minimo valor 
for jj = 1:S*minimo,
    
    % se UEj pertence ao setor 1, ent�o: 
    if find(sector1 == jj),
        angDownTild_Esp(1,jj) = atand((H_BS - H_UE)./(norm(vtUePos(jj) - vtBsSetor(1))));
        angDownTild_Esp(2,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector2(jj)) - vtBsSetor(2))));
        angDownTild_Esp(3,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector3(jj)) - vtBsSetor(3))));
    
    % se UEj pertence ao setor 2, ent�o:     
    elseif find(sector2 == jj),
        angDownTild_Esp(1,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector1(jj)) - vtBsSetor(1))));
        angDownTild_Esp(2,jj) = atand((H_BS - H_UE)./(norm(vtUePos(jj) - vtBsSetor(2))));
        angDownTild_Esp(3,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector3(jj)) - vtBsSetor(3))));
        
    % se UEj pertence ao setor 3, ent�o:     
    elseif find(sector3 == jj),
        angDownTild_Esp(1,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector1(jj)) - vtBsSetor(1))));
        angDownTild_Esp(2,jj) = atand((H_BS - H_UE)./(norm(vtUePos(sector2(jj)) - vtBsSetor(2))));
        angDownTild_Esp(3,jj) = atand((H_BS - H_UE)./(norm(vtUePos(jj) - vtBsSetor(3))));
    end
    
    
end

% TENSOR para calcular o padr�o de RADIA��O VERTICAL da antena para cada ANGULO \theta_3dB
Av_Esp = [];

% la�o percorrendo cada angulo theta_3dB
for ii = 1:length(vtTheta3dB_Esp),
    
    Av_Esp(:,:,ii) = -min(12.*(((mtThetaUE(:,1:minimo) - angDownTild_Esp)./vtTheta3dB_Esp(ii)).^2), SLA);  

end

% tensor para calcular o padr�o de radia��o horizontal da antena para cada angulo \theta_3dB
Ah_Esp = [];

% diferen�a do angulo azimutal de cada UE p/ cada angulo de steering de cada setor
mtdifAnguloEsp_d = mtdifAngulos_d;   % linhas: setores, colunas: UEs

% la�o percorrendo cada UE's
for jj = 1:numUE,
    
    % mudando apenas a diferen�a de angulos na horizontal p/ o setor de cada UE
    mtdifAnguloEsp_d( I(jj), jj) = 0;
    
end


% la�o percorrendo cada angulo fi_3dB
for ii = 1:length(vtFi3dB_Esp),

    % padr�o de radia��o HORIZONTAL
    Ah_Esp(:,:,ii) = -min(12.*((mtdifAnguloEsp_d./vtFi3dB_Esp(ii)).^2), Am);
    
end

% Padr�o de radia��o TOTAL
A_ESP = -min(-(Ah_Esp + Av_Esp), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_ESP = G_BS + A_ESP + mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_esp = db2lin(H_ESP);

% SINR
Y_ESP = [];   % SINR

% la�o percorrendo cada dimensao do tensor
for jj = 1:length(vtFi3dB_Esp),
    
    % la�o percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        aux = sum(h_esp(:, ii, jj)) - h_esp(I(ii), ii, jj);
        
        % calculo da SNIR, onde cada linha ser� Y_ESP p/ Fi3dB e theta3dB
        Y_ESP(jj, ii) = (Pot*h_esp(I(ii), ii, jj))/(Pot*aux + PN);
        
    end
    
end

% SINR em dB
YESP_dB = 10*log10(Y_ESP);

% hold on;
% figure;
cdfplot(YESP_dB(1, :))
% cdfplot(YESP_dB(2, :))
% cdfplot(YESP_dB(3, :))
legend('Conventional', 'UE especifica - (\theta_{3dB}, \phi_{3dB) = (70�, 10�)}');


%% BEAMFORMING GRUPO-ESPECIFICO

% valores tirados do artigo 
vtFi3dB_gr = [70 10 5];            % largura de feixe de 3 dB na horizontal [GRAUS]
vtTheta3dB_gr = [10 10 5];         % largura de feixe de 3 dB na vertical   [GRAUS]
B = 16;                            % numero de padr�es de feixes (ou grupos)
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
mtAngsFiSt_gr = mtAngsFiSt_gr';                  % linha: setor, coluna: \fi_st 's p/ cada setor

% angDownTild_gr = linspace(-90, 90, Bv);
angDownTild_gr = [4 8];

% matriz de angulos \theta_thild p/ uma celula;
mtAngsThetaTild_gr = repmat(angDownTild_gr, FatorSetor*M, 1); % linha: setor, coluna: \theta_tild 's p/ cada setor

% diferen�a do angulo azimutal de cada UE de grupo
mtdifAngsHor_gr = zeros(S, numUE);

% la�o percorrendo todos UE's
for jj= 1:numUE,

    % angulo azimutal de cada UE em [-180, 180]
    anguloUE180 = (180/pi).*angle(vtUePos(jj));  
   
    % angulo azimutal de cada UE em [0, 360]
    anguloUE360 = wrapTo360(anguloUE180);
    
    % calculo o angulo de steering para cada setor entre [0, 360]
    angStrTo360 = (180/pi).*vtAngST;
    
    % calculo o angulo de steering para cada setor entre [-180, 180]
    angStrTo180 = wrapTo180(angStrTo360);
    
    % la�o percorrendo cada setor
    for ii = 1:S,
        
        % se o setor da itera��o for igual ao correspondente setor do UE, ent�o
        if ii == I(jj) 
            % calcula a diferen�a absoluta entre o UE e \fi_st mais proximo 
            vtAux = sqrt((anguloUE360 - mtAngsFiSt_gr(ii, :)).^2);
            mtdifAngsHor_gr(ii, jj) = min(vtAux);    
        else
            % calcula a diferen�a entre angulos
            dif1 = sqrt(((anguloUE360 - angStrTo360(ii)).^2));
            dif2 = sqrt(((anguloUE180 - angStrTo180(ii)).^2));
            difAngulo = min(dif1, dif2);
            mtdifAngsHor_gr(ii, jj) = difAngulo;
        end    
    end  
end

% diferen�a de angulos de eleva��o de cada UE de grupo
mtdifAngsVer_gr = zeros(S, numUE);

% la�o percorrendo todos UE's
for jj= 1:numUE,
    
    % la�o percorrendo cada setor
    for ii = 1:S,
        
        % diferen�a entre o angulos de downtild
        mtdifAngsVer_gr(ii, jj) = min(abs(mtThetaUE(ii,jj) - mtAngsThetaTild_gr(ii,:)));
    end
    
end

% tensor para calcular o padr�o de radia��o vertical da antena para cada angulo \theta_3dB
Ah_gr = [];

% calculando o padr�o de radia��o vertical
for ii = 1:length(vtFi3dB_gr),
    
    % Ah para cada \fi_3dB
    Ah_gr(:,:,ii) = -min(12.*((mtdifAngsHor_gr./vtFi3dB_gr(ii)).^2), Am);
    
end

% tensor para calcular o padr�o de radia��o vertical da antena para cada angulo \theta_3dB
Av_gr = [];

% calculando o padr�o de radia��o vertical
for ii = 1:length(vtTheta3dB_gr),
    
    % Av para cada \theta_3dB
    Av_gr(:,:,ii) = -min(12.*((mtdifAngsVer_gr./vtTheta3dB_gr(ii)).^2), SLA);
      
end


% Padr�o de radia��o TOTAL
A_GR = -min(-(Av_gr + Ah_gr), Am);

% coeficientes do canal ao quadrado em dB (sem fast-fading)
H_GR = G_BS + A_GR + mtPL + mtNormal;   % [dB]

% coeficientes do canal ao quadrado em escala linear (sem fast-fading)
h_gr = db2lin(H_GR);

% SINR
Y_GR = [];   % SINR

% la�o percorrendo cada dimensao do tensor
for jj = 1:length(vtFi3dB_Esp),
  
    % la�o percorrendo cada UE para calcular a SINR
    for ii = 1:numUE,
        
        % soma todo os coeficientes de canal at� o UE menos o coef. do
        % canal do setor que ele se encontra
        aux = sum(h_gr(:, ii, jj)) - h_gr(I(ii), ii, jj);
        
        % calculo da SNIR, onde cada linha ser� Y_ESP p/ Fi3dB e theta3dB
        Y_GR(jj, ii) = (Pot*h_gr(I(ii), ii, jj))/(Pot*aux + PN);     % linha: SNIR p/ cada \theta_3dB e \fi_3dB
        
    end
end

% SINR em dB
YGR_dB = 10*log10(Y_GR);
cdfplot(YGR_dB(1,:))
legend('Conventional', 'UE especifica', '16 grupo');


