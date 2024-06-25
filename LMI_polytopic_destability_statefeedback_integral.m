% Controle de LMIs
% José Sergio Cruz Dantas Junior

% Tarefa 1: statefeedback control + integral action + D-stability + polytopic uncertain

clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% Modelo do Processo
% H = tf(9,[1 4.8 9]); % processo estável
H = tf(9,[1 -4.8 9]); % processo instável
[num,den] = tfdata(H,'v');
[A,B,C,D] = tf2ss(num,den);

delta = [0.2 0.4
         0.1  0.2];
A1 = A + delta;
A2 = A - delta;

% return
%% Sistema Aumentado
AA1 = [A1 zeros(length(A1),1); -C 0];
AA2 = [A2 zeros(length(A2),1); -C 0];
BB = [B; -D];

G1 = ss(A1,B,C,D);
OL_poles_G1 = eig(A1);
G2 = ss(A2,B,C,D);
OL_poles_G2 = eig(A2);

% return
%% Parâmetros de Desempenho
t_min = 2; % min settling time
t_max = 5; % max settling time
UP = 0.1; % overshoot

alfa = -4/t_max; % Vertical Strip
beta = -4/t_min; % Vertical Strip
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
theta = acos(zeta); % Conic Sector angle

Gtf = tf([4],[1 4*zeta 4]);

% return
%% LMI
Z = sdpvar(1,length(AA1));
P = sdpvar(length(AA1));

% LMI0 = [P >= 0, Z >= 0];
LMI0 = [P >= 0];
% LMIs para A1
LMI1 = [[sin(theta)*(AA1*P + BB*Z + P*AA1' + Z'*BB'),cos(theta)*(AA1*P + BB*Z - P*AA1' - Z'*BB');
       cos(theta)*(P'*AA1' + Z'*BB' - AA1*P - BB*Z),sin(theta)*(AA1*P + BB*Z + P*AA1' + Z'*BB')] <= 0];
LMI2 = [AA1*P + P*AA1' + BB*Z + Z'*BB' - 2*alfa*P <=0];
LMI3 = [AA1*P + P*AA1' + BB*Z + Z'*BB' - 2*beta*P >= 0];
% LMIs para A2
LMI4 = [[sin(theta)*(AA2*P + BB*Z + P*AA2' + Z'*BB'),cos(theta)*(AA2*P + BB*Z - P*AA2' - Z'*BB');
       cos(theta)*(P'*AA2' + Z'*BB' - AA2*P - BB*Z),sin(theta)*(AA2*P + BB*Z + P*AA2' + Z'*BB')] <= 0];
LMI5 = [AA2*P + P*AA2' + BB*Z + Z'*BB' - 2*alfa*P <=0];
LMI6 = [AA2*P + P*AA2' + BB*Z + Z'*BB' - 2*beta*P >= 0];
LMI = [LMI0,LMI1,LMI2,LMI3,LMI4,LMI5,LMI6]; 
optimize(LMI);
checkset(LMI);

Z = value(Z);
P = value(P);
K = Z*inv(P); % controlador
Kx = K(1:end-1)
Ki = K(end)

Amf1 = AA1 + BB*K;
Amf2 = AA2 + BB*K;

% Gmf1 = ss(Amf1,B,C,D);
CL_poles_G1 = eig(Amf1);
% Gmf2 = ss(Amf2,B,C,D);
CL_poles_G2 = eig(Amf2);

% U1 = ss(Amf1,B,K,0);
% U2 = ss(Amf2,B,K,0);

% return
%% Figuras
t1 = -10:0.01:10;
t2 = -10:0.01:0;
% l1 = -theta*t2*2.162;
l1 = -tan(theta)*t2;
% l2 = theta*t2*1.454;
l2 = tan(theta)*t2;
h1=alfa*ones(size(t1));
h2=beta*ones(size(t1));

% area d-estabilidade (intercessoes entre as linhas)
aux1 = beta*tan(theta);
aux2 = alfa*tan(theta);
pgon = polyshape([beta beta alfa alfa],[aux1 -aux1 -aux2 aux2]);

% Parametros de Simulaçao
Tsim = 20;
ref = 1;
Tref = 0;
Ts = 1;

A = A1;
out = sim('lqi_simulink.slx');

figure(1)
hold on
subplot(3,2,1)
plot(out.time,out.yOL,'b','LineWidth',2)
% step(G1)
title('$G_1$ open loop')

subplot(3,2,3)
plot(out.time,out.yCL,'b','LineWidth',2)
title('$G_1$ closed loop')

subplot(3,2,5)
plot(out.time,out.u,'b','LineWidth',2)
title('$U_1$ control signal')

subplot(3,2,[2 4 6])
hold on
scatter(real(OL_poles_G1),imag(OL_poles_G1),'bx','LineWidth',2)
scatter(real(CL_poles_G1),imag(CL_poles_G1),'rx','LineWidth',2)
% scatter(real(eig(Gtf)),imag(eig(Gtf)),'gx','LineWidth',2)
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
plot(t2,l1,'k--','LineWidth', 2)
plot(t2,l2,'k--','LineWidth',2)
plot(h1,t1,'k--','LineWidth',2)
plot(h2,t1,'k--','LineWidth',2)
rlocus(tf(1,1))
axis([-3 3 -3 3])
legend('Polos em MA','Polos em MF','D-stability area')
grid on

A = A2;
out = sim('lqi_simulink.slx');

figure(2)
hold on
subplot(3,2,1)
plot(out.time,out.yOL,'b','LineWidth',2)
% step(G2)
title('$G_2$ open loop')

subplot(3,2,3)
plot(out.time,out.yCL,'b','LineWidth',2)
title('$G_2$ closed loop')

subplot(3,2,5)
plot(out.time,out.u,'b','LineWidth',2)
title('$U_2$ control signal')

subplot(3,2,[2 4 6])
hold on
scatter(real(OL_poles_G2),imag(OL_poles_G2),'bx','LineWidth',2)
scatter(real(CL_poles_G2),imag(CL_poles_G2),'rx','LineWidth',2)
% scatter(real(eig(Gtf)),imag(eig(Gtf)),'gx','LineWidth',2)
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
plot(t2,l1,'k--','LineWidth', 2)
plot(t2,l2,'k--','LineWidth',2)
plot(h1,t1,'k--','LineWidth',2)
plot(h2,t1,'k--','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MA','Polos em MF','D-stability area')
axis([-3 3 -3 3])
grid on


%fim