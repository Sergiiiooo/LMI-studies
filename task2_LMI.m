% ANÁLISE E CONTROLE DE SISTEMAS LINEARES POR DESIGUALDADES MATRICIAIS LINEARES (LMIS)
% Aluno: José Sergio Cruz Dantas Junior

% Tarefa 2: statefeedback control + integral action + D-stability +
% polytopic uncertain + H2 control + Hinf control + H2/Hinf control

% Pitch control for wind turbine systems
% http://dx.doi.org/10.1016/j.renene.2016.01.057

clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% Modelo do Processo
% Wind turbine parameters (WT1)
a2 = 2.426; a1 = -4.6345; a0 = -147.3;
b4 = 1; b3 = 4.857; b2 = 126.2; b1 = 266.4; b0 = 3659;
A1 = [-b3 -b2 -b1 -b0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;];
Bu1 = [1 0 0 0]';
Cy1 = [0 a2 a1 a0];
Dy1 = 0;
G1 = ss(A1,Bu1,Cy1,Dy1);

% deltaA = [[-b3 -b2 -b1 -b0]*0.1;
%     0 0 0 0;
%     0 0 0 0;
%     0 0 0 0;];

% deltaC = [0 a2 a1 a0]*0.1;

% Wind turbine parameters (WT2)
a2 = -0.6219*-2; a1 = -8.7165; a0 = -2911*0.1;
b4 = 1; b3 = 5.018; b2 = 691.3*0.5; b1 = 1949*0.1; b0 = 1.15e5*0.05;
A2 = [-b3 -b2 -b1 -b0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;];
% A2 = A1+deltaA;
Bu2 = [1 0 0 0]';
Cy2 = [0 a2 a1 a0];
% Cy2 = Cy1+deltaC;
Dy2 = 0;
G2 = ss(A2,Bu2,Cy2,Dy2);

% Wind turbine parameters (WT3)
a2 = -0.2545*-1; a1 = -0.0647*30; a0 = 0.9384*-200;
b4 = 1; b3 = 2.28; b2 = 878.5*0.5; b1 = 437.7*0.8; b0 = 7.7e4*0.06;
A3 = [-b3 -b2 -b1 -b0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;];
% A3 = A1+2*deltaA;
Bu3 = [1 0 0 0]';
Cy3 = [0 a2 a1 a0];
% Cy3 = Cy1-deltaC;
Dy3 = 0;
G3 = ss(A3,Bu3,Cy3,Dy3);

% Polos em Malha Aberta
polos_MA_WT1 = eig(A1)
polos_MA_WT2 = eig(A2)
polos_MA_WT3 = eig(A3)

figure
step(G1,G2,G3)
title('Malha aberta')
xlabel('Tempo')
ylabel('Amplitude')
legend('WT1','WT2','WT3')

figure
hold on
scatter(real(polos_MA_WT1),imag(polos_MA_WT1),'bx','LineWidth',2)
scatter(real(polos_MA_WT2),imag(polos_MA_WT2),'gx','LineWidth',2)
scatter(real(polos_MA_WT3),imag(polos_MA_WT3),'mx','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MA WT1','Polos em MA WT2','Polos em MA WT3')
grid on

return
%% Lyapunov stability
% P = sdpvar(length(A),length(A));
% 
% F1 = [P>=0];
% F2 = [A'*P+P*A<=0];
% F = [F1,F2]
% optimize(F)
% checkset(F)
% 
% P = value(P);
% return
%% State-feedback control
P = sdpvar(length(A1),length(A1));
Y = sdpvar(1,length(A1));

F1 = [P>=0]; 
F2 = [P*A1'+Y'*Bu1'+A1*P+Bu1*Y<=0];
F3 = [P*A2'+Y'*Bu2'+A2*P+Bu2*Y<=0];
F4 = [P*A3'+Y'*Bu3'+A3*P+Bu3*Y<=0];
F = [F1,F2,F3,F4];
optimize(F);
checkset(F)

P = value(P);
Y = value(Y);
K = Y*inv(P);

Gmf_1 = ss(A1+Bu1*K,Bu1,Cy1,Dy1);
Gmf_2 = ss(A2+Bu2*K,Bu2,Cy2,Dy2);
Gmf_3 = ss(A3+Bu3*K,Bu3,Cy3,Dy3);

polos_MF_LQ_WT1 = eig(A1+Bu1*K)
polos_MF_LQ_WT2 = eig(A2+Bu2*K)
polos_MF_LQ_WT3 = eig(A3+Bu3*K)

figure
step(Gmf_1,Gmf_2,Gmf_3)
title('Malha fechada')
xlabel('Tempo')
ylabel('Amplitude')

figure
hold on
scatter(real(polos_MF_LQ_WT1),imag(polos_MF_LQ_WT1),'bx','LineWidth',2)
scatter(real(polos_MF_LQ_WT2),imag(polos_MF_LQ_WT2),'gx','LineWidth',2)
scatter(real(polos_MF_LQ_WT3),imag(polos_MF_LQ_WT3),'mx','LineWidth',2)
% scatter(real(polos_MA),imag(polos_MA),'bx','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MF_LQ WT1','Polos em MF_LQ WT2','Polos em MF_LQ WT3')
grid on

return
%% State-feedback + integral control
AA1 = [A1 zeros(size(A1,1),size(Cy1,1)); -Cy1 zeros(size(Cy1,1),size(Cy1,1))];
BBu1 = [Bu1; -Dy1];
AA2 = [A2 zeros(size(A2,1),size(Cy2,1)); -Cy2 zeros(size(Cy2,1),size(Cy2,1))];
BBu2 = [Bu2; -Dy2];
AA3 = [A3 zeros(size(A3,1),size(Cy3,1)); -Cy3 zeros(size(Cy3,1),size(Cy3,1))];
BBu3 = [Bu3; -Dy3];

Q = sdpvar(length(AA1),length(AA1),'symmetric');
Y = sdpvar(1,length(BBu1));

LMI1 = [Q>=0];
LMI2 = [Q*AA1'+AA1*Q+Y'*BBu1'+BBu1*Y<=0];
LMI3 = [Q*AA2'+AA2*Q+Y'*BBu2'+BBu2*Y<=0];
LMI4 = [Q*AA3'+AA3*Q+Y'*BBu3'+BBu3*Y<=0];
LMI = [LMI1,LMI2,LMI3,LMI4];
optimize(LMI);
checkset(LMI)

Q = value(Q);
Y = value(Y);

K = Y*inv(Q);
Kx = K(1:length(A1));
Ki = K(length(A1)+1:end);

Amf_1 = AA1 + BBu1*K;
Amf_2 = AA2 + BBu2*K;
Amf_3 = AA3 + BBu3*K;

polos_MF_LQI_WT1 = eig(Amf_1)
polos_MF_LQI_WT2 = eig(Amf_2)
polos_MF_LQI_WT3 = eig(Amf_3)

figure
hold on
scatter(real(polos_MF_LQI_WT1),imag(polos_MF_LQI_WT1),'rx','LineWidth',2)
scatter(real(polos_MA_WT1),imag(polos_MA_WT1),'bx','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MF_LQI','Polos em MA')
grid on
title('Polos WT1')

figure
hold on
scatter(real(polos_MF_LQI_WT2),imag(polos_MF_LQI_WT2),'rx','LineWidth',2)
scatter(real(polos_MA_WT2),imag(polos_MA_WT2),'bx','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MF_LQI','Polos em MA')
grid on
title('Polos WT2')

figure
hold on
scatter(real(polos_MF_LQI_WT3),imag(polos_MF_LQI_WT3),'rx','LineWidth',2)
scatter(real(polos_MA_WT3),imag(polos_MA_WT3),'bx','LineWidth',2)
rlocus(tf(1,1))
legend('Polos em MF_LQI','Polos em MA')
grid on
title('Polos WT3')
return
%% Simulação
% Parametros de Simulaçao
Tsim = 20;
ref = 1;
Tref = 0;
Ts = .1;
Tpert = Tsim/2;
pert = 10;

% G1
A=A1;B=Bu1;C=Cy1;D=Dy1;
out = sim('task1_simulink.slx');

% figure
% plot(out.time,out.yOL,'b','LineWidth',1)
% title('Malha Aberta WT1')
% xlabel('Tempo')
% ylabel('Amplitude')

figure
subplot(1,3,1)
plot(out.time,out.yCL,'r','LineWidth',2)
title('Malha Fechada WT1')
xlabel('Tempo')
ylabel('Amplitude')

% figure
% plot(out.time,out.u,'g','LineWidth',2)
% title('Sinal de controle WT1')
% xlabel('Tempo')
% ylabel('Amplitude')

% G2
A=A2;B=Bu2;C=Cy2;D=Dy2;
out = sim('task1_simulink.slx');

% figure
% plot(out.time,out.yOL,'b','LineWidth',1)
% title('Malha Aberta WT2')
% xlabel('Tempo')
% ylabel('Amplitude')

% figure
subplot(1,3,2)
plot(out.time,out.yCL,'r','LineWidth',2)
title('Malha Fechada WT2')
xlabel('Tempo')
ylabel('Amplitude')

% figure
% plot(out.time,out.u,'g','LineWidth',2)
% title('Sinal de controle WT2')
% xlabel('Tempo')
% ylabel('Amplitude')

% G3
A=A3;B=Bu3;C=Cy3;D=Dy3;
out = sim('task1_simulink.slx');

% figure
% plot(out.time,out.yOL,'b','LineWidth',1)
% title('Malha Aberta WT3')
% xlabel('Tempo')
% ylabel('Amplitude')

% figure
subplot(1,3,3)
plot(out.time,out.yCL,'r','LineWidth',2)
title('Malha Fechada WT3')
xlabel('Tempo')
ylabel('Amplitude')

% figure
% plot(out.time,out.u,'g','LineWidth',2)
% title('Sinal de controle WT3')
% xlabel('Tempo')
% ylabel('Amplitude')
sgtitle("H2/H∞ control + Statefeedback control + Integral action + D-stability")
% return
%% Parâmetros de Desempenho
t_min = 0.1 % min settling time
t_max = 2 % max settling time
UP = 0.20 % overshoot

alfa = -4/t_max % Vertical Strip
beta = -4/t_min % Vertical Strip
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
theta = acos(zeta) % Conic Sector angle

% Linhas D-estabilidade
t1 = -100:0.01:100;
t2 = -100:0.01:0;
l1 = -tan(theta)*t2;
l2 = tan(theta)*t2;
h1=alfa*ones(size(t1));
h2=beta*ones(size(t1));

% Regiao D-estabilidade (intercessoes entre as linhas)
aux1 = beta*tan(theta);
aux2 = alfa*tan(theta);
pgon = polyshape([beta beta alfa alfa],[aux1 -aux1 -aux2 aux2]);

figure
hold on
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
rlocus(tf(1,1))
plot(t2,l1,'k','LineWidth', 3)
plot(t2,l2,'k','LineWidth',3)
plot(h1,t1,'k','LineWidth',3)
plot(h2,t1,'k','LineWidth',3)
axis([-45 10 -100 100])
grid on
title('Regiao de d-estabilidade')

% return
%% State-feedback + integral + D-stability control
AA1 = [A1 zeros(size(A1,1),size(Cy1,1)); -Cy1 zeros(size(Cy1,1),size(Cy1,1))];
BBu1 = [Bu1; -Dy1];
AA2 = [A2 zeros(size(A2,1),size(Cy2,1)); -Cy2 zeros(size(Cy2,1),size(Cy2,1))];
BBu2 = [Bu2; -Dy2];
AA3 = [A3 zeros(size(A3,1),size(Cy3,1)); -Cy3 zeros(size(Cy3,1),size(Cy3,1))];
BBu3 = [Bu3; -Dy3];

P = sdpvar(size(AA1,1),size(AA1,1),'symmetric');
Z = sdpvar(size(BBu1,2),size(AA1,1));

LMI0 = [P >= 0];
% WT1
LMI1 = [[sin(theta)*(AA1*P + BBu1*Z + P*AA1' + Z'*BBu1'),cos(theta)*(AA1*P + BBu1*Z - P*AA1' - Z'*BBu1');
       cos(theta)*(P*AA1' + Z'*BBu1' - AA1*P - BBu1*Z),sin(theta)*(AA1*P + BBu1*Z + P*AA1' + Z'*BBu1')] <= 0];
LMI2 = [AA1*P + P*AA1' + BBu1*Z + Z'*BBu1' - 2*alfa*P <= 0]; % Vertical restriction alfa
LMI3 = [AA1*P + P*AA1' + BBu1*Z + Z'*BBu1' - 2*beta*P >= 0]; % Vertical restriction beta
% WT2
LMI4 = [[sin(theta)*(AA2*P + BBu2*Z + P*AA2' + Z'*BBu2'),cos(theta)*(AA2*P + BBu2*Z - P*AA2' - Z'*BBu2');
       cos(theta)*(P*AA2' + Z'*BBu2' - AA2*P - BBu2*Z),sin(theta)*(AA2*P + BBu2*Z + P*AA2' + Z'*BBu2')] <= 0];
LMI5 = [AA2*P + P*AA2' + BBu2*Z + Z'*BBu2' - 2*alfa*P <= 0]; % Vertical restriction alfa
LMI6 = [AA2*P + P*AA2' + BBu2*Z + Z'*BBu2' - 2*beta*P >= 0]; % Vertical restriction beta
% WT3
LMI7 = [[sin(theta)*(AA3*P + BBu3*Z + P*AA3' + Z'*BBu3'),cos(theta)*(AA3*P + BBu3*Z - P*AA3' - Z'*BBu3');
       cos(theta)*(P*AA3' + Z'*BBu3' - AA3*P - BBu3*Z),sin(theta)*(AA3*P + BBu3*Z + P*AA3' + Z'*BBu3')] <= 0];
LMI8 = [AA3*P + P*AA3' + BBu3*Z + Z'*BBu3' - 2*alfa*P <= 0]; % Vertical restriction alfa
LMI9 = [AA3*P + P*AA3' + BBu3*Z + Z'*BBu3' - 2*beta*P >= 0]; % Vertical restriction beta

LMI = [LMI0,LMI1,LMI2,LMI3,LMI4,LMI5,LMI6,LMI7,LMI8,LMI9]; 
optimize(LMI);
checkset(LMI)

Z = value(Z);
P = value(P);

K = Z*inv(P); % controlador
Kx = K(1:end-1);
Ki = K(end);

Amf_1 = AA1 + BBu1*K;
Amf_2 = AA2 + BBu2*K;
Amf_3 = AA3 + BBu3*K;

polos_MF_LQI_WT1 = eig(Amf_1)
polos_MF_LQI_WT2 = eig(Amf_2)
polos_MF_LQI_WT3 = eig(Amf_3)

figure
hold on
scatter(real(polos_MF_LQI_WT1),imag(polos_MF_LQI_WT1),'rx','LineWidth',2)
scatter(real(polos_MA_WT1),imag(polos_MA_WT1),'bx','LineWidth',2)
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
rlocus(tf(1,1))
plot(t2,l1,'k--','LineWidth', 2)
plot(t2,l2,'k--','LineWidth',2)
plot(h1,t1,'k--','LineWidth',2)
plot(h2,t1,'k--','LineWidth',2)
legend('Polos em MF_LQI','Polos em MA','Regiao D-estabilidade')
axis([-45 5 -80 80])
grid on
title('Polos WT1')

figure
hold on
scatter(real(polos_MF_LQI_WT2),imag(polos_MF_LQI_WT2),'rx','LineWidth',2)
scatter(real(polos_MA_WT2),imag(polos_MA_WT2),'bx','LineWidth',2)
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
rlocus(tf(1,1))
plot(t2,l1,'k--','LineWidth', 2)
plot(t2,l2,'k--','LineWidth',2)
plot(h1,t1,'k--','LineWidth',2)
plot(h2,t1,'k--','LineWidth',2)
legend('Polos em MF_LQI','Polos em MA','Regiao D-estabilidade')
axis([-45 5 -80 80])
grid on
title('Polos WT2')

figure
hold on
scatter(real(polos_MF_LQI_WT3),imag(polos_MF_LQI_WT3),'rx','LineWidth',2)
scatter(real(polos_MA_WT3),imag(polos_MA_WT3),'bx','LineWidth',2)
plot(pgon,'FaceColor','green','FaceAlpha',0.1)
rlocus(tf(1,1))
plot(t2,l1,'k--','LineWidth', 2)
plot(t2,l2,'k--','LineWidth',2)
plot(h1,t1,'k--','LineWidth',2)
plot(h2,t1,'k--','LineWidth',2)
legend('Polos em MF_LQI','Polos em MA','Regiao D-estabilidade')
axis([-45 5 -80 80])
grid on
title('Polos WT3')

%% Sistema com modelagem de perturbaçao
% \dot{x} = A x + B_w w
% z = C_z x + D_{zw} w

% WT1
Bw1 = [0.5; 0; 0; 0]
Cz1 = [0 1 0 0]
Dzu1 = [0.2]
Dzw1 = [0.1]
% WT2
Bw2 = 1.5*[0.5; 0; 0; 0]
Cz2 = 1.5*[0 1 0 0]
Dzu2 = 5e-2*[0.2]
Dzw2 = [0.1]
% WT3
Bw3 = 0.5*[0.5; 0; 0; 0]
Cz3 = 0.5*[0 1 0 0]
Dzu3 = 1e-2*[0.2]
Dzw3 = [0.1]

% Augmented system
AA1 = [A1 zeros(length(A1),1); -Cy1 zeros(1,1)];
AA2 = [A2 zeros(length(A2),1); -Cy2 zeros(1,1)];
AA3 = [A3 zeros(length(A3),1); -Cy3 zeros(1,1)];
BBu1 = [Bu1; -zeros(size(Cy1,1), size(Bu1,2))];
BBu2 = [Bu2; -zeros(size(Cy2,1), size(Bu2,2))];
BBu3 = [Bu3; -zeros(size(Cy3,1), size(Bu3,2))];
BBw1 = [Bu1; -zeros(size(Cy1,1), size(Bw1,2))];
BBw2 = [Bu2; -zeros(size(Cy2,1), size(Bw2,2))];
BBw3 = [Bu3; -zeros(size(Cy3,1), size(Bw3,2))];
CCz1 = [Cz1, zeros(size(Cz1,1), size(Cy1,1))];
CCz2 = [Cz2, zeros(size(Cz2,1), size(Cy2,1))];
CCz3 = [Cz3, zeros(size(Cz3,1), size(Cy3,1))];

%% H2 control + state feedback + integral action + D-stability + polytopic uncertain
% Parâmetros de Desempenho
t_min = 0.1 % min settling time
t_max = 2 % max settling time
UP = 0.20 % overshoot

alfa = -4/t_max % Vertical Strip
beta = -4/t_min % Vertical Strip
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
theta = acos(zeta) % Conic Sector angle

% Decision variables
J = sdpvar(size(BBw1,2),size(BBw1,2));
N = sdpvar(size(AA1,1), size(AA1,1));
M = sdpvar(size(BBu1,2), size(AA1,1));

LMI_integ = [N >= 0;
            [N*AA1' + M'*BBu1' + AA1*N + BBu1*M] <= 0;
            [N*AA2' + M'*BBu2' + AA2*N + BBu2*M] <= 0;
            [N*AA3' + M'*BBu3' + AA3*N + BBu3*M] <= 0];

LMI_H2 = [J >= 0;
          N >= 0;
         [J BBw1';
          BBw1 N] >= 0;
         [N*AA1' + M'*BBu1' + AA1*N + BBu1*M, N*CCz1' + M'*Dzu1';
          CCz1*N + Dzu1*M, -eye(size(CCz1,1))] <= 0;
         [J BBw2';
          BBw2 N] >= 0;
         [N*AA2' + M'*BBu2' + AA2*N + BBu2*M, N*CCz2' + M'*Dzu2';
          CCz2*N + Dzu2*M, -eye(size(CCz2,1))] <= 0;
         [J BBw3';
          BBw3 N] >= 0;
         [N*AA3' + M'*BBu3' + AA3*N + BBu3*M, N*CCz3' + M'*Dzu3';
          CCz3*N + Dzu3*M, -eye(size(CCz3,1))] <= 0];

LMI_Dstab = [N >= 0;
            [sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1'),cos(theta)*(AA1*N + BBu1*M - N*AA1' - M'*BBu1');
             cos(theta)*(N*AA1' + M'*BBu1' - AA1*N - BBu1*M),sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1')] <= 0;
            [AA1*N + N*AA1' + BBu1*M + M'*BBu1' - 2*alfa*N] <= 0;
            [-AA1*N - N*AA1' - BBu1*M - M'*BBu1' + 2*beta*N] <= 0;
            [sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2'),cos(theta)*(AA2*N + BBu2*M - N*AA2' - M'*BBu2');
             cos(theta)*(N*AA2' + M'*BBu2' - AA2*N - BBu2*M),sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2')] <= 0;
            [AA2*N + N*AA2' + BBu2*M + M'*BBu2' - 2*alfa*N] <= 0;
            [-AA2*N - N*AA2' - BBu2*M - M'*BBu2' + 2*beta*N] <= 0;
            [sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3'),cos(theta)*(AA3*N + BBu3*M - N*AA3' - M'*BBu3');
             cos(theta)*(N*AA3' + M'*BBu3' - AA3*N - BBu3*M),sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3')] <= 0;
            [AA3*N + N*AA3' + BBu3*M + M'*BBu3' - 2*alfa*N] <= 0;
            [-AA3*N - N*AA3' - BBu3*M - M'*BBu3' + 2*beta*N] <= 0];

LMI = [LMI_integ, LMI_H2, LMI_Dstab];

optimize(LMI,trace(J));
checkset(LMI);

K = value(M)*inv(value(N));
Kx = K(1:length(A1));
Ki = K(length(A1)+1:end);

%% H2 plots

w = 1:0.01:10^3;

% Tzw - Open Loop
Tzw1 = ss(AA1,BBw1,CCz1,Dzw1);
Tzw2 = ss(AA2,BBw2,CCz2,Dzw2);
Tzw3 = ss(AA3,BBw3,CCz3,Dzw3);
[mag1,~,~] = bode(Tzw1,w);
[mag2,~,~] = bode(Tzw2,w);
[mag3,~,~] = bode(Tzw3,w);
Mag1=squeeze(20*log10(mag1(1,1,:)));
Mag2=squeeze(20*log10(mag2(1,1,:)));
Mag3=squeeze(20*log10(mag3(1,1,:)));


% Tzw - Closed Loop
Tzw1_cl = ss(AA1+BBu1*K,BBw1,CCz1+Dzu1*K,Dzw1);
Tzw2_cl = ss(AA2+BBu2*K,BBw2,CCz2+Dzu2*K,Dzw3);
Tzw3_cl = ss(AA3+BBu3*K,BBw3,CCz3+Dzu3*K,Dzw3);
[mag1,~,wout1] = bode(Tzw1_cl,w);
[mag2,~,wout2] = bode(Tzw2_cl,w);
[mag3,~,wout3] = bode(Tzw3_cl,w);
Mag1_cl=squeeze(20*log10(mag1(1,1,:)));
Mag2_cl=squeeze(20*log10(mag2(1,1,:)));
Mag3_cl=squeeze(20*log10(mag3(1,1,:)));

figure
subplot(1,3,1)
semilogx(wout1,Mag1,wout1,Mag1_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT1",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,2)
semilogx(wout2,Mag2,wout2,Mag2_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT2",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,3)
semilogx(wout3,Mag3,wout3,Mag3_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT3",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
sgtitle("H2 control",'Interpreter','latex')

%% Hinf control + state feedback + integral action + D-stability + polytopic uncertain
% Parâmetros de Desempenho
t_min = 0.1 % min settling time
t_max = 2 % max settling time
UP = 0.20 % overshoot

alfa = -4/t_max % Vertical Strip
beta = -4/t_min % Vertical Strip
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
theta = acos(zeta) % Conic Sector angle

% Decision variables
N = sdpvar(size(AA1,1), size(AA1,1));
M = sdpvar(size(BBu1,2), size(AA1,1));
delta = sdpvar(1,1);

LMI_integ = [N >= 0;
            [N*AA1' + M'*BBu1' + AA1*N + BBu1*M] <= 0;
            [N*AA2' + M'*BBu2' + AA2*N + BBu2*M] <= 0;
            [N*AA3' + M'*BBu3' + AA3*N + BBu3*M] <= 0];

LMI_Hinf = [N >= 0;
           [N*AA1' + M'*BBu1' + AA1*N + BBu1*M, BBw1, N*CCz1' + M'*Dzu1';
            BBw1', -delta*eye(size(BBw1,2),size(BBw1,2)), Dzw1';
            CCz1*N + Dzu1*M, Dzw1, -eye(size(Dzw1,1),size(Dzw1,1))] <= 0;
           [N*AA2' + M'*BBu2' + AA2*N + BBu2*M, BBw2, N*CCz2' + M'*Dzu2';
            BBw2', -delta*eye(size(BBw2,2),size(BBw2,2)), Dzw2';
            CCz2*N + Dzu2*M, Dzw2, -eye(size(Dzw2,1),size(Dzw2,1))] <= 0;
           [N*AA3' + M'*BBu3' + AA3*N + BBu3*M, BBw3, N*CCz3' + M'*Dzu3';
            BBw3', -delta*eye(size(BBw3,2),size(BBw3,2)), Dzw3';
            CCz3*N + Dzu3*M, Dzw3, -eye(size(Dzw3,1),size(Dzw3,1))] <= 0];

LMI_Dstab = [N >= 0;
            [sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1'),cos(theta)*(AA1*N + BBu1*M - N*AA1' - M'*BBu1');
             cos(theta)*(N*AA1' + M'*BBu1' - AA1*N - BBu1*M),sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1')] <= 0;
            [AA1*N + N*AA1' + BBu1*M + M'*BBu1' - 2*alfa*N] <= 0;
            [-AA1*N - N*AA1' - BBu1*M - M'*BBu1' + 2*beta*N] <= 0;
            [sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2'),cos(theta)*(AA2*N + BBu2*M - N*AA2' - M'*BBu2');
             cos(theta)*(N*AA2' + M'*BBu2' - AA2*N - BBu2*M),sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2')] <= 0;
            [AA2*N + N*AA2' + BBu2*M + M'*BBu2' - 2*alfa*N] <= 0;
            [-AA2*N - N*AA2' - BBu2*M - M'*BBu2' + 2*beta*N] <= 0;
            [sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3'),cos(theta)*(AA3*N + BBu3*M - N*AA3' - M'*BBu3');
             cos(theta)*(N*AA3' + M'*BBu3' - AA3*N - BBu3*M),sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3')] <= 0;
            [AA3*N + N*AA3' + BBu3*M + M'*BBu3' - 2*alfa*N] <= 0;
            [-AA3*N - N*AA3' - BBu3*M - M'*BBu3' + 2*beta*N] <= 0];

LMI = [LMI_integ, LMI_Hinf, LMI_Dstab];

optimize(LMI,delta);
checkset(LMI);

K = value(M)*inv(value(N));
Kx = K(1:length(A1));
Ki = K(length(A1)+1:end);

%% Hinf plots

w = 1:0.01:10^3;

% Tzw - Open Loop
Tzw1 = ss(AA1,BBw1,CCz1,Dzw1);
Tzw2 = ss(AA2,BBw2,CCz2,Dzw2);
Tzw3 = ss(AA3,BBw3,CCz3,Dzw3);
[mag1,~,~] = bode(Tzw1,w);
[mag2,~,~] = bode(Tzw2,w);
[mag3,~,~] = bode(Tzw3,w);
Mag1=squeeze(20*log10(mag1(1,1,:)));
Mag2=squeeze(20*log10(mag2(1,1,:)));
Mag3=squeeze(20*log10(mag3(1,1,:)));


% Tzw - Closed Loop
Tzw1_cl = ss(AA1+BBu1*K,BBw1,CCz1+Dzu1*K,Dzw1);
Tzw2_cl = ss(AA2+BBu2*K,BBw2,CCz2+Dzu2*K,Dzw3);
Tzw3_cl = ss(AA3+BBu3*K,BBw3,CCz3+Dzu3*K,Dzw3);
[mag1,~,wout1] = bode(Tzw1_cl,w);
[mag2,~,wout2] = bode(Tzw2_cl,w);
[mag3,~,wout3] = bode(Tzw3_cl,w);
Mag1_cl=squeeze(20*log10(mag1(1,1,:)));
Mag2_cl=squeeze(20*log10(mag2(1,1,:)));
Mag3_cl=squeeze(20*log10(mag3(1,1,:)));

figure
subplot(1,3,1)
semilogx(wout1,Mag1,wout1,Mag1_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT1",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,2)
semilogx(wout2,Mag2,wout2,Mag2_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT2",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,3)
semilogx(wout3,Mag3,wout3,Mag3_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT3",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
sgtitle("H∞ control",'Interpreter','latex')

%% H2/Hinf control + state feedback + integral action + D-stability + polytopic uncertain
% Parâmetros de Desempenho
t_min = 0.1 % min settling time
t_max = 2 % max settling time
UP = 0.20 % overshoot

alfa = -4/t_max % Vertical Strip
beta = -4/t_min % Vertical Strip
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
theta = acos(zeta) % Conic Sector angle

% Decision variables
J = sdpvar(size(BBw1,2),size(BBw1,2));
N = sdpvar(size(AA1,1), size(AA1,1));
M = sdpvar(size(BBu1,2), size(AA1,1));
delta = sdpvar(1,1);

LMI_integ = [N >= 0;
            [N*AA1' + M'*BBu1' + AA1*N + BBu1*M] <= 0;
            [N*AA2' + M'*BBu2' + AA2*N + BBu2*M] <= 0;
            [N*AA3' + M'*BBu3' + AA3*N + BBu3*M] <= 0];

LMI_H2 = [J >= 0;
          N >= 0;
         [J BBw1';
          BBw1 N] >= 0;
         [N*AA1' + M'*BBu1' + AA1*N + BBu1*M, N*CCz1' + M'*Dzu1';
          CCz1*N + Dzu1*M, -eye(size(CCz1,1))] <= 0;
         [J BBw2';
          BBw2 N] >= 0;
         [N*AA2' + M'*BBu2' + AA2*N + BBu2*M, N*CCz2' + M'*Dzu2';
          CCz2*N + Dzu2*M, -eye(size(CCz2,1))] <= 0;
         [J BBw3';
          BBw3 N] >= 0;
         [N*AA3' + M'*BBu3' + AA3*N + BBu3*M, N*CCz3' + M'*Dzu3';
          CCz3*N + Dzu3*M, -eye(size(CCz3,1))] <= 0];

LMI_Hinf = [N >= 0;
           [N*AA1' + M'*BBu1' + AA1*N + BBu1*M, BBw1, N*CCz1' + M'*Dzu1';
            BBw1', -delta*eye(size(BBw1,2),size(BBw1,2)), Dzw1';
            CCz1*N + Dzu1*M, Dzw1, -eye(size(Dzw1,1),size(Dzw1,1))] <= 0;
           [N*AA2' + M'*BBu2' + AA2*N + BBu2*M, BBw2, N*CCz2' + M'*Dzu2';
            BBw2', -delta*eye(size(BBw2,2),size(BBw2,2)), Dzw2';
            CCz2*N + Dzu2*M, Dzw2, -eye(size(Dzw2,1),size(Dzw2,1))] <= 0;
           [N*AA3' + M'*BBu3' + AA3*N + BBu3*M, BBw3, N*CCz3' + M'*Dzu3';
            BBw3', -delta*eye(size(BBw3,2),size(BBw3,2)), Dzw3';
            CCz3*N + Dzu3*M, Dzw3, -eye(size(Dzw3,1),size(Dzw3,1))] <= 0];

LMI_Dstab = [N >= 0;
            [sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1'),cos(theta)*(AA1*N + BBu1*M - N*AA1' - M'*BBu1');
             cos(theta)*(N*AA1' + M'*BBu1' - AA1*N - BBu1*M),sin(theta)*(AA1*N + BBu1*M + N*AA1' + M'*BBu1')] <= 0;
            [AA1*N + N*AA1' + BBu1*M + M'*BBu1' - 2*alfa*N] <= 0;
            [-AA1*N - N*AA1' - BBu1*M - M'*BBu1' + 2*beta*N] <= 0;
            [sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2'),cos(theta)*(AA2*N + BBu2*M - N*AA2' - M'*BBu2');
             cos(theta)*(N*AA2' + M'*BBu2' - AA2*N - BBu2*M),sin(theta)*(AA2*N + BBu2*M + N*AA2' + M'*BBu2')] <= 0;
            [AA2*N + N*AA2' + BBu2*M + M'*BBu2' - 2*alfa*N] <= 0;
            [-AA2*N - N*AA2' - BBu2*M - M'*BBu2' + 2*beta*N] <= 0;
            [sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3'),cos(theta)*(AA3*N + BBu3*M - N*AA3' - M'*BBu3');
             cos(theta)*(N*AA3' + M'*BBu3' - AA3*N - BBu3*M),sin(theta)*(AA3*N + BBu3*M + N*AA3' + M'*BBu3')] <= 0;
            [AA3*N + N*AA3' + BBu3*M + M'*BBu3' - 2*alfa*N] <= 0;
            [-AA3*N - N*AA3' - BBu3*M - M'*BBu3' + 2*beta*N] <= 0];

LMI = [LMI_integ, LMI_H2, LMI_Hinf, LMI_Dstab];

optimize(LMI,trace(J));
checkset(LMI);

K = value(M)*inv(value(N));
Kx = K(1:length(A1));
Ki = K(length(A1)+1:end);

%% H2Hinf plots

w = 1:0.01:10^3;

% Tzw - Open Loop
Tzw1 = ss(AA1,BBw1,CCz1,Dzw1);
Tzw2 = ss(AA2,BBw2,CCz2,Dzw2);
Tzw3 = ss(AA3,BBw3,CCz3,Dzw3);
[mag1,~,~] = bode(Tzw1,w);
[mag2,~,~] = bode(Tzw2,w);
[mag3,~,~] = bode(Tzw3,w);
Mag1=squeeze(20*log10(mag1(1,1,:)));
Mag2=squeeze(20*log10(mag2(1,1,:)));
Mag3=squeeze(20*log10(mag3(1,1,:)));


% Tzw - Closed Loop
Tzw1_cl = ss(AA1+BBu1*K,BBw1,CCz1+Dzu1*K,Dzw1);
Tzw2_cl = ss(AA2+BBu2*K,BBw2,CCz2+Dzu2*K,Dzw3);
Tzw3_cl = ss(AA3+BBu3*K,BBw3,CCz3+Dzu3*K,Dzw3);
[mag1,~,wout1] = bode(Tzw1_cl,w);
[mag2,~,wout2] = bode(Tzw2_cl,w);
[mag3,~,wout3] = bode(Tzw3_cl,w);
Mag1_cl=squeeze(20*log10(mag1(1,1,:)));
Mag2_cl=squeeze(20*log10(mag2(1,1,:)));
Mag3_cl=squeeze(20*log10(mag3(1,1,:)));

figure
subplot(1,3,1)
semilogx(wout1,Mag1,wout1,Mag1_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT1",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,2)
semilogx(wout2,Mag2,wout2,Mag2_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT2",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
subplot(1,3,3)
semilogx(wout3,Mag3,wout3,Mag3_cl,'LineWidth',2);grid on;
axis tight
ylabel("Magnitude(dB)")
title("WT3",'Interpreter','latex')
legend("Open Loop","Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
sgtitle("H2/H∞ control",'Interpreter','latex')

% fim