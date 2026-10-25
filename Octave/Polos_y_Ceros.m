clear all; close all; clc

Rd1  = 2200;     % 2.2 kΩ
Rd2  = 1500;     % 1.5 kΩ
Rl   = 1000;     % Resistencia de carga Rl [Ω]
Rs1  = 270;
Rs2  = 100;

rds1 = 2793.71;
rdso1 = 616.25;
rds2 = 751.97;
rdso2 = 275.36;

gm1 = 21.4264e-3;
gmo1 = 21.4264e-3;
gm2 = 31.9865e-3;
gmo2 = 31.9865e-3;

mu1 = rds1 * gm1;
muo1 = rdso1 * gmo1;
mu2 = rds2 * gm2;
muo2 = rdso2 * gmo2;

Co = 10e-6;
Cc = 470e-6;

%   ===========================
%   ANÁLISIS EN BAJA FRECUENCIA
%   ===========================

% FRECUENCIAS DE CORTE DE LOS POLOS
% Chf está siempre abierto
% se asume Cs1 como el polo dominante

%Cálculos para evaluar Cc
Ra1 = rds1 + rdso1*(mu1+1) + Rs1*(mu1+1)*(muo1+1);
Ra2 = rdso2 + Rs2*(muo2+1);
Ra3 = (rds2+Rd2)/(mu2+1);
RCc1 = 1/(1/Rd1+1/Ra1)+1/(1/Ra2+1/Ra3) %Considero Co abierto

Ra4 = (rds2 + 1/(1/Rd2+1/Rl))/(mu2+1);
RCc2 = 1/(1/Rd1+1/Ra1)+1/(1/Ra2+1/Ra4) %Considero Co cortado

%Cálculos para evaluar Co
Ra5 = rds1+(mu1+1)*(rdso1+Rs1*(muo1+1));
Ra6 = 1/(1/Rd1+1/Ra5+1/Ra2);
RCo1 = Rl + 1/(1/Rd2+1/(Ra6*(mu2+1)+rds2)) %Cosnidero Cc cortado

Ra7 = rds2+rdso2*(mu2+1)+Rs2*(mu2+1)*(muo2+1);
RCo2 = Rl +1/(1/Rd2+1/Ra7) %Considero Cc abierto

disp('');
disp('FRECUENCIAS DE CORTE EN BAJA');
disp('Cuando Co está cortado y Cc abierto:');
fprintf('Wco %.4f  rad/s\n', 1/(Co*RCo2));
fprintf('Wcc %.4f  rad/s\n', 1/(Cc*RCc2));

disp('Cuando Cc está cortado y Co abierto:');
fprintf('Wco %.4f  rad/s\n', 1/(Co*RCo1));
fprintf('Wcc %.4f  rad/s\n', 1/(Cc*RCc1));

wcc = 1/(Cc*RCc1);
wco = 1/(Co*RCo1);

% Co está abierto cuando Cc está cortando
% Para un Cs1 dominante su w debe ser mayor o igual a 400 rad/s
disp('');
disp('Cálculo de Cs1');
Ra8 = (rds2 + (Rd2*Rl)/(Rd2+Rl))/(mu2+1);
Ra9 = rdso2+Rs2*(muo2+1);
Ra10 = 1/(1/Rd1 + 1/Ra9 + 1/Ra8);
Ra11 = rdso1+Rs1*(muo1+1);
RCs1 = 1/(1/Ra10+1/Ra11);

Wps1 = 400; %rad/s
fprintf('Cs1 = %.4f uF\n',(1/(RCs1*Wps1))*1e6);

%   ===========================
%   ANÁLISIS EN ALTA FRECUENCIA
%   ===========================

Ciss = 50e-12;
Coss = 25e-12;
Crss = 5e-12;
Cgd = Crss;

Cds = Coss - Cgd;
Cgs = Ciss - Cgd;

Req = 1/(1/rds2+1/Rd2+1/Rl);
wphf = 20e3;

Ceq = 1/(Req * wphf);

Chf = Ceq - Cgd;

wzhf = 1 / (Ceq * Rd2);
wpg2 = wzhf;

disp('');
disp('CALCULO DE Chf');
fprintf('Chf %.4f  nF\n', Chf*1e9);
fprintf('Wp %.4f  rad/s\n', wphf);
fprintf('Wz %.4f  rad/s\n', wzhf);

Rg2 = 2/3 * 1e6;
Rb1 = rdso2 + Rs2 * (muo2 +1);
Rth = 1/(1/rds1 + 1/Rd1 + 1/Rb1);

Rs = Rth * (mu2+1);
Rgs = Rs + Rg2;
Req = (Rg2 * Rs)/(Rg2+Rs);

wpgs = 1/(Cgs*Rgs)
Ceq2 = 1/(wpg2*Req);
Cg2 = Ceq - Cgd
wzg2 = 1/(Ceq2*Rg2)

%   =========================
%   FUNCIÓN DE TRANSFERENCIA
%   =========================

pkg load control
s = tf('s');

Av =24.6879;

num = s^4*(1 + s/wzhf)*(1 + s/wzg2);
den = (s + wco)*(s + wcc)*(s + wphf)*(1 + s/wpg2)*(1 + s/wphf)*(1 + s/wpgs);
fdt = Av *num/den;

figure(1);
bode(fdt);

##figure(2);
##pzmap(fdt);
##grid on;


