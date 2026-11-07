clear all; close all; clc

Rd1  = 2200;     % 2.2 kΩ
Rd2  = 1500;     % 1.5 kΩ
Rl   = 1000;     % Resistencia de carga Rl [Ω]
Rs1  = 270;
Rs2  = 100;
Rl   = 1000;

rds1  = 2793.71;
rdso1 = 616.25;
rds2  = 751.97;
rdso2 = 275.36;

gm1   = 21.4264e-3;
gmo1  = 21.4264e-3;
gm2   = 31.9865e-3;
gmo2  = 31.9865e-3;

mu1  = rds1 * gm1;
muo1 = rdso1 * gmo1;
mu2  = rds2 * gm2;
muo2 = rdso2 * gmo2;

Co = 10e-6;
Ci = 10e-6;
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

wps1 = 400; %rad/s
Cs1 = 1/(RCs1*wps1);
fprintf('Cs1 = %.4f uF\n',(1/(RCs1*wps1))*1e6);

%   ===========================
%   ANÁLISIS EN ALTA FRECUENCIA
%   ===========================

Ciss = 50e-12;
Coss = 25e-12;
Crss = 5e-12;
Cgd = Crss;

Cds = Coss - Cgd;
Cgs = Ciss - Cgd;

Req = 1/(1/Rd2+1/Rl);
wphf = 2*3.14*20e3;

Ceq = 1/(Req * wphf);

Chf = Ceq - Cgd;
Rd2l = (Rd2*Rl)/(Rd2+Rl);

wzhf = 1 / (Ceq * Rd2l);
wpg2 = wzhf;

disp('');
disp('CALCULO DE Chf');
fprintf('Chf %.4f  nF\n', Chf*1e9);
fprintf('Wp %.4f  rad/s\n', wphf);
fprintf('Wz %.4f  rad/s\n', wzhf);

Rg2 = 2/3 * 1e6;
Rb1 = rdso2 + Rs2 * (muo2 +1);
Rth = 1/(+ 1/Rd1 + 1/Rb1);

Rs = Rth * (mu2+1);
Rgs = Rs + Rg2;
Req2 = (Rg2 * Rs)/(Rg2+Rs);

wpgs = 1/(Cgs*Rgs)
Ceq2 = 1/(wpg2*Req2);
Cg2 = Ceq - Cgd
wzg2 = 1/(Ceq2*Rg2)


polosinChf = 1/(Cgd*Rth);

%   =========================
%   FUNCIÓN DE TRANSFERENCIA
%   =========================

pkg load control
s = tf('s');

Av =24.6879;

numL = s^3;
denL = (s + wco)*(s + wcc)*(s + wps1);
fdtL = Av *numL/denL;

fdtH = Av/(1+s/wphf);

fdt = Av/(1+s/wphf) * numL/denL;

fdtsinChf = Av/(1+s/polosinChf) * numL/denL;

figure(1);
bode(fdt);

figure(2);
bode(fdtsinChf);

##polo = 1/(Ceq*Req)
##1/polo
##
##figure(3);
##bode(fdtH);

%% ETAPA 1

Av1 = (-mu1*Rd1)/(rds1+Rd1);

wi = 1/(Ci*Rb1);

Rs1_eq = 1/( 1/(rdso1 + Rs1 * (muo1 + 1)) + (mu1 + 1)/(rds1 + Rd1));

ws1 = 1/(Cs1*Rs1);

Rc = Rd1 + rds1 + (mu1 + 1)*(rdso1 + Rd1*(muo1 + 1));

wc = 1/(Cc*Rc);

FdTE1 = (s^3 * Av1)/( (s + wi) * (s + ws1) * (s + wc));

figure(3);
bode(FdTE1);

%% ETAPA 2 CON CHF

Rdl = (Rd1*Rl)/(Rd1+Rl);

Av2 = (-mu2*Rdl)/(Rdl + rds2);

Rcc = 1/( 1/(rdso2 + (muo2 + 1)*Rs2) + 1/(rds2 + Rd2));

wcc_e2 = 1/(Rcc*Cc);

Ro = Rl + (Rd2*rds2)/(Rd2+rds2);

wo = 1/(Ro*Co);

FdTE2 = (Av2*s^2)/( (s+wcc) * (s+wo) * (1+s/wphf));

figure(4);
bode(FdTE2);

% Corrección en el denominador de FdTE2:
FdTE2 = (Av2 * s^2) / ((s + wcc) * (s + wo) * (1 + s/wphf));

freqs = [10, 100, 1000, 10000, 20000, 100000];
w_vec = 2 * pi * freqs;

[mag1, ~] = bode(FdTE1, w_vec);
[mag2, ~] = bode(FdTE2, w_vec);
[mag_tot, ~] = bode(fdt, w_vec);

mag1 = squeeze(mag1);
mag2 = squeeze(mag2);
mag_tot = squeeze(mag_tot);

fprintf('\n========================================================\n');
fprintf(' f (Hz)    |  |Etapa 1| (V/V) |  |Etapa 2| (V/V) |  |Total| (V/V)\n');
fprintf('--------------------------------------------------------\n');
for k = 1:length(freqs)
    fprintf('%8d   | %13.4f | %13.4f | %11.4f\n', ...
        freqs(k), mag1(k), mag2(k), mag_tot(k));
end
fprintf('========================================================\n');
