clear all; close all; clc;

% Opción para silenciar los mensajes de fsolve
opts = optimset('Display', 'off');

% --- Parámetros fijos ---
K    = 0.056;     % Constante del MOSFET [A/V^2]
Vt   = 2;         % Tensión de umbral [V]
Vdd  = 12;        % Tensión de alimentación [V]

Rg1a = 2e6;       % 2 MΩ
Rg2a = 1e6;       % 1 MΩ
Rd1  = 2200;      % 2.2 kΩ

Rg1b = 2e6;       % 2 MΩ
Rg2b = 1e6;       % 1 MΩ
Rd2  = 1500;      % 1.5 kΩ
Rl   = 10e3;      % Resistencia de carga (10 kΩ)

Rref = 300;
Rs1  = 270;
Rs2  = 100;

Vg1  = (Vdd * Rg2a) / (Rg1a + Rg2a); % 4 V
Rg1  = (Rg1a * Rg2a) / (Rg1a + Rg2a);
Vg2  = (Vdd * Rg2b) / (Rg1b + Rg2b); % 4 V

Vi   = 1;         % Tensión de prueba (1 V)

% =========================================================================
% POLARIZACIÓN Y PARÁMETROS DC
% =========================================================================

% 1. Circuito de Referencia (Espejo)
sistemap = @(x) [
    x(1) - K*(x(2) - Vt)^2;
    Vdd - x(1)*Rref - x(2)
];
[xp, ~, info_p] = fsolve(sistemap, [0.001, 2.5], opts);

Vp   = xp(2);
gmp  = 2 * K * (Vp - Vt);
rdsp = (Vp) / xp(1); % Tensión Vds del MOSFET p
mup  = gmp * rdsp;

% 2. Etapa 1
sistema1 = @(x) [
    x(1) - K * (x(2) - Vt)^2;
    x(1) - K * (x(3) - Vt)^2;
    x(2) - (Vg1 - (x(5) + x(1)*Rs1));
    x(3) - (Vp - x(1)*Rs1);
    Vdd - x(1)*(Rd1 + Rs1) - x(4) - x(5)
];
[x1, ~, info_1] = fsolve(sistema1, [0.001, 2.5, 2.5, 3, 3], opts);

gm1   = 2 * K * (x1(2) - Vt);
gmo1  = 2 * K * (x1(3) - Vt);
rds1  = x1(4) / x1(1);
rdso1 = x1(5) / x1(1);
mu1   = gm1 * rds1;
muo1  = gmo1 * rdso1;

% 3. Etapa 2
sistema2 = @(x) [
    x(1) - K * (x(2) - Vt)^2;
    x(1) - K * (x(3) - Vt)^2;
    x(2) - (Vg2 - (x(5) + x(1)*Rs2));
    x(3) - (Vp - x(1)*Rs2);
    Vdd - x(1)*(Rd2 + Rs2) - x(4) - x(5)
];
[x2, ~, info_2] = fsolve(sistema2, [0.001, 2.5, 2.5, 3, 3], opts);

gm2   = 2 * K * (x2(2) - Vt);
gmo2  = 2 * K * (x2(3) - Vt);
rds2  = x2(4) / x2(1);
rdso2 = x2(5) / x2(1);
mu2   = gm2 * rds2;
muo2  = gmo2 * rdso2;

% =========================================================================
% CÁLCULO DE GANANCIA (PEQUEÑA SEÑAL)
% =========================================================================
A = zeros(14, 14); b = zeros(14, 1);
A(1,1) = 1/Rg1; A(1,10) = 1;
A(2,2) = 1/rds1; A(2,3) = -1/rds1; A(2,11) = 1;
A(3,2) = -1/rds1; A(3,3) = (1/rds1)+(1/Rd1)+(1/rdso2); A(3,6) = -1/rdso2; A(3,13) = 1;
A(4,4) = 1/rds2; A(4,5) = -1/rds2; A(4,13) = -1;
A(5,4) = -1/rds2; A(5,5) = (1/rds2)+(1/Rd2)+(1/Rl);
A(6,3) = -1/rdso2; A(6,6) = 1/rdso2; A(6,12) = -1;
A(7,7) = 1/Rs2; A(7,12) = 1;
A(8,8) = (1/Rref)+(1/rdsp); A(8,9) = -1/rdsp;
A(9,8) = -1/rdsp; A(9,9) = 1/rdsp; A(9,14) = 1;

A(10,1) = 1; b(10) = Vi;
A(11,2) = 1; A(11,1) = -mu1;
A(12,6) = -1; A(12,7) = 1 + muo2; A(12,8) = -muo2;
A(13,3) = 1 + mu2; A(13,4) = -1;
A(14,8) = -mup; A(14,9) = 1;

X = A \ b;
Ganancia = abs(X(5) / Vi);

% =========================================================================
% CÁLCULO DE POTENCIAS EN POLARIZACIÓN (DC)
% =========================================================================
% Corrientes y tensiones extraídas
Idsp = xp(1); Vdsp = Vdd - Vp;
Ids1 = x1(1); Vds1_1 = x1(4); Vds1_2 = x1(5);
Ids2 = x2(1); Vds2_1 = x2(4); Vds2_2 = x2(5);

% Potencias en Resistencias [mW]
P_Rref = (Idsp^2 * Rref) * 1000;
P_Rs1  = (Ids1^2 * Rs1) * 1000;
P_Rd1  = (Ids1^2 * Rd1) * 1000;
P_Rs2  = (Ids2^2 * Rs2) * 1000;
P_Rd2  = (Ids2^2 * Rd2) * 1000;

P_Rg1a = ((Vdd - Vg1)^2 / Rg1a) * 1000;
P_Rg2a = (Vg1^2 / Rg2a) * 1000;
P_Rg1b = ((Vdd - Vg2)^2 / Rg1b) * 1000;
P_Rg2b = (Vg2^2 / Rg2b) * 1000;

% Potencias en Transistores [mW] (P = Vds * Ids)
P_Mp   = (Vdsp * Idsp) * 1000;
P_M1_1 = (Vds1_1 * Ids1) * 1000;
P_M1_2 = (Vds1_2 * Ids1) * 1000;
P_M2_1 = (Vds2_1 * Ids2) * 1000;
P_M2_2 = (Vds2_2 * Ids2) * 1000;

P_total_res = P_Rref + P_Rs1 + P_Rd1 + P_Rs2 + P_Rd2 + P_Rg1a + P_Rg2a + P_Rg1b + P_Rg2b;
P_total_mos = P_Mp + P_M1_1 + P_M1_2 + P_M2_1 + P_M2_2;

fprintf('\n================== RESULTADOS DE GANANCIA ==================\n');
fprintf('Ganancia |Av|: %.4f (%.2f dB)\n', Ganancia, 20*log10(Ganancia));

% =========================================================================
% IMPRESIÓN DE PARÁMETROS PEQUEÑA SEÑAL (rds, gm, mu)
% =========================================================================
fprintf('\n================ PARÁMETROS PEQUEÑA SEÑAL (rds, gm, mu) ================\n');
fprintf('%-25s | %-12s | %-12s | %-10s\n', 'Transistor', 'rds [Ω]', 'gm [mS]', 'mu');
fprintf('------------------------------------------------------------------------\n');
fprintf('%-25s | %-12.2f | %-12.4f | %-10.2f\n', 'Referencia (M_p)',          rdsp, gmp*1000, mup);
fprintf('%-25s | %-12.2f | %-12.4f | %-10.2f\n', 'Etapa 1 - Superior (M1_1)', rds1, gm1*1000, mu1);
fprintf('%-25s | %-12.2f | %-12.4f | %-10.2f\n', 'Etapa 1 - Inferior (M1_2)', rdso1, gmo1*1000, muo1);
fprintf('%-25s | %-12.2f | %-12.4f | %-10.2f\n', 'Etapa 2 - Superior (M2_1)', rds2, gm2*1000, mu2);
fprintf('%-25s | %-12.2f | %-12.4f | %-10.2f\n', 'Etapa 2 - Inferior (M2_2)', rdso2, gmo2*1000, muo2);
fprintf('========================================================================\n');

fprintf('\n================ DISIPACIÓN DE POTENCIA [mW] ================\n');
fprintf('1. RESISTENCIAS:\n');
fprintf('   P_Rref = %.4f mW\n', P_Rref);
fprintf('   P_Rs1  = %.4f mW | P_Rd1 = %.4f mW\n', P_Rs1, P_Rd1);
fprintf('   P_Rs2  = %.4f mW | P_Rd2 = %.4f mW\n', P_Rs2, P_Rd2);
fprintf('   Divisor A: P_Rg1a = %.4f mW | P_Rg2a = %.4f mW\n', P_Rg1a, P_Rg2a);
fprintf('   Divisor B: P_Rg1b = %.4f mW | P_Rg2b = %.4f mW\n', P_Rg1b, P_Rg2b);
fprintf('   Total en Resistencias: %.4f mW\n', P_total_res);

fprintf('\n2. TRANSISTORES (P = Vds * Ids):\n');
fprintf('   Referencia (M_p) : %.4f mW  (Vds = %.2f V, Ids = %.2f mA)\n', P_Mp, Vdsp, Idsp*1000);
fprintf('   Etapa 1 - M1_1   : %.4f mW  (Vds = %.2f V, Ids = %.2f mA)\n', P_M1_1, Vds1_1, Ids1*1000);
fprintf('   Etapa 1 - M1_2   : %.4f mW  (Vds = %.2f V, Ids = %.2f mA)\n', P_M1_2, Vds1_2, Ids1*1000);
fprintf('   Etapa 2 - M2_1   : %.4f mW  (Vds = %.2f V, Ids = %.2f mA)\n', P_M2_1, Vds2_1, Ids2*1000);
fprintf('   Etapa 2 - M2_2   : %.4f mW  (Vds = %.2f V, Ids = %.2f mA)\n', P_M2_2, Vds2_2, Ids2*1000);
fprintf('   Total en Transistores: %.4f mW\n', P_total_mos);

fprintf('\nPOTENCIA TOTAL CONSUMIDA DE VDD: %.4f mW\n', P_total_res + P_total_mos);
fprintf('============================================================\n');

% =========================================================================
% GRAFICACIÓN DE INTERSECCIÓN (CARACTERÍSTICA vs RECTA DE POLARIZACIÓN)
% =========================================================================
vgs_axis = linspace(Vt, Vdd/2, 500);
ids_parabola = K * (vgs_axis - Vt).^2;

% --- Grafico 1: Referencia ---
figure(1,'Name', 'Puntos de Operacion Q Referencia', 'NumberTitle', 'off');
ids_recta_ref = (Vdd - vgs_axis) / Rref;
plot(vgs_axis, ids_parabola*1000, 'b-', 'LineWidth', 2); hold on;
plot(vgs_axis, ids_recta_ref*1000, 'r--', 'LineWidth', 1.5);
plot(xp(2), xp(1)*1000, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
ylim([0 200]);
grid on;
title('Referencia (MOSFET p)');
xlabel('Vgs [V]'); ylabel('Ids [mA]');
legend('Parábola MOSFET', 'Recta Carga', 'Punto Q', 'Location', 'southeast');

% --- Grafico 2: Etapa 1 ---
figure(2,'Name', 'Puntos de Operacion Q Etapa 1', 'NumberTitle', 'off');
ids_recta_e1 = (Vp - vgs_axis) / Rs1;
plot(vgs_axis, ids_parabola*1000, 'b-', 'LineWidth', 2); hold on;
plot(vgs_axis, ids_recta_e1*1000, 'r--', 'LineWidth', 1.5);
plot(x1(2), x1(1)*1000, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
ylim([0 10]);
xlim([0 4]);
grid on;
title('Etapa 1 (Transistor Superior)');
xlabel('Vgs [V]'); ylabel('Ids [mA]');
legend('Parábola MOSFET', 'Recta de Polarización', 'Punto Q', 'Location', 'northeast');

% --- Grafico 3: Etapa 2 ---
figure(3,'Name', 'Puntos de Operacion Q Etapa 2', 'NumberTitle', 'off');
ids_recta_e2 = (Vp - vgs_axis) / Rs2;
plot(vgs_axis, ids_parabola*1000, 'b-', 'LineWidth', 2); hold on;
plot(vgs_axis, ids_recta_e2*1000, 'r--', 'LineWidth', 1.5);
plot(x2(2), x2(1)*1000, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
ylim([0 10]);
xlim([0 4]);
grid on;
title('Etapa 2 (Transistor Superior)');
xlabel('Vgs [V]'); ylabel('Ids [mA]');
legend('Parábola MOSFET', 'Recta de Polarización', 'Punto Q', 'Location', 'northeast');

% =========================================================================
%   RECTAS DE CARGA EN CONTINUA (DC) Y EN ALTERNA (AC)
% =========================================================================

fprintf('\n Impedancias de cada etapa \n');

Rl = 1000;

% Primer Etapa

Raux1 = rds1+(rdso1+Rs1*(muo1+1))*(mu1+1);
Zout1 = 1/(1/Rd1+1/Raux1)

Zin1 = (Rg1a * Rg2a) / (Rg1a + Rg2a)

% Segunda Etapa

Raux1 = rds2+rdso2*(mu2+1)+Rs2*(mu2+1)*(muo2+1);
Zout2 = 1/(1/Rd2+1/Raux1)

Raux2 = rdso2+Rs2/(muo2+1);
Raux3 = (rds2+1/(1/Rd2+1/Rl))*(mu2+1);
Zin2 = 1/(1/Raux2+1/Raux3)

% Resistencias de carga DC para cada etapa
Rdc1 = Rd1 + Rs1;
Rdc2 = Rd2 + Rs2;

% Resistencias de carga AC efectivas vistas en el drenador de cada etapa
% Etapa 1 ve hacia arriba la resistencia Rd1 en paralelo con la impedancia de la etapa siguiente
Rac1 = (Rd1 * Zin2) / (Rd1 + Zin2);
% Etapa 2 (salida) ve Rd2 en paralelo con Rl
Rac2 = (Rd2 * Rl) / (Rd2 + Rl);

vds_vec = linspace(0, Vdd, 500);

% --- Recta de Carga DC y AC para Transistor Inferior Etapa 1 (M1_1) ---
figure(4, 'Name', 'Rectas de Carga - Etapa 1 (M1_1)', 'NumberTitle', 'off');
% Ecuación DC: Ids = Ids1 + (Vds1_1 - vds_vec) / Rdc1
ids_dc_e1 = Ids1 + (Vds1_1 - vds_vec) / Rdc1;
% Ecuación AC: Ids = Ids1 + (Vds1_1 - vds_vec) / Rac1
ids_ac_e1 = Ids1 + (Vds1_1 - vds_vec) / Rac1;

plot(vds_vec, ids_dc_e1 * 1000, 'r-', 'LineWidth', 2); hold on;
plot(vds_vec, ids_ac_e1 * 1000, 'b--', 'LineWidth', 2);
plot(Vds1_1, Ids1 * 1000, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
grid on;
title('Rectas de Carga DC y AC - Etapa 1 (M1_1)');
xlabel('Vds [V]'); ylabel('Ids [mA]');
legend('Recta de Carga DC', 'Recta de Carga AC', 'Punto Q', 'Location', 'northeast');
ylim([0 max(ids_dc_e1)*1.2*1000]); xlim([0 Vdd]);

% --- Recta de Carga DC y AC para Transistor Inferior Etapa 2 (M2_1) ---
figure(5, 'Name', 'Rectas de Carga - Etapa 2 (M2_1)', 'NumberTitle', 'off');
% Ecuación DC: Ids = Ids2 + (Vds2_1 - vds_vec) / Rdc2
ids_dc_e2 = Ids2 + (Vds2_1 - vds_vec) / Rdc2;
% Ecuación AC: Ids = Ids2 + (Vds2_1 - vds_vec) / Rac2
ids_ac_e2 = Ids2 + (Vds2_1 - vds_vec) / Rac2;

plot(vds_vec, ids_dc_e2 * 1000, 'r-', 'LineWidth', 2); hold on;
plot(vds_vec, ids_ac_e2 * 1000, 'b--', 'LineWidth', 2);
plot(Vds2_1, Ids2 * 1000, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
grid on;
title('Rectas de Carga DC y AC - Etapa 2 (M2_1)');
xlabel('Vds [V]'); ylabel('Ids [mA]');
legend('Recta de Carga DC', 'Recta de Carga AC', 'Punto Q', 'Location', 'northeast');
ylim([0 max(ids_dc_e2)*1.2*1000]); xlim([0 Vdd]);
