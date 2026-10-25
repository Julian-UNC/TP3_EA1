clear all; close all; clc

Vi   = 1.0;      % Tensión de entrada de prueba [V]
Rg1  = 666666.67;% 0.66 MΩ
%rds  = 3.2      % Encontrado ed datasheet
rds  = 7.5;    % gemini está loco
Rd1  = 2200;     % 1.5 kΩ
Rs2  = 6500;     % Resistencia de fuente Rs2 [Ω]
Rd2  = 1500;     % 1.5 kΩ
Rl   = 1200;     % Resistencia de carga Rl [Ω]
%mu1  = 0.036464;       % Factor de amplificación del MOSFET 1 (mu = gm * rds)
%mu2  = 0.45912;       % Factor de amplificación del MOSFET 2
gm1 = gm2 = 0.320; % tipico según DS
%gm1 = 0.0138;   %calculados según la polarización
%gm2 = 0.01007;
mu1 = gm1*rds
mu2 = gm2*rds


%     ARMADO DE LA MATRIZ

% Vector de incógnitas X = [V_N001; V_N003; V_N005; V_N007]
% (Nota: V_N002, V_N004 y V_N006 ya están expresados en función de estos nodos)

% Dimensionamos la matriz A y el vector B (4x4)
A = zeros(4, 4);
B = zeros(4, 1);

% --- Fila 1: Nodo N001 ---
% V(N001) = Vi
A(1, 1) = 1;
B(1, 1) = Vi;

% --- Fila 2: Nodo N003 (KCL en Nodo 3) ---
% Conectado a: R2 (hacia N002), Rd1 (a 0), R4 (hacia N006), y rama V4-R6
% Recordando sustituciones:
% V_N002 = -mu1 * V_N001
% V_N004 = (1 + mu2) * V_N003
% V_N006 = (1 + mu2) * V_N007
G_rds = 1 / rds;
G_Rd1 = 1 / Rd1;

% Coeficientes para [V_N001, V_N003, V_N005, V_N007]
A(2, 1) = G_rds * mu1;                                     % Derivado de V_N002
A(2, 2) = G_rds + G_Rd1 + G_rds + G_rds*(1 + mu2);         % Términos en V_003
A(2, 3) = -G_rds;                                          % Desde R6 (hacia N005)
A(2, 4) = -G_rds;                                          % Desde R4 (hacia N006)
B(2, 0+1) = 0;

% --- Fila 3: Nodo N005 (Nodo de Salida - KCL) ---
% Conectado a: R6 (hacia N004), Rd2 (a 0), Rl (a 0)
G_Rd2 = 1 / Rd2;
G_Rl  = 1 / Rl;

A(3, 2) = -G_rds * (1 + mu2);                             % Derivado de V_N004
A(3, 3) = G_rds + G_Rd2 + G_Rl;                           % Términos en V_N005
A(3, 1) = 0;
A(3, 4) = 0;
B(3, 1) = 0;

% --- Fila 4: Nodo N007 (KCL en Nodo 7) ---
% Conectado a: Rs2 (a 0) y la fuente V3 conectada a N006
G_Rs2 = 1 / Rs2;
% De V3: V_N007 - V_N006 = -mu2 * V_N007  =>  V_N006 = (1 + mu2) * V_N007
% La corriente por la fuente V3 hacia N006 se equilibra con la rama R4
A(4, 4) = G_Rs2 + (1 + mu2) * G_rds;
A(4, 2) = -G_rds;                                          % Desde R4
A(4, 1) = 0;
A(4, 3) = 0;
B(4, 1) = 0;

%   SOLUCIÓN DEL SISTEMA

X = A \ B;

% Extracción de tensiones nodales principales
V_N001 = X(1);
V_N003 = X(2);
V_N005 = X(3);  % <-- SALIDA DEL AMPLIFICADOR
V_N007 = X(4);

% Tensiones secundarias dependientes
V_N002 = -mu1 * V_N001;
V_N004 = (1 + mu2) * V_N003;
V_N006 = (1 + mu2) * V_N007;

% Cálculo de la Ganancia de Tensión (Av = V_out / Vi)
Ganancia_V = V_N005 / Vi;
Ganancia_dB = 20 * log10(abs(Ganancia_V));

fprintf('\n================ RESULTADOS DEL ANÁLISIS NODAL ================\n');
fprintf(' Tensión de Entrada (Vi)   : %.4f V\n', Vi);
fprintf(' Tensión de Salida (Nod 5) : %.4f V\n', V_N005);
fprintf('---------------------------------------------------------------\n');
fprintf(' Ganancia de Tensión (Av)  : %.4f V/V\n', Ganancia_V);
fprintf(' Ganancia en dB            : %.4f dB\n', Ganancia_dB);
fprintf('===============================================================\n\n');
