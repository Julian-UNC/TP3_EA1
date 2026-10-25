% =========================================================================
% Optimización del Amplificador por Pasos en Octave (Corregido y Estructurado)
% =========================================================================
clear; clc; close all;

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

Vg1  = (Vdd * Rg2a) / (Rg1a + Rg2a); % 4 V
Rg1  = (Rg1a * Rg2a) / (Rg1a + Rg2a);
Vg2  = (Vdd * Rg2b) / (Rg1b + Rg2b); % 4 V

Vi   = 1;         % Tensión de prueba (1 V)

% --- Rangos de búsqueda ---
vec_Rref = 100:100:5000; % Rref: 100 Ω a 5 kΩ
vec_Rs2  = 10:10:470;    % Rs2: 10 Ω a 470 Ω
vec_Rs1  = 10:10:470;    % Rs1: 10 Ω a 470 Ω

opts = optimset('Display', 'off');

% =========================================================================
% PASO 1: Iterar Rref y Rs2 para maximizar la Ganancia (|Av| > 6) y rds > 70
% =========================================================================
fprintf('PASO 1: Buscando Rref y Rs2 óptimos...\n');

max_ganancia = 0;
opt_Rref = 0;
opt_Rs2  = 0;
opt_Vp   = 0;
opt_x2   = [];
opt_xp   = [];

for Rref = vec_Rref
    % Sistema de referencia (Espejo)
    sistemap = @(x) [
        x(1) - K*(x(2) - Vt)^2;
        Vdd - x(1)*Rref - x(2)
    ];
    [xp, ~, info_p] = fsolve(sistemap, [0.001, 2.5], opts);
    if info_p <= 0 || xp(2) <= Vt, continue; end

    Vp = xp(2);
    % Vp - x(1)*Rref es aproximadamente el Vds del transistor espejo
    rdsp_val = (Vp) / xp(1);
    if rdsp_val < 70, continue; end

    for Rs2 = vec_Rs2
        % Sistema Etapa 2
        sistema2 = @(x) [
            x(1) - K * (x(2) - Vt)^2;
            x(1) - K * (x(3) - Vt)^2;
            x(2) - (Vg2 - (x(5) + x(1)*Rs2));
            x(3) - (Vp - x(1)*Rs2);
            Vdd - x(1)*(Rd2 + Rs2) - x(4) - x(5)
        ];
        [x2, ~, info_2] = fsolve(sistema2, [0.001, 2.5, 2.5, 3, 3], opts);

        if info_2 <= 0, continue; end

        % Condición de saturación para Etapa 2
        if (x2(2) <= Vt || x2(3) <= Vt || x2(4) < (x2(2)-Vt) || x2(5) < (x2(3)-Vt))
            continue;
        end

        % Cálculo de rds físicos (Vds / Ids)
        rds2_val  = x2(4) / x2(1);
        rdso2_val = x2(5) / x2(1);

        % Restricción: rds debe ser mayor a 70 ohms
        if (rds2_val < 70 || rdso2_val < 70), continue; end

        gm2  = 2 * K * (x2(2) - Vt);
        gmo2 = 2 * K * (x2(3) - Vt);
        mu2  = gm2 * rds2_val;
        muo2 = gmo2 * rdso2_val;

        % Estimación provisional de mu1 para testear ganancia global
        mu1_est = mu2;

        A = zeros(14, 14); b = zeros(14, 1);
        A(1,1) = 1/Rg1; A(1,10) = 1;
        A(2,2) = 1/rds2_val; A(2,3) = -1/rds2_val; A(2,11) = 1; % Aproximación inicial rds1=rds2
        A(3,2) = -1/rds2_val; A(3,3) = (1/rds2_val)+(1/Rd1)+(1/rdso2_val); A(3,6) = -1/rdso2_val; A(3,13) = 1;
        A(4,4) = 1/rds2_val; A(4,5) = -1/rds2_val; A(4,13) = -1;
        A(5,4) = -1/rds2_val; A(5,5) = (1/rds2_val)+(1/Rd2)+(1/Rl);
        A(6,3) = -1/rdso2_val; A(6,6) = 1/rdso2_val; A(6,12) = -1;
        A(7,7) = 1/Rs2; A(7,12) = 1;
        A(8,8) = (1/Rref)+(1/rdsp_val); A(8,9) = -1/rdsp_val;
        A(9,8) = -1/rdsp_val; A(9,9) = 1/rdsp_val; A(9,14) = 1;

        A(10,1) = 1; b(10) = Vi;
        A(11,2) = 1; A(11,1) = -mu1_est;
        A(12,6) = -1; A(12,7) = 1 + muo2; A(12,8) = -muo2;
        A(13,3) = 1 + mu2; A(13,4) = -1;
        gmp_val = 2 * K * (Vp - Vt);
        mup_val = gmp_val * rdsp_val;
        A(14,8) = -mup_val; A(14,9) = 1;

        X = A \ b;
        Ganancia = abs(X(5) / Vi);

        if (Ganancia > 6 && Ganancia > max_ganancia)
            max_ganancia = Ganancia;
            opt_Rref = Rref;
            opt_Rs2  = Rs2;
            opt_Vp   = Vp;
            opt_x2   = x2;
            opt_xp   = xp;
        end
    end
end

if opt_Rref == 0
    error('No se encontraron valores que cumplan |Av| > 6 y rds > 70 Ω.');
end

fprintf('Paso 1 Completado:\n  Rref óptimo = %.2f Ω\n  Rs2 óptimo  = %.2f Ω\n  Ganancia preliminar = %.2f\n\n', ...
        opt_Rref, opt_Rs2, max_ganancia);

% =========================================================================
% PASO 2: Iterar Rs1 para lograr Máxima Excursión Simétrica en Etapa 1
% =========================================================================
fprintf('PASO 2: Buscando Rs1 para máxima excursión simétrica (con rds > 70 Ω)...\n');

opt_Rs1 = 0;
min_error_excursion = inf;
mejores_datos1 = [];

for Rs1 = vec_Rs1
    sistema1 = @(x) [
        x(1) - K * (x(2) - Vt)^2;
        x(1) - K * (x(3) - Vt)^2;
        x(2) - (Vg1 - (x(5) + x(1)*Rs1));
        x(3) - (opt_Vp - x(1)*Rs1);
        Vdd - x(1)*(Rd1 + Rs1) - x(4) - x(5)
    ];
    [x1, ~, info_1] = fsolve(sistema1, [0.001, 2.5, 2.5, 3, 3], opts);

    if info_1 <= 0, continue; end

    Vds1_1 = x1(4);
    Vds1_2 = x1(5);
    Vov1   = x1(2) - Vt;

    % Saturación
    if (x1(2) <= Vt || x1(3) <= Vt || Vds1_1 < Vov1 || Vds1_2 < (x1(3)-Vt))
        continue;
    end

    % Restricción rds > 70 Ω en Etapa 1
    rds1_val  = Vds1_1 / x1(1);
    rdso1_val = Vds1_2 / x1(1);
    if (rds1_val < 70 || rdso1_val < 70), continue; end

    % Criterio de máxima excursión simétrica
    Vds_ideal = (Vdd - Vov1) / 2;
    error_excursion = abs(Vds1_1 - Vds_ideal);

    if error_excursion < min_error_excursion
        min_error_excursion = error_excursion;
        opt_Rs1 = Rs1;
        mejores_datos1 = x1;
    end
end

if isempty(mejores_datos1)
    error('No se encontró un Rs1 válido que cumpla rds > 70 Ω para la Etapa 1.');
end

% =========================================================================
% RESULTADOS FINALES Y CONDICIONES DE POLARIZACIÓN
% =========================================================================
fprintf('\n================ RESUMEN DE DISEÑO ÓPTIMO ================\n');
fprintf('Rref : %.2f Ω\n', opt_Rref);
fprintf('Rs2  : %.2f Ω\n', opt_Rs2);
fprintf('Rs1  : %.2f Ω\n', opt_Rs1);
fprintf('----------------------------------------------------------\n');

datos_ref = opt_xp;
datos_etapa2 = opt_x2;

fprintf('\n================ CONDICIONES DE POLARIZACIÓN ================\n');
fprintf('1. Rama de Referencia:\n');
fprintf('   Idsp = %.4f mA | Vgsp = %.4f V | rds_p = %.2f Ω\n', ...
        datos_ref(1)*1000, datos_ref(2), (datos_ref(2))/datos_ref(1));

fprintf('2. Etapa 1:\n');
fprintf('   Ids1  = %.4f mA\n', mejores_datos1(1)*1000);
fprintf('   M1_1: Vgs1_1 = %.4f V | Vds1_1 = %.4f V | rds1_1 = %.2f Ω\n', ...
        mejores_datos1(2), mejores_datos1(4), mejores_datos1(4)/mejores_datos1(1));
fprintf('   M1_2: Vgs1_2 = %.4f V | Vds1_2 = %.4f V | rds1_2 = %.2f Ω\n', ...
        mejores_datos1(3), mejores_datos1(5), mejores_datos1(5)/mejores_datos1(1));

fprintf('3. Etapa 2:\n');
fprintf('   Ids2  = %.4f mA\n', datos_etapa2(1)*1000);
fprintf('   M2_1: Vgs2_1 = %.4f V | Vds2_1 = %.4f V | rds2_1 = %.2f Ω\n', ...
        datos_etapa2(2), datos_etapa2(4), datos_etapa2(4)/datos_etapa2(1));
fprintf('   M2_2: Vgs2_2 = %.4f V | Vds2_2 = %.4f V | rds2_2 = %.2f Ω\n', ...
        datos_etapa2(3), datos_etapa2(5), datos_etapa2(5)/datos_etapa2(1));
fprintf('=============================================================\n');

% --- APARTE: CÁLCULO E IMPRESIÓN DE rds Y mu ---
gmp_val  = 2 * K * (datos_ref(2) - Vt);
rdsp_val = (Vdd - datos_ref(2)) / datos_ref(1);
mup_val  = gmp_val * rdsp_val;

gm1_1    = 2 * K * (mejores_datos1(2) - Vt);
gm1_2    = 2 * K * (mejores_datos1(3) - Vt);
rds1_1   = mejores_datos1(4) / mejores_datos1(1);
rds1_2   = mejores_datos1(5) / mejores_datos1(1);
mu1_1    = gm1_1 * rds1_1;
mu1_2    = gm1_2 * rds1_2;

gm2_1    = 2 * K * (datos_etapa2(2) - Vt);
gm2_2    = 2 * K * (datos_etapa2(3) - Vt);
rds2_1   = datos_etapa2(4) / datos_etapa2(1);
rds2_2   = datos_etapa2(5) / datos_etapa2(1);
mu2_1    = gm2_1 * rds2_1;
mu2_2    = gm2_2 * rds2_2;

fprintf('\n================ PARÁMETROS DE PEQUEÑA SEÑAL (rds Y mu) ================\n');
fprintf('1. Transistor de Referencia (p-MOSFET):\n');
fprintf('   rds = %.2f Ω | gm = %.4f mS | mu = %.2f\n', rdsp_val, gmp_val*1000, mup_val);

fprintf('2. Etapa 1:\n');
fprintf('   Transistor Superior (M1_2): rds = %.2f Ω | gm = %.4f mS | mu = %.2f\n', rds1_1, gm1_1*1000, mu1_1);
fprintf('   Transistor Inferior (M1_1): rds = %.2f Ω | gm = %.4f mS | mu = %.2f\n', rds1_2, gm1_2*1000, mu1_2);

fprintf('3. Etapa 2:\n');
fprintf('   Transistor Superior (M2_2): rds = %.2f Ω | gm = %.4f mS | mu = %.2f\n', rds2_1, gm2_1*1000, mu2_1);
fprintf('   Transistor Inferior (M2_1): rds = %.2f Ω | gm = %.4f mS | mu = %.2f\n', rds2_2, gm2_2*1000, mu2_2);
fprintf('=======================================================================\n');

% =========================================================================
% TABLA DE TENSIONES DE TERMINALES (Vg, Vd, Vs) E Ids
% =========================================================================

% --- 1. Transistor de Referencia (p-MOSFET: Mp) ---
Vs_p  = Vdd - datos_ref(1)*Rref;
Vg_p  = Vp;
Vd_p  = Vp;
Ids_p = datos_ref(1);

% --- 2. Etapa 1 ---
% Transistor Inferior (M1_1)
Vs_11  = 0;
Vg_11  = Vg1;
Vd_11  = mejores_datos1(4); % Vds1_1
Ids_11 = mejores_datos1(1);

% Transistor Superior / Carga Cascode (M1_2)
Vs_12  = Vd_11;
Vg_12  = Vp;
Vd_12  = Vdd - Ids_11 * Rd1;
Ids_12 = mejores_datos1(1);

% --- 3. Etapa 2 ---
% Transistor Inferior (M2_1)
Vs_21  = datos_etapa2(1) * opt_Rs2;
Vg_21  = Vg2;
Vd_21  = Vs_21 + datos_etapa2(4); % Vs + Vds2_1
Ids_21 = datos_etapa2(1);

% Transistor Superior / Carga Cascode (M2_2)
Vs_22  = Vd_21;
Vg_22  = Vp;
Vd_22  = Vdd - Ids_21 * Rd2;
Ids_22 = datos_etapa2(1);

% --- IMPRESIÓN DE LA TABLA ---
fprintf('\n================== TABLA DE TERMINALES DC (Vg, Vd, Vs, Ids) ==================\n');
fprintf('%-25s | %-10s | %-10s | %-10s | %-12s\n', 'Transistor', 'Vg [V]', 'Vd [V]', 'Vs [V]', 'Ids [mA]');
fprintf('------------------------------------------------------------------------------\n');
fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f | %-12.4f\n', 'Referencia (Mp)',           Vg_p,  Vd_p,  Vs_p,  Ids_p * 1000);
fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f | %-12.4f\n', 'Etapa 1 - Inferior (M1_1)', Vg_11, Vd_11, Vs_11, Ids_11 * 1000);
fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f | %-12.4f\n', 'Etapa 1 - Superior (M1_2)', Vg_12, Vd_12, Vs_12, Ids_12 * 1000);
fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f | %-12.4f\n', 'Etapa 2 - Inferior (M2_1)', Vg_21, Vd_21, Vs_21, Ids_21 * 1000);
fprintf('%-25s | %-10.4f | %-10.4f | %-10.4f | %-12.4f\n', 'Etapa 2 - Superior (M2_2)', Vg_22, Vd_22, Vs_22, Ids_22 * 1000);
fprintf('==============================================================================\n');
