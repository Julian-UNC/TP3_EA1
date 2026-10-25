clear all; close all; clc;

K   = 0.056;  % Constante del MOSFET [A/V^2]
Vt  = 2;         % Tensión de umbral [V]
Vdd = 12;        % Tensión de alimentación [V]
Rg1 = 2000000;   % 2 MΩ
Rg2 = 1000000;   % 1 MΩ
Rd  = 2200;      % 2.2 kΩ

Vg1 = (Vdd * Rg2) / (Rg1 + Rg2);  % Vg1 = 4 V


% Rangos de búsqueda para Vp y Rs
v_Vp = 0:0.25:12;       % Rango de Vp [V]
v_Rs = 100:200:10000;    % Rango de Rs [Ω]

mejor_error = inf;
puntos_optimos = [];

% Bucle de iteración
for Vp = v_Vp
    for Rs = v_Rs

        % Definición del sistema de ecuaciones: x = [Ids, Vgs1, Vgs2, Vds1, Vds2]
        sistema = @(x) [
            x(1) - K * (x(2) - Vt)^2;                                % Ecuación 1: Ids1
            x(1) - K * (x(3) - Vt)^2;                                % Ecuación 2: Ids2
            x(2) - (Vg1 - (x(5) - x(1)*Rs));                         % Ecuación 3: Vgs1 = Vg1 - Vd2
            x(3) - (Vp - x(1)*Rs);                                   % Ecuación 4: Vgs2 = Vp - Vs2
            Vdd - x(1)*(Rd + Rs) - x(4) - x(5)                      % Ecuación 5: Malla Vdd
        ];

        % Semilla / Valor inicial de prueba [Ids, Vgs1, Vgs2, Vds1, Vds2]
        x0 = [0.001, 3, 3, 4, 3];

        % Opciones del resolvedor para silenciar mensajes en consola
        options = optimset('Display', 'off');

        % Resolver el sistema de ecuaciones
        [x, fval, info] = fsolve(sistema, x0, options);

        % Desempaquetar resultados
        Ids  = x(1);
        Vgs1 = x(2);
        Vgs2 = x(3);
        Vds1 = x(4);
        Vds2 = x(5);

        % Validar que la solución sea físicamente válida y opere en saturación:
        % 1. Corriente positiva
        % 2. Ambos MOSFETs por encima del umbral (Vgs > Vt)
        % 3. Ambos MOSFETs en región de saturación (Vds >= Vgs - Vt)
        if (info == 1) && (Ids > 0) && (Vgs1 > Vt) && (Vgs2 > Vt) && ...
           (Vds1 >= (Vgs1 - Vt)) && (Vds2 >= (Vgs2 - Vt))

            % Tensión óptima ideal para Máxima Excursión Simétrica en Vds1
            %Vds1_opt = (Vdd - Ids*(Rd + Rs) - Vds2) / 2;
            %Vds1_opt = Vdd/2;
            Vds1_opt = (Vdd-Vds2)/2;

            % Error respecto a la máxima excursión
            error_excursion = abs(Vds1 - Vds1_opt);

            % Guardar la combinación con menor error
            if error_excursion < mejor_error
                mejor_error = error_excursion; %actualiza el mejor error
                puntos_optimos = [Vp, Rs, Ids, Vgs1, Vgs2, Vds1, Vds2, Vds1_opt];
            endif
        endif

    endfor
endfor

% =========================================================================
% MOSTRAR RESULTADOS EN PANTALLA
% =========================================================================
if !isempty(puntos_optimos)
    fprintf('\n=== RESULTADOS MÁXIMA EXCURSIÓN SIMÉTRICA ETAPA 1 ===\n');
    fprintf('Vp óptimo : %.3f V\n', puntos_optimos(1));
    fprintf('Rs óptimo : %.2f Ohm\n', puntos_optimos(2));
    fprintf('--------------------------------------------\n');
    fprintf('Ids  : %.4f mA \n', puntos_optimos(3)*1000);
    fprintf('Vgs1 : %.3f V\n', puntos_optimos(4));
    fprintf('Vgs2 : %.3f V\n', puntos_optimos(5));
    fprintf('Vds1 : %.3f V (Punto de operación)\n', puntos_optimos(6));
    fprintf('Vds2 : %.3f V\n', puntos_optimos(7));
    fprintf('Vds1_optimo (Mitad de rama): %.3f V\n', puntos_optimos(8));
    fprintf('Error en excursión: %.6f V\n', mejor_error);
else
    fprintf('No se encontró una solución válida dentro de los rangos especificados.\n');
endif

% Valores óptimos encontrados para Máxima Excursión Simétrica
Vp_opt = puntos_optimos(1);   % [V]
Rs_opt = puntos_optimos(2);  % [Ohm]

sistema = @(x) [
    x(1) - K * (x(2) - Vt)^2;         % Ids - K*(Vgs1-Vt)^2 = 0
    x(1) - K * (x(3) - Vt)^2;         % Ids - K*(Vgs2-Vt)^2 = 0
    x(2) - (Vg1 - (x(5) - x(1)*Rs_opt)); % Vgs1 = Vg1 - Vd2
    x(3) - (Vp_opt - x(1)*Rs_opt);    % Vgs2 = Vp - Vs2
    Vdd - x(1)*(Rd + Rs_opt) - x(4) - x(5) % Malla Vdd
];

x0 = [0.002, 2.5, 2.5, 3.0, 4.0];
[x, ~, info] = fsolve(sistema, x0);

if info == 1
    Ids_q  = x(1);  % Corriente en Amperios (si K está en mA/V^2, Ids estará en mA)
    Vgs1_q = x(2);
    Vgs2_q = x(3);
    Vds1_q = x(4);
    Vds2_q = x(5);
else
    error('El sistema no convergió con los valores óptimos provistos.');
endif

gm1 = 2*K*(Vgs1_q-Vt);
gm2 = 2*K*(Vgs2_q-Vt);

% =========================================================================
% 3. CÁLKULO E IMPRESIÓN DE LOS PARÁMETROS FINALES DE CADA TRANSISTOR
% =========================================================================
fprintf('\n======================================================\n');
fprintf('  PARÁMETROS OPERATIVOS FINALES (PUNTO Q) DE LOS MOSFETs \n');
fprintf('======================================================\n');
fprintf(' Parámetros de entrada seleccionados:\n');
fprintf('   Vp = %.2f V\n', Vp_opt);
fprintf('   Rs = %.2f Ohm\n', Rs_opt);
fprintf('------------------------------------------------------\n');
fprintf(' TRANSISTOR 1 (M1):\n');
fprintf('   Id1  = %.4f mA\n', Ids_q*1000);
fprintf('   Vgs1 = %.3f V\n', Vgs1_q);
fprintf('   Vds1 = %.3f V\n', Vds1_q);
fprintf('   gm1 = %.3f mS\n', gm1*1000);
fprintf('------------------------------------------------------\n');
fprintf(' TRANSISTOR 2 (M2):\n');
fprintf('   Id2  = %.4f mA\n', Ids_q*1000);
fprintf('   Vgs2 = %.3f V\n', Vgs2_q);
fprintf('   Vds2 = %.3f V\n', Vds2_q);
fprintf('   gm2 = %.3f mS\n', gm2*1000);
fprintf('======================================================\n\n');

% Vector para graficar la parábola característica
Vgs_vec = linspace(Vt, Vt + 2, 300);
Ids_vec = K * (Vgs_vec - Vt) .^ 2;

%%%   CURVAS CARACTERÍSTICAS DE LOS TRANSITORES

% --- Gráfica Transistor 1 ---
figure(1);
plot(Vgs_vec, Ids_vec, 'b-', 'LineWidth', 2);
hold on;
line([Vgs1_q Vgs1_q], [0 Ids_q], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line([Vt Vgs1_q], [Ids_q Ids_q], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
plot(Vgs1_q, Ids_q, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% --- LÍMITES DE LOS EJES (Rango Reducido) ---
xlim([Vgs1_q - 0.2, Vgs1_q + 0.2]);   % Rango X centrado en Vgs1_q (ej. ±0.2V)
ylim([0, Ids_q * 2.5]);                % Rango Y desde 0 hasta un poco más del Id de reposo

grid on;
xlabel('V_{gs1} [V]');
ylabel('I_{ds} [A]');
title(sprintf('Punto Q - Transistor 1 (V_{gs1}=%.3fV, Id=%.3fmA)', Vgs1_q, Ids_q));
legend('Curva I_{ds}', 'Intersección V_{gs1}', 'Intersección Id', 'Location', 'northwest');
hold off;

% --- Gráfica Transistor 2 ---
figure(2);
plot(Vgs_vec, Ids_vec, 'b-', 'LineWidth', 2);
hold on;
line([Vgs2_q Vgs2_q], [0 Ids_q], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line([Vt Vgs2_q], [Ids_q Ids_q], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
plot(Vgs2_q, Ids_q, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% --- LÍMITES DE LOS EJES (Rango Reducido) ---
xlim([Vgs2_q - 0.2, Vgs2_q + 0.2]);   % Rango X centrado en Vgs2_q
ylim([0, Ids_q * 2.5]);                % Rango Y adaptado a la corriente de reposo

grid on;
xlabel('V_{gs2} [V]');
ylabel('I_{ds} [A]');
title(sprintf('Punto Q - Transistor 2 (V_{gs2}=%.3fV, Id=%.3fmA)', Vgs2_q, Ids_q));
legend('Curva I_{ds}', 'Intersección V_{gs2}', 'Intersección Id', 'Location', 'northwest');
hold off;

%%%     RECTAS DE CARGA DE LOS TRANSISTORES

Zin2 = 556.79;

figure(3);

Vds_1 = linspace(0, Vdd-Vds2_q, 100);
Ids_1 = (Vdd-Vds2_q)/(Rd+Rs_opt) - Vds_1/(Rd+Rs_opt);
Rca1  = Rs_opt + 1/(1/Zin2+1/Rd);
ids_1 = Ids_q + (Vds1_q-Vds_1)/Rca1;

plot(Vds_1, Ids_1, 'm', 'LineWidth', 1.5);
hold on;
plot(Vds_1, ids_1, 'm', 'LineWidth', 1  );
hold on;
plot(Vds1_q, Ids_q, 'mo', 'markersize', 8, 'markerfacecolor', 'm');
title('Recta de Carga CC - CA - M1 (Primer Etapa)');
xlabel('Vds [V]'); ylabel('Id [A]');
grid on;
hold off;

figure(4);

Vds_2 = linspace(0, Vdd-Vds1_q, 100);
Ids_2 = (Vdd-Vds1_q)/(Rd+Rs_opt) - Vds_2/(Rd+Rs_opt);
Rca2  = Rs_opt + 1/(1/Zin2+1/Rd);
ids_2 = Ids_q + (Vds2_q-Vds_2)/Rca1;

plot(Vds_2, Ids_2, 'm', 'LineWidth', 1.5);
hold on;
plot(Vds_1, ids_1, 'm', 'LineWidth', 1.5);
hold on;
plot(Vds2_q, Ids_q, 'mo', 'markersize', 8, 'markerfacecolor', 'm');
title('Recta de Carga CC - CA - Mo1 (Primer Etapa)');
xlabel('Vds [V]'); ylabel('Id [A]');
grid on;
hold off;

disp('-- IMPEDANCIAS --');
rds = 7.5;
mu1 = gm1*rds
mu2 = gm2*rds

Raux1 = rds+(rds+Rs_opt*(mu2+1))*(mu1+1);
Zout = 1/(1/Rd+1/Raux1)

Zin = Rg1


