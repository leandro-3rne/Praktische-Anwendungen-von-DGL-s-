clear all;
close all;
clc;

% Parameter
m = 1.0;    % Masse in kg
L = 1.0;    % Länge in m
g = 9.81;   % Erdbeschleunigung
b = 0.5;    % Dämpfungskoeffizient [b_kritisch = 2 * m * sqrt(g/L)]

% Abgeleitete Grössen
omega0 = sqrt(g/L);                 % Eigenfrequenz (ungedämpft)
delta = b / (2*m);                  % Dämpfungskonstante
omegaD = sqrt(omega0^2 - delta^2);  % gedämpfte Eigenfrequenz

% Anfangsbedingungen
theta0 = pi/3;  % Startwinkel
dtheta0 = 0;    % Startgeschwindigkeit

% Zeitachse
t = linspace(0, 20, 2000);


% -- Analytische Lösung (homogen, Kleinwinkel) --

% Koeffizienten aus Anfangsbedingungen
A = theta0;
B = (dtheta0 + delta*theta0) / omegaD;

% Analytische homogene Lösung
theta_hom = exp(-delta*t) .* (A*cos(omegaD*t) + B*sin(omegaD*t));


% -- Numerische Lösung mit ode45 --

z0 = [theta0; dtheta0];
[t_num, z_num] = ode45(@(t,z) pendel_ode(t, z, b, m, L, g), t, z0);
theta_num = z_num(:,1);


% -- Plotten der DGL-Lösungen --

figure(1);
set(gcf, 'Position', [100, 100, 700, 500]);
plot(t, theta_hom, 'b', 'LineWidth', 2); hold on; % blau
plot(t_num, theta_num, 'r--', 'LineWidth', 1.5);    % rot
xlabel('Zeit t [s]');
ylabel('Winkel \theta [rad]');
title('Gedämpftes Pendel - Vergleich analytisch vs. numerisch');
legend('Analytisch (homogen, Kleinwinkel)', 'Numerisch (ode45, exakt)');
grid on;


% -- Animation --

figure(2);
set(gcf, 'Position', [820, 100, 700, 500]);
for i = 1:5:length(t_num)
    figure(2);
    clf;
    % Pendelposition berechnen
    x = L * sin(theta_num(i));
    y = -L * cos(theta_num(i));

    % Zeichnen
    plot([0, x], [0, y], 'k-', 'LineWidth', 3); hold on;
    plot(0, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(x, y, 'ro', 'MarkerSize', 20, 'MarkerFaceColor', 'r');

    axis equal;
    xlim([-L*1.2, L*1.2]);
    ylim([-L*1.2, 0.3]);
    title(sprintf('t = %.2f s,  theta = %.2f rad', t_num(i), theta_num(i)));
    grid on;
    drawnow;
end


% -- Hilfsfunktion ode45 --

% Umwandlugn DGL 2.Ordnung in zwei DGLs 1.Ordnung

function dz = pendel_ode(t, z, b, m, L, g)
    theta = z(1);
    dtheta = z(2);
    dz = zeros(2, 1);
    dz(1) = dtheta;
    dz(2) = -(b/m)*dtheta - (g/L)*sin(theta);
end

