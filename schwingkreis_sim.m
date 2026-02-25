clear all;
close all;
clc;

% Parameter
R = 2.0;    % Widerstand in Ohm  [R_kritisch = 2*sqrt(L/C)]
L = 1.0;    % Induktivität in H
C = 0.1;    % Kapazität in F
U0 = 3.0;   % Amplitude der Antriebsspannung in V (0 = freie Schwingung)
Omega = 1.0;  % Antriebsfrequenz in rad/s

% Abgeleitete Grössen
omega0 = 1 / sqrt(L*C);         % Eigenfrequenz (ungedämpft)
delta  = R / (2*L);             % Dämpfungskonstante

% Anfangsbedingungen
Q0  = 0.01;   % Startladung am Kondensator [C]
dQ0 = 0;      % Startstromstärke I(0) = dQ/dt [A]

% Zeitachse
t = linspace(0, 20, 2000);

% -- Analytische homogene Lösung (drei Fälle) --

if delta < omega0                           % (a) Unterdämpfung
    omegaD = sqrt(omega0^2 - delta^2);
    A = Q0;
    B = (dQ0 + delta*Q0) / omegaD;
    Q_hom = exp(-delta*t) .* (A*cos(omegaD*t) + B*sin(omegaD*t));

elseif abs(delta - omega0) < 1e-6           % (b) Kritische Dämpfung
    c1 = Q0;
    c2 = dQ0 + delta*Q0;
    Q_hom = (c1 + c2*t) .* exp(-delta*t);

else                                        % (c) Überdämpfung
    mu = sqrt(delta^2 - omega0^2);
    c1 = (Q0*(delta + mu) + dQ0) / (2*mu);
    c2 = (Q0*(mu - delta) - dQ0) / (2*mu);
    Q_hom = c1*exp((-delta+mu)*t) + c2*exp((-delta-mu)*t);
end

% -- Analytische partikuläre Lösung (inhomogener Fall, U0 ~= 0) --

Q_anal = Q_hom;  % Standardfall: homogene Lösung (U0=0)

if U0 ~= 0
    a0   = 1 / (L*C);
    a1   = R / L;
    c_inh = U0 / L;

    % Prüfe Resonanzkatastrophe: a0 = Omega^2 und R = 0
    if abs(a0 - Omega^2) < 1e-9 && abs(R) < 1e-9
        % Resonanzkatastrophe: Qp(t) = U0/(2*L*omega0) * t * sin(omega0*t)
        Qp   = (U0 / (2*L*omega0)) * t .* sin(omega0*t);
        dQp0 = U0 / (2*L*omega0);   % Qp'(0) = U0/(2*L*omega0)
    else
        % Normalfall: Qp(t) = (c_inh/r) * cos(Omega*t - zeta)
        re_part = a0 - Omega^2;
        im_part = a1 * Omega;
        r    = sqrt(re_part^2 + im_part^2);
        zeta = atan2(im_part, re_part);   % angle(a0 - Omega^2 + i*a1*Omega)
        Qp   = (c_inh / r) * cos(Omega*t - zeta);
        dQp0 = (c_inh / r) * (-Omega) * sin(-zeta);   % Qp'(0)
    end

    % Anfangsbedingungen des homogenen Anteils anpassen:
    % Q_h(0)  = Q0 - Qp(0)
    % Q_h'(0) = dQ0 - Qp'(0)
    Q0_h  = Q0  - Qp(1);
    dQ0_h = dQ0 - dQp0;

    % Homogenen Anteil mit angepassten AB neu berechnen
    if delta < omega0
        omegaD = sqrt(omega0^2 - delta^2);
        Ah = Q0_h;
        Bh = (dQ0_h + delta*Q0_h) / omegaD;
        Q_hom_adj = exp(-delta*t) .* (Ah*cos(omegaD*t) + Bh*sin(omegaD*t));
    elseif abs(delta - omega0) < 1e-6
        c1h = Q0_h;
        c2h = dQ0_h + delta*Q0_h;
        Q_hom_adj = (c1h + c2h*t) .* exp(-delta*t);
    else
        mu = sqrt(delta^2 - omega0^2);
        c1h = (Q0_h*(delta + mu) + dQ0_h) / (2*mu);
        c2h = (Q0_h*(mu - delta) - dQ0_h) / (2*mu);
        Q_hom_adj = c1h*exp((-delta+mu)*t) + c2h*exp((-delta-mu)*t);
    end

    % Vollständige analytische Lösung: allgemeine = homogen + partikulär
    Q_anal = Q_hom_adj + Qp;
end

% -- Numerische Lösung mit ode45 --

z0 = [Q0; dQ0];
[t_num, z_num] = ode45(@(t,z) rlc_ode(t, z, R, L, C, U0, Omega), t, z0);
Q_num = z_num(:,1);
I_num = z_num(:,2);    % Stromstärke I = dQ/dt

% Plotten der DGL-Lösungen
figure(1);
set(gcf, 'Position', [300, 485, 500, 400]);
plot(t, Q_anal, 'b', 'LineWidth', 2); hold on;
plot(t_num, Q_num, 'r--', 'LineWidth', 1.5);
xlabel('Zeit t [s]');
ylabel('Ladung Q [C]');
if U0 ~= 0
    title('Diff-Gleichung analytisch (hom+part) vs. numerisch');
    legend('Analytisch (Q_h + Q_p)', 'Numerisch (ode45, exakt)');
else
    title('Diff-Gleichung analytisch vs. numerisch');
    legend('Analytisch (homogen)', 'Numerisch (ode45, exakt)');
end
grid on;

% ---- Animationen & Plots ----

% Energie vorberechnen
% E_C = Q^2 / (2C)  -->  im Kondensator gespeicherte elektrische Energie
% E_L = L*I^2 / 2   -->  in der Spule gespeicherte magnetische Energie
E_C = Q_num.^2 / (2*C);
E_L = 0.5 * L * I_num.^2;
E_tot = E_C + E_L;
E_tot_max = max(E_tot) + 1e-12;
I_max = max(abs(I_num)) + 1e-12;

% Partikel-Setup fuer Stromanimation
% Schaltkreis-Rechteck (Uhrzeigersinn = positiver Strom):
%   unten  (0  - 8 ): x: 1->9,  y=2
%   rechts (8  - 14): x=9,      y: 2->8
%   oben   (14 - 22): x: 9->1,  y=8   (durch R)
%   links  (22 - 28): x=1,      y: 8->2  (durch C, dann L)
L_pfad = 32;
n_partikel = 8;
partikel_offset = 0;
dt_anim = mean(diff(t_num));
geschwindigkeit_max = 12.0;   % Einheiten pro Sekunde bei I_max

% Figure-Fenster positionieren (1920x1200)
figure(2);
set(gcf, 'Position', [800, 485, 500, 400]);
figure(3);
set(gcf, 'Position', [300, 20, 1000, 400]);

% Gemeinsame Animationsschleife
for i = 1:5:length(t_num)

    I_now = I_num(i);
    Q_now = Q_num(i);

    % Stromfarbe: blau = positiv, rot = negativ
    if I_now >= 0
        stromfarbe   = [0.1 0.35 1.0];
    else
        stromfarbe   = [1.0 0.15 0.15];
    end

    % Partikel-Position aktualisieren
    geschwindigkeit = (I_now / I_max) * geschwindigkeit_max;
    partikel_offset = mod(partikel_offset + geschwindigkeit * dt_anim * 2, L_pfad);

    % Energieverlauf
    figure(2);
    clf;
    plot(t_num(1:i), E_C(1:i), 'b', 'LineWidth', 1.5); hold on;
    plot(t_num(1:i), E_L(1:i), 'r', 'LineWidth', 1.5);
    plot(t_num(1:i), E_tot(1:i), 'g', 'LineWidth', 2);
    xlim([0, t_num(end)]);
    ylim([0, E_tot_max * 1.15]);
    xlabel('Zeit t [s]');
    ylabel('Energie [J]');
    title(sprintf('Energieverlauf - t = %.2f s', t_num(i)));
    legend('E_C  (Kondensator)', 'E_L  (Spule)', 'E_{total}', ...
           'Location', 'northeast', 'FontSize', 9);
    grid on;

    % Schaltkreis + Stromverlauf
    figure(3);
    clf;

    % Linke Haelfte: Schaltkreis
    subplot(1,2,1);
    cla; hold on; axis off;
    xlim([-2 13]); ylim([-1.5 11.5]);
    title(sprintf('t = %.2f s', t_num(i)), 'FontSize', 10);

    lf = [0.5 0.5 0.5];   % Leitungsfarbe grau

    % Leitung oben links (zu R)
    plot([1 2.5], [9 9], 'Color', lf, 'LineWidth', 2.5);
    % Leitung oben rechts (von R)
    plot([7.5 9], [9 9], 'Color', lf, 'LineWidth', 2.5);
    % Leitung rechts (oben -> unten, keine Komponente)
    plot([9 9], [1 9], 'Color', lf, 'LineWidth', 2.5);
    % Leitung unten links (zu L)
    plot([1 2.5], [1 1], 'Color', lf, 'LineWidth', 2.5);
    % Leitung unten rechts (von L)
    plot([7.5 9], [1 1], 'Color', lf, 'LineWidth', 2.5);
    % Leitung links oben (von oben-links zu C)
    plot([1 1], [6.2 9], 'Color', lf, 'LineWidth', 2.5);
    % Leitung links unten (von C zu unten-links)
    plot([1 1], [1 3.8], 'Color', lf, 'LineWidth', 2.5);

    % Widerstand R (oben, Zickzack von x=2.5 bis x=7.5, y=9)
    n_zick = 9;
    rx = linspace(2.5, 7.5, n_zick*2+1);
    ry = 9.0 * ones(1, n_zick*2+1);
    for k = 2:2:n_zick*2
        ry(k) = 9.45;
    end
    for k = 3:2:n_zick*2-1
        ry(k) = 8.55;
    end
    plot(rx, ry, 'Color', [0.85 0.3 0.0], 'LineWidth', 2.5);
    text(5.0, 9.7, sprintf('R=%.1f\\Omega', R), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.85 0.3 0.0]);

    % Kondensator C (links, zwei horizontale Platten bei x=1, y=4 bis y=6)
    plot([1 1], [3.8 4.5], 'Color', [0.1 0.3 0.9], 'LineWidth', 2.5);
    plot([1 1], [5.5 6.2], 'Color', [0.1 0.3 0.9], 'LineWidth', 2.5);
    plot([0.0 2.0], [4.5 4.5], 'Color', [0.1 0.3 0.9], 'LineWidth', 5);
    plot([0.0 2.0], [5.5 5.5], 'Color', [0.1 0.3 0.9], 'LineWidth', 5);
    text(-1.5, 5.0, sprintf('C=%.2fF', C), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.1 0.3 0.9], 'Rotation', 90);
    text(2.5, 5.0, sprintf('Q=%.3fC', Q_now), ...
        'FontSize', 8, 'Color', [0.3 0.3 0.8], 'HorizontalAlignment', 'left');

    % Spule L (unten, Bögen von x=2.5 bis x=7.5, y=1)
    n_boegen = 5;
    bogen_breite = 1.0;
    bogen_start_x = 3.0;
    for k = 1:n_boegen
        ang = linspace(pi, 0, 30);   % Bögen nach unten
        bx  = (bogen_start_x + (k-1)*bogen_breite) + bogen_breite/2 * cos(ang);
        by  = 1.0 - bogen_breite/2 * sin(ang);
        plot(bx, by, 'Color', [0.0 0.55 0.15], 'LineWidth', 2.5);
    end
    text(5.0, -1.2, sprintf('L=%.1fH', L), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.0 0.55 0.15]);

    % Spannungsquelle (rechts, gegenüber vom Kondensator), nur wenn U0 ~= 0
    if U0 ~= 0
        U_now = U0 * cos(Omega * t_num(i));   % momentane Spannung

        % Leitung rechts wird durch Spannungsquelle unterbrochen
        % Spannungsquelle zwischen y=3.8 und y=6.2 bei x=9
        % Verbindungsleitungen oben/unten zur Quelle (bereits gezeichnet als Leitung rechts)
        % -> überschreibe den mittleren Teil der rechten Leitung mit Spannungsquellensymbol

        % Kreis der Spannungsquelle
        theta_kreis = linspace(0, 2*pi, 60);
        r_kreis = 1.2;
        xk = 9 + r_kreis * cos(theta_kreis);
        yk = 5.0 + r_kreis * sin(theta_kreis);
        % Hintergrund weiss übermalen (damit Leitungslinie nicht durchscheint)
        fill(xk, yk, 'w', 'EdgeColor', 'none');
        plot(xk, yk, 'Color', [0.7 0.0 0.7], 'LineWidth', 2.5);

        % + und - Zeichen (+ oben, - unten, Konvention: + an oberem Anschluss)
        text(9, 6.0, '+', 'HorizontalAlignment', 'center', ...
            'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.7 0.0 0.7]);
        text(9, 4.0, '–', 'HorizontalAlignment', 'center', ...
            'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.7 0.0 0.7]);

        % Beschriftung: kompakt rechts neben dem Kreis
        text(10.5, 5.8, sprintf('U_0=%.1fV', U0), ...
            'HorizontalAlignment', 'left', 'FontSize', 9, ...
            'FontWeight', 'bold', 'Color', [0.7 0.0 0.7]);
        text(10.5, 5.0, sprintf('\\Omega=%.2f', Omega), ...
            'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', [0.7 0.0 0.7]);
        text(10.5, 4.2, sprintf('U=%.2fV', U_now), ...
            'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', [0.7 0.0 0.7]);

        % Sinuskurve im Kreis als Symbol
        t_sym = linspace(0, 2*pi, 40);
        x_sym = 9 + 0.9 * (t_sym - pi) / pi;
        y_sym = 5.0 + 0.55 * sin(t_sym);
        plot(x_sym, y_sym, 'Color', [0.7 0.0 0.7], 'LineWidth', 1.5);
    end

    % Fliessende Strom-Partikel
    for p = 1:n_partikel
        s = mod(partikel_offset + (p-1) * L_pfad/n_partikel, L_pfad);
        [px, py] = pfad_zu_xy(s);
        plot(px, py, 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', stromfarbe, 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
    end

    % Stromrichtungsinfo oben – kompakt, zwei Zeilen
    if I_now >= 0
        richtung_kurz = 'I > 0 (Uhrzeigersinn)';
    else
        richtung_kurz = 'I < 0 (Gegenuhrzeigersinn)';
    end
    text(5.0, 11.0, sprintf('%s  |I|=%.3fA', richtung_kurz, abs(I_now)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, ...
        'FontWeight', 'bold', 'Color', stromfarbe);

    % Rechte Haelfte: Stromverlauf
    subplot(1,2,2);
    plot(t_num(1:i), I_num(1:i), 'Color', [0.1 0.35 1.0], 'LineWidth', 1.5); hold on;
    plot(t_num(i), I_num(i), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', stromfarbe, 'MarkerEdgeColor', stromfarbe);
    yline(0, 'k--', 'LineWidth', 1);
    xlim([0, t_num(end)]);
    ylim([-I_max*1.2, I_max*1.2]);
    xlabel('Zeit t [s]');
    ylabel('Stromstärke I [A]');
    title(sprintf('Stromstärke I(t) = dQ/dt - t = %.2f s', t_num(i)));
    grid on;

    drawnow;
end

% Hilfsfunktion: Bogenlänge s -> (x,y) auf Rechteck-Pfad (Uhrzeigersinn)
% Uhrzeigersinn: unten-links -> unten-rechts -> oben-rechts -> oben-links -> zurück
%   0- 8: unten  x: 1->9, y=1
%   8-16: rechts y: 1->9, x=9
%   16-24: oben  x: 9->1, y=9
%   24-32: links y: 9->1, x=1
% Gesamtumfang = 8+8+8+8 = 32
function [x, y] = pfad_zu_xy(s)
    s = mod(s, 32);
    if s < 8                        % unten: links nach rechts  (durch L)
        x = 1 + s;
        y = 1;
    elseif s < 16                   % rechts: unten nach oben
        x = 9;
        y = 1 + (s - 8);
    elseif s < 24                   % oben: rechts nach links   (durch R)
        x = 9 - (s - 16);
        y = 9;
    else                            % links: oben nach unten    (durch C)
        x = 1;
        y = 9 - (s - 24);
    end
end

% Hilfsfunktion ode45
% Umwandlung DGL 2.Ordnung in zwei DGLs 1.Ordnung
function dz = rlc_ode(t, z, R, L, C, U0, Omega)
    Q  = z(1);
    dQ = z(2);
    dz = zeros(2,1);
    dz(1) = dQ;
    dz(2) = -(R/L)*dQ - (1/(L*C))*Q + (U0/L)*cos(Omega*t);
end