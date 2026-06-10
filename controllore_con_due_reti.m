function Kfinal = controllore_con_due_reti(mu1,mu2,mu3)
close all;
clc;

s = tf('s');

% mu1=0.1981;
% mu2=1.5550;
% mu3=3.2470;


    mu1=0.2679
    mu2=1.0000
    mu3=3.7321

%% Plant dopo feedback linearization

P = 1/s^2;

%% Autovalori della Laplaciana
lambda1 = mu1;
lambda2 = mu2;
lambda3= mu3;
% Il caso critico è rappresentato dall'autovalore più alto quindi lambda2

L1 = lambda1 * P;
L2 = lambda2 * P;
L3 = lambda3 * P;

figure;
bode(L1, L2, L3);
grid on;
legend('\lambda=mu1','\lambda=mu2');
title('Open-loop non compensato');

%% Frequenza di crossover desiderata
wc = 1; % rad/s

[mag1, ph1] = bode(L1, wc);
[mag2, ph2] = bode(L2, wc);
[mag3, ph3] = bode(L3, wc);
ph1 = ph1(:);
ph2 = ph2(:);
ph3 = ph3(:);

%% Calcolo fase necessaria - DUE RETI ANTICIPATRICI
PM_des = 70;
PM_att = 180 + ph3; % caso peggiore lambda3=3, ph2=-180 -> PM_att=0
phi_m_tot = PM_des - PM_att; % fase totale da aggiungere

% Margine di sicurezza
phi_m_tot = phi_m_tot + 10;

% Con due reti anticipatrici, divido la fase a metà
phi_m_each = phi_m_tot / 2;

fprintf('Fase totale da aggiungere: %.1f gradi\n', phi_m_tot);
fprintf('Fase per ogni rete: %.1f gradi\n', phi_m_each);

%% Calcolo parametri di ogni rete
phi_m_rad = deg2rad(phi_m_each);
alpha = (1 - sin(phi_m_rad)) / (1 + sin(phi_m_rad));
tau = 1 / (wc * sqrt(alpha));

fprintf('alpha = %.4f\n', alpha);
fprintf('tau   = %.4f s\n', tau);

% Singola rete anticipatrice
Klead_single = (1 + s*tau) / (1 + s*alpha*tau)

% Due reti anticipatrici in cascata
Klead = Klead_single^2

%% Calcolo guadagno K per imporre wc come frequenza di crossover
K0 = 1;
L_test = lambda3 * K0 * Klead * P;
[mag_wc, ~] = bode(L_test, wc);
mag_wc = mag_wc(:);

K = 1 / mag_wc;
Kfin = K * Klead

fprintf('Guadagno K = %.4f\n', K);

%INSERIMENTO DELLA
%RITARDATRICE---------------------------------------------

%% Rete ritardatrice tarata sul loop rho1 con le lead

L1_preLag = lambda1 * Kfin * P;

% --- wc1_target scelto a 0.3 rad/s ---
% la fase di rho1 li' vale -130 deg => PM = 50 deg, accettabile
wc1_target = 0.3;%0.25

% --- beta: quanto guadagno serve per portare il crossover a wc1_target ---
% vuoi |N_lag(j*wc1_target) * L1_preLag(j*wc1_target)| = 1
% a bassa frequenza N_lag ~ beta, quindi:
% beta * |L1_preLag(j*wc1_target)| = 1
mag_L1 = abs(evalfr(L1_preLag, 1j * wc1_target));
beta = 1 / mag_L1;
fprintf('|L1_preLag| a wc1_target: %.4f\n', mag_L1);
fprintf('Beta necessario: %.4f\n', beta);

% --- tau1: zero della lag almeno una decade sotto wc1_target ---
% cosi' a wc1_target la lag si comporta gia' come guadagno puro beta
% e non sporca la fase
tau1 = 10 / wc1_target;
tau2 = tau1 / beta;
fprintf('tau1 = %.4f s\n', tau1);
fprintf('tau2 = %.4f s\n', tau2);

% --- rete lag ---
N_lag = (1 + tau1*s) / (1 + tau2*s);

% --- controllore finale ---
%Kg=1.5;
Kfinal = N_lag * Kfin;

rho=[lambda1,lambda2,lambda3];
w0=0.4;

%% Verifica finale con lag
fprintf('\n=== VERIFICA CON LAG ===\n\n');
for i = 1:length(rho)
    Li = rho(i) * Kfinal * P;
    [Gmi, Pmi, ~, Wcpi] = margin(Li);

    Si   = 1 / (1 + Li);
    Ti   = feedback(Li, 1);
    info = stepinfo(Ti);
    Si_at_w0 = 20*log10(abs(evalfr(Si, 1j*w0)));

    fprintf('--- Modo rho = %.4f ---\n', rho(i));
    fprintf('  PM       = %.2f deg   (spec: >= 45 deg)\n', Pmi);
    fprintf('  wc       = %.4f rad/s\n', Wcpi);
    fprintf('  |S(jw0)| = %.2f dB    (spec: < -20 dB)\n', Si_at_w0);
    fprintf('  Sovraelong = %.1f%%   (spec: <= 20%%)\n', info.Overshoot);
    fprintf('  T_assest   = %.2f s   (spec: <= 25 s)\n\n', info.SettlingTime);
end


%% Loop compensati
L1c = lambda1 * Kfinal * P;
L2c = lambda2 * Kfinal * P;
L3c = lambda3 * Kfinal * P;

%% Bode con margini
figure;
margin(L1c);
hold on;
margin(L2c);
hold on;
margin(L3c);
grid on;
legend('\lambda=1','\lambda=2');
title('Loop compensato - due reti anticipatrici');

%% Verifica numerica margini
[Gm1, Pm1, ~, Wcp1] = margin(L1c);
[Gm2, Pm2, ~, Wcp2] = margin(L2c);

fprintf('\n--- Margini di stabilità ---\n');
fprintf('lambda=1: GM=%.2f dB, PM=%.2f deg, wc=%.4f rad/s\n', ...
        20*log10(Gm1), Pm1, Wcp1);

fprintf('lambda=2: GM=%.2f dB, PM=%.2f deg, wc=%.4f rad/s\n', ...
        20*log10(Gm2), Pm2, Wcp2);

%% Risposta al gradino in anello chiuso
T1 = feedback(L1c, 1);
T2 = feedback(L2c, 1);
T3 = feedback(L3c, 1);

figure;
step(T1, T2, T3);
grid on;
legend('\lambda=mu1','\lambda=mu2','\lambda=mu3');
title('Closed-loop response - due reti anticipatrici');

%% Funzione di sensitività
S1 = 1 / (1 + L1c);
S2 = 1 / (1 + L2c);
S3 = 1 / (1 + L3c);

figure;
bodemag(S1, S2, S3);
grid on;
legend('\lambda=1','\lambda=2');
title('Sensitivity function');

%% Nyquist - zoom sul punto critico
figure;
nyquist(L1c, L2c);
axis([-2 2 -2 2]);
grid on;
legend('\lambda=1','\lambda=2');
title('Nyquist plot - due reti anticipatrici');

fprintf('PM lambda=2: %.2f deg\n', Pm2);

[~, ph_L2c_at_wc] = bode(L2c, wc);
[~, ph_L1c_at_wc] = bode(L1c, wc);

ph_L2c_at_wc = ph_L2c_at_wc(:);
ph_L1c_at_wc = ph_L1c_at_wc(:);

fprintf('\n--- Diagnostica fase a wc=%.2f ---\n', wc);
fprintf('Fase L2c a wc: %.2f deg  (PM contribuito: %.2f)\n', ...
        ph_L2c_at_wc, 180+ph_L2c_at_wc);

fprintf('Fase L1c a wc: %.2f deg  (PM contribuito: %.2f)\n', ...
        ph_L1c_at_wc, 180+ph_L1c_at_wc);

fprintf('Crossover effettivo L2c: %.4f rad/s\n', Wcp2);
fprintf('Crossover effettivo L1c: %.4f rad/s\n', Wcp1);

fprintf('Frequenza massima fase rete: %.4f rad/s\n', ...
        1/(tau*sqrt(alpha)));

fprintf('\n=== VERIFICA CONTROLLORE SU TUTTI I MODI ===\n\n');

rho=[lambda1,lambda2,lambda3];

for i = 1:length(rho)
    Li = rho(i) * Kfinal * P;
    [Gmi, Pmi, ~, Wcpi] = margin(Li);
    
    % Sensitività a omega0
    Si = 1/(1 + Li);
    w0 = 0.4;
    Si_at_w0 = abs(evalfr(Si, 1j*w0));
    
    % Risposta al gradino
    Ti = feedback(Li, 1);
    info = stepinfo(Ti);
    
    fprintf('--- Modo rho = %.4f ---\n', rho(i));
    fprintf('  GM     = %.2f dB      (spec: >= 6 dB)\n', 20*log10(Gmi));
    fprintf('  PM     = %.2f deg     (spec: >= 45 deg)\n', Pmi);
    fprintf('  wc     = %.4f rad/s\n', Wcpi);
    fprintf('  |S(jw0)| = %.2f dB   (spec: < -20 dB)\n', 20*log10(Si_at_w0));
    fprintf('  Sovraelong = %.1f%%   (spec: <= 20%%)\n', info.Overshoot);
    fprintf('  T_assest   = %.2f s  (spec: <= 25 s)\n\n', info.SettlingTime);

end