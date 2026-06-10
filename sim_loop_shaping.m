%% Initialization

% 'on'  => u_tot = u_ff + u_fb, R(s) cancellato dall'errore
% 'off' => solo feedback, R(s) deve essere attenuato da S_i^pin
use_feedforward = 'off';

% Switch disturbo sull'uscita 
% 'none'         => nessun disturbo
% 'common'       => d1=d2=d3, non deforma la formazione
% 'differential' => d_i diversi, eccita i modi di formazione
dist_type = 'none';

% Switch rumore sensore
use_noise = 'off';             % 'on', 'off'
sigma     = 0.01;             % deviazione standard [m] — GPS ~1cm

%Parametri di simulazione per il loop
T_sim = 300; 
dt = 0.005; 
t = 0:dt:T_sim; 
Nt = length(t);

% triangular formation
h1 = [0;   0  ]; 
h2 = [1;   0  ]; 
h3 = [0.5; sqrt(3)/2];
h  = {h1, h2, h3};

% offset wtr to node 1
deltah = [ h2 - h1 ; h3 - h1 ];
H0 = [ 0 ; deltah ];
% D1 matrix and M matrix
N = 3; D1 = [ - ones(N-1,1) , eye(N-1) ];
M = [ zeros(1,N-1) ; eye(N-1) ];

% type of reference trajectory
ref_type = 'circular';
r = zeros(2,Nt); rdot = zeros(2,Nt); rddot = zeros(2,Nt);

switch ref_type
    case 'constant'
        v_ref = [0.3; 0.2];
        for k = 1:Nt
            r(:,k)    = v_ref * t(k);
            rdot(:,k) = v_ref;
        end

    case 'circular'
        Rc = 2;  om = 0.15;
        for k = 1:Nt
            r(:,k)     = Rc * [ cos(om*t(k));  sin(om*t(k))];
            rdot(:,k)  = Rc*om * [-sin(om*t(k));  cos(om*t(k))];
            rddot(:,k) = Rc*om^2 * [-cos(om*t(k)); -sin(om*t(k))];
        end

    case 'sinusoidal'
        As = 1.5;  oms = 0.15;...
        for k = 1:Nt
            r(:,k)     = [As*sin(oms*t(k));  0.3*t(k)];
            rdot(:,k)  = [As*oms*cos(oms*t(k));  0.3];
            rddot(:,k) = [-As*oms^2*sin(oms*t(k));  0];
        end
end

% desired trajectory 

y_star = zeros(2,3,Nt);

% initial conditions
% y*_i(t) = r(t) + h_i - h_1
offset = [0.5; 0.4];
X = zeros(4,3,Nt);
for i = 1:3
    for k=1:Nt
        y_star(:,i,k) = r(:,k) + h{i} - h1;
    end
    % p0 = y_star(:,i,1) + offset;
    % X(:,i,1) = [p0(1); 0; p0(2); 0];
end

% % GRAFO 
% source = [1,2];
% dest = [2,3]; 
% G = graph([1,2],[2,3]); 
%plot(G)
% A = adjacency(G)'; 
% D = diag(sum(A,1)); 
% L = D - A; 

% source = [1,2];
% dest = [2,3]; 
% G = graph(source,dest); 
% A = full(adjacency(G))'; 
A =   [0 1 1;
       1 0 0;
       1 0 0];
D = diag(sum(A,2));
L = D - A;

%PINNING DI L
Pi = diag([1 0 0]);
gamma = 1; 
Lp = L + gamma*Pi;

%% Loop Shaping Design

% Lp diagonalization
% Vp eigenvectors matrix , mu eigenvalues diagonal matrix
[Vp, mu] = eig(Lp)
% extract the eigenvalues in a vector
rho = diag(mu);
% sort the eigenvalues, rho(1) = 0 etc... , 
% idxp contains the indexes of the sorted eigenvalues
[rho, idxp] = sort(rho);
% sort Vp columns in the same order of rho, otherwise the
% eigenvalues and the eigenvectors will be disaligned
Vp = Vp(:,idxp);

%% Controller design

% Specifiche calcolate da Fra :
% Margine di fase tra 45-60 
% Margine di ampiezza >= 6 db
% Larghezza di banda 0.5 - 3 rad/s
% tempo di assestamento ~25 s
% sovraelongazione <= 20%

s = tf('s');

% Plant dopo feedback linearization
P = 1/s^2;

lambda1 = mu(1,1)
lambda2 = mu(2,2)
lambda3 = mu(3,3)

% Il caso critico è rappresentato dall'autovalore più alto quindi lambda2

% L1 = lambda1 * P;
% L2 = lambda2 * P;
% L3 = lambda3 * P;

K=controllore_con_due_reti(lambda1,lambda2,lambda3);
close all;

S_pin = cell(3,1);
for i = 1:3
    S_pin{i} = 1 / (1 + mu(i,i) * K * P);
end

% L1c = lambda1 * K * P;
% L2c = lambda2 * K * P;
% L3c = lambda3 * K * P;
% 
% 
% T1 = feedback(L1c, 1);
% T2 = feedback(L2c, 1);
% T3 = feedback(L3c, 1);
% 
% 
% S1 = 1 / (1 + L1c);
% S2 = 1 / (1 + L2c);
% S3 = 1 / (1 + L3c);

% Verifica |S_i^pin(j*w0)| alla frequenza del riferimento
w0 = om;  % 0.4 rad/s
fprintf('\n|S_pin_i(jw0)| alla frequenza del riferimento (w0=%.2f rad/s):\n', w0);
for i = 1:3
    val = abs(evalfr(S_pin{i}, 1j*w0));
    fprintf('  rho_%d = %.4f:  |S| = %.4f (%.2f dB)\n', ...
        i, rho(i), val, 20*log10(val));
end


% --- Parametri disturbo sull'uscita ---
A_dist   = 0.1;                    % ampiezza [m]
om_dist  = 5;                     % stessa frequenza del riferimento
phi_dist = [0; pi/3; 2*pi/3];     % sfasamenti per caso differenziale

%% Guadagno dei loop modali a omega0
omega0 = om_dist;
fprintf('\nGuadagno loop modali a omega0=%.2f rad/s:\n', omega0);
for i = 1:3
    Li = rho(i) * K * P;
    L_at_w0 = abs(evalfr(Li, 1j*omega0));
    fprintf('  rho_%d: |L_i(jw0)| = %.4f (%.2f dB)\n', ...
        i, L_at_w0, 20*log10(L_at_w0));
end

% Verifica attenuazione riferimento circolare senza feedforward
% w0 = 0.4;  % rad/s moto circolare
% S1_at_w0 = abs(evalfr(S1, 1j*w0));
% S2_at_w0 = abs(evalfr(S2, 1j*w0));
% fprintf('|S1(jw0)| = %.4f (%.2f dB)\n', S1_at_w0, 20*log10(S1_at_w0));
% fprintf('|S2(jw0)| = %.4f (%.2f dB)\n', S2_at_w0, 20*log10(S2_at_w0));



%% --- Condizioni iniziali ---
% Agenti partono con offset rispetto alla posizione desiderata
offset = [0.5; 0.4];
X = zeros(4, 3, Nt);   % stato [px; vx; py; vy] per ogni agente
for i = 1:3
    p0          = y_star(:,i,1) + offset;
    X(:,i,1)    = [p0(1); 0; p0(2); 0];
end
%% Simulation Loop
% ============================================================
%  SIMULATION LOOP
%  Stato: X(:,i,k) = [px; vx; py; vy] per agente i al passo k
%  Post feedback-linearization => ogni agente e' un doppio integratore
%  nu_i = virtual input (asse x e y separati)
% ============================================================

% --- Stato interno del controllore K(s) per ogni agente e asse ---
% K e' dinamico => serve memoria tra i passi
Kd  = c2d(ss(K), dt, 'tustin');
Kss = Kd;
nk  = size(Kss.A, 1);          % ordine del controllore
xk  = zeros(nk, 3, 2);         % (ordine_K, n_agenti=3, n_assi=2)

% --- Preallocazione ---
Y  = zeros(2, 3, Nt);          % posizioni [px; py] per ogni agente e passo
FE = zeros(2, 3, Nt);          % formation error per ogni agente e passo
NU = zeros(2, 3, Nt);          % virtual input nu per ogni agente e passo

% --- Condizioni iniziali ---
for i = 1:3
    Y(:, i, 1) = X([1,3], i, 1);   % estrai posizione iniziale dallo stato
end

% ============================================================
%  MAIN LOOP
% ============================================================
for k = 1:Nt-1

    %  DISTURBO SULL'USCITA 
    %  y_i misurata = y_i + d_i
    % ----------------------------------------------------------
    switch dist_type
        case 'none'
            y_meas = Y(:,:,k);

        case 'common'
            % d_1 = d_2 = d_3 => eccita solo modo consenso
            % => non deforma la formazione (D1*1_N = 0, slide 28)
            d_c = A_dist * sin(om_dist * t(k));
            y_meas = Y(:,:,k) + d_c * ones(2,3);

        % case 'differential'
        %     % d_i diversi => eccita i modi di formazione
        %     % => deve essere attenuato da |S_i^pin(jw)| (slide 30)
        %     y_meas = Y(:,:,k);
        %     for i = 1:3
        %         d_i = A_dist * sin(om_dist * t(k) + phi_dist(i));
        %         y_meas(:,i) = y_meas(:,i) + [d_i; d_i];
        %     end

            case 'differential'
                y_meas = Y(:,:,k);
                for i = 1:3
                % disturbo diverso su x e y per ogni agente
                d_x = A_dist * sin(om_dist * t(k) + phi_dist(i));
                d_y = A_dist * sin(om_dist * t(k) + phi_dist(i) + pi/2);
                y_meas(:,i) = y_meas(:,i) + [d_x; d_y];
                end
    end

    % 2. Rumore sensore gaussiano (indipendente dal disturbo)
    switch use_noise
        case 'on'
            y_meas = y_meas + sigma * randn(2,3);
        % case 'off': nessuna aggiunta
    end

    % ----------------------------------------------------------
    %  1. ERRORE DI TRACKING LOCALE
    %     f_i(t) = y_i(t) - r(t) - (h_i - h_1)
    % ----------------------------------------------------------
    eta=zeros(2,3);
    for i = 1:3
        eta(:,i) = y_meas(:,i) - r(:, k) - (h{i} - h1);
    end

    % Formation error (senza disturbo, per monitoraggio)
    for i = 1:3
        FE(:,i,k) = Y(:,i,k) - r(:,k) - (h{i} - h1);
    end

    % ----------------------------------------------------------
    %  2. ERRORE PINNATO
    %     eps_i = sum_j a_ij [(y_j - y_i) - (h_j - h_i)]
    %             - gamma * pi_i * f_i(t)
    %     pi_i = 1 solo per nodo 1 (pinning), 0 altrimenti
    % ----------------------------------------------------------
    eps = zeros(2, 3);
    for i = 1:3
        % contributo dai vicini
        for j = 1:3
                eps(:,i) = eps(:,i) -Lp(i,j) * eta(:,j);
                    % ( (Y(:,j,k) - Y(:,i,k)) - (h{j} - h{i}) );
        end
    end
        % termine di pinning (solo agente 1)
    %     if i == 1
    %         eps(:,1) = eps(:,1) - gamma * FE(:,1,k);
    %     end
    % end

    % ----------------------------------------------------------
    %  3. CONTROLLORE K(s) 
    %     Per ogni agente i e ogni asse (x=1, y=2):
    %       nu  = C*xk + D*eps
    %       xk' = A*xk + B*eps
    % ----------------------------------------------------------
for i = 1:3
    for ax = 1:2
        e_in         = eps(ax, i);
        % uscita: y = C*x + D*u
        NU(ax, i, k) = Kss.C * xk(:,i,ax) + Kss.D * e_in;
        % aggiornamento stato discreto: x[k+1] = A*x[k] + B*u[k]
        xk(:,i,ax)   = Kss.A * xk(:,i,ax) + Kss.B * e_in;
    end
end
    % ----------------------------------------------------------
    %  4. FEEDFORWARD
    %  Con FF: u_tot = u_ff + u_fb
    %          u_ff = 1_N * U_ff,r  tale che Gp*1_N*U_ff,r = 1_N*R(s)
    %          => per Gp = 1/s^2  =>  u_ff,r = r_ddot(t)
    %  Senza FF: solo u_fb, R(s) entra nella dinamica dell'errore
    % ----------------------------------------------------------
    switch use_feedforward
        case 'on'
            for i = 1:3
                NU(:,i,k) = NU(:,i,k) + rddot(:,k);
            end
        case 'off'
            % nessuna aggiunta — R(s) deve essere attenuato da S_i^pin
    end
    % ----------------------------------------------------------
    %  5. INTEGRAZIONE — doppio integratore (post FL)
    %     px' = vx       vx' = nu_x
    %     py' = vy       vy' = nu_y
    %  (sostituire con dinamica nonlineare + legge FL per il caso reale)
    % ----------------------------------------------------------
    for i = 1:3
        nu_i = NU(:, i, k);

        X(1,i,k+1) = X(1,i,k) + dt * X(2,i,k);          % px
        X(2,i,k+1) = X(2,i,k) + dt * nu_i(1);            % vx
        X(3,i,k+1) = X(3,i,k) + dt * X(4,i,k);          % py
        X(4,i,k+1) = X(4,i,k) + dt * nu_i(2);            % vy

        Y(:, i, k+1) = X([1,3], i, k+1);
    end

end

% formation error all'ultimo passo
for i = 1:3
    FE(:, i, Nt) = Y(:, i, Nt) - r(:, Nt) - (h{i} - h1);
end




%% ============================================================
%  ANALISI FORMATION ERROR
% ============================================================

k_reg = round(0.8 * Nt);    % inizio regime: ultimi 20% della simulazione

fprintf('\n===== Formation Error a regime =====\n');
for i = 1:3
    fe_norm = squeeze(sqrt(FE(1,i,:).^2 + FE(2,i,:).^2));
    fprintf('Agente %d: max=%.4f m | rms=%.4f m\n', ...
        i, max(fe_norm(k_reg:end)), rms(fe_norm(k_reg:end)));
end

% --- Plot formation error nel tempo ---
figure('Name','Formation Error');
tiledlayout(2,1);

nexttile;
hold on;
colors_fe = {'b','r','g'};
for i = 1:3
    fe_norm = squeeze(sqrt(FE(1,i,:).^2 + FE(2,i,:).^2));
    plot(t, fe_norm, 'Color', colors_fe{i}, 'LineWidth', 1.2);
end
legend('Agente 1','Agente 2','Agente 3');
ylabel('||f_i(t)|| [m]'); xlabel('t [s]');
title('Formation error — intero transitorio');
grid on;

nexttile;
hold on;
for i = 1:3
    fe_norm = squeeze(sqrt(FE(1,i,:).^2 + FE(2,i,:).^2));
    plot(t(k_reg:end), fe_norm(k_reg:end), 'Color', colors_fe{i}, 'LineWidth', 1.2);
end
legend('Agente 1','Agente 2','Agente 3');
ylabel('||f_i(t)|| [m]'); xlabel('t [s]');
title('Formation error — regime (ultimi 20%)');
grid on;


%% ============================================================
%  ANALISI TRAIETTORIE
% ============================================================

figure('Name','Traiettorie');
hold on;
colors_fe = {'b','r','g'};

% --- traiettorie reali ---
for i = 1:3
    plot(squeeze(Y(1,i,:)), squeeze(Y(2,i,:)), ...
         'Color', colors_fe{i}, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('A%d reale',i));
end

% --- traiettorie desiderate: ricalcolate esplicitamente ---
for i = 1:3
    p_des_x = r(1,:) + h{i}(1) - h1(1);
    p_des_y = r(2,:) + h{i}(2) - h1(2);
    plot(p_des_x, p_des_y, '--', ...
         'Color', colors_fe{i}, 'LineWidth', 1.0, ...
         'DisplayName', sprintf('A%d desiderato',i));
end

% --- traiettoria del leader r(t) ---
plot(r(1,:), r(2,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'r(t)');

legend('Location','best');
grid on; axis equal;
title('Traiettorie reali vs desiderate');

%% ============================================================
%  SNAPSHOT FORMAZIONE
%  Verifica visiva che il triangolo sia mantenuto
% ============================================================

figure('Name','Snapshot formazione');
k_snaps = [round(Nt*0.25), round(Nt*0.5), round(Nt*0.75), Nt];
labels_snap = {'t=T/4','t=T/2','t=3T/4','t=T'};
colors_fe = {'b','r','g'};

for s = 1:4
    subplot(2,2,s); hold on;
    k_s = k_snaps(s);

    % --- posizioni reali ---
    for i = 1:3
        plot(Y(1,i,k_s), Y(2,i,k_s), 'o', ...
             'MarkerSize', 10, ...
             'MarkerFaceColor', colors_fe{i}, ...
             'MarkerEdgeColor', 'k', ...
             'DisplayName', sprintf('A%d reale',i));
    end

    % --- posizioni desiderate: r(t) + h_i - h_1 ---
    for i = 1:3
        p_des = r(:,k_s) + h{i} - h1;   % calcolo esplicito, non da y_star
        plot(p_des(1), p_des(2), 'x', ...
             'MarkerSize', 12, ...
             'LineWidth', 2, ...
             'Color', colors_fe{i}, ...
             'DisplayName', sprintf('A%d desiderato',i));
    end

    % --- triangolo reale ---
    idx = [1 2 3 1];
    xr = squeeze(Y(1,idx,k_s));
    yr = squeeze(Y(2,idx,k_s));
    plot(xr, yr, 'k-', 'DisplayName', 'triangolo reale');

    % --- triangolo desiderato ---
    xd = arrayfun(@(i) r(1,k_s) + h{i}(1) - h1(1), 1:3);
    yd_arr = arrayfun(@(i) r(2,k_s) + h{i}(2) - h1(2), 1:3);
    plot([xd xd(1)], [yd_arr yd_arr(1)], 'k--', 'DisplayName', 'triangolo desiderato');

    % --- centro di riferimento r(t) ---
    plot(r(1,k_s), r(2,k_s), 'k+', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', 'r(t)');

    title(sprintf('%s  (t=%.1fs)', labels_snap{s}, (k_s-1)*dt));
    grid on; axis equal;
    legend('Location','best','FontSize',7)
end





%% Formation errors


%% Bode plot

% try
%     figure('Name','[LS] Open-Loop Bode','Position',[50 50 700 500]);
%     margin(Ls);
%     title('[Loop Shaping]  Open-loop  L(s) = K(s) \cdot G(s)','FontSize',12);
%     grid on;
% catch
%     fprintf('(Bode plot skipped — Control System Toolbox not available)\n');
% end
% 
% %% Plots
% colors = {'#1f77b4','#ff7f0e','#2ca02c'};
% lw = 1.6;
% 
% % --- Figure: Trajectories ---
% figure('Name','[LS] Trajectories','Position',[50 50 700 600]);
% hold on; grid on; axis equal;
% title(sprintf('[Loop Shaping]  Trajectories — %s reference', ref_type),'FontSize',13);
% xlabel('p_x [m]');  ylabel('p_y [m]');
% for i = 1:3
%     plot(squeeze(Y(1,i,:)), squeeze(Y(2,i,:)), ...
%         'Color',colors{i},'LineWidth',lw,'DisplayName',sprintf('Agent %d',i));
%     plot(Y(1,i,1),  Y(2,i,1),  'o','Color',colors{i},'MarkerSize',8,'HandleVisibility','off');
%     plot(Y(1,i,end),Y(2,i,end),'s','Color',colors{i},'MarkerSize',8,'HandleVisibility','off');
% end
% plot(r(1,:),r(2,:),'k--','LineWidth',1,'DisplayName','Reference r(t)');
% tx = [Y(1,1,end) Y(1,2,end) Y(1,3,end) Y(1,1,end)];
% ty = [Y(2,1,end) Y(2,2,end) Y(2,3,end) Y(2,1,end)];
% plot(tx,ty,'k-','LineWidth',2,'HandleVisibility','off');
% legend('Location','best');
% 
% % --- Figure: Position tracking ---
% figure('Name','[LS] Position Tracking','Position',[50 50 900 500]);
% subplot(2,1,1); hold on; grid on;
% title('x-position tracking');  xlabel('t [s]');  ylabel('p_x [m]');
% for i = 1:3
%     plot(t,squeeze(Y(1,i,:)),'Color',colors{i},'LineWidth',lw,'DisplayName',sprintf('Agent %d',i));
%     plot(t,squeeze(y_star(1,i,:)),'--','Color',colors{i},'LineWidth',1,'HandleVisibility','off');
% end
% legend('Location','best');
% subplot(2,1,2); hold on; grid on;
% title('y-position tracking');  xlabel('t [s]');  ylabel('p_y [m]');
% for i = 1:3
%     plot(t,squeeze(Y(2,i,:)),'Color',colors{i},'LineWidth',lw,'DisplayName',sprintf('Agent %d',i));
%     plot(t,squeeze(y_star(2,i,:)),'--','Color',colors{i},'LineWidth',1,'HandleVisibility','off');
% end
% legend('Location','best');
% 
% % --- Figure: Formation errors ---
% figure('Name','[LS] Formation Errors','Position',[50 50 900 500]);
% subplot(2,2,1); plot(t,delta2(1,:),'Color',colors{2},'LineWidth',lw); grid on;
% title('\delta_{2,x}');  xlabel('t [s]');  ylabel('[m]');
% subplot(2,2,2); plot(t,delta2(2,:),'Color',colors{2},'LineWidth',lw); grid on;
% title('\delta_{2,y}');  xlabel('t [s]');  ylabel('[m]');
% subplot(2,2,3); plot(t,delta3(1,:),'Color',colors{3},'LineWidth',lw); grid on;
% title('\delta_{3,x}');  xlabel('t [s]');  ylabel('[m]');
% subplot(2,2,4); plot(t,delta3(2,:),'Color',colors{3},'LineWidth',lw); grid on;
% title('\delta_{3,y}');  xlabel('t [s]');  ylabel('[m]');
% sgtitle('[Loop Shaping]  Formation errors  \delta_i = (y_i - y_1) - (h_i - h_1)','FontSize',12);
% 
% % --- Figure: Physical inputs ---
% figure('Name','[LS] Physical Inputs','Position',[50 50 900 600]);
% for i = 1:3
%     subplot(3,2,2*i-1);
%     plot(t,squeeze(U(1,i,:)),'Color',colors{i},'LineWidth',lw);
%     grid on;  title(sprintf('Agent %d  —  u_x',i));  xlabel('t [s]');  ylabel('u_x');
%     subplot(3,2,2*i);
%     plot(t,squeeze(U(2,i,:)),'Color',colors{i},'LineWidth',lw);
%     grid on;  title(sprintf('Agent %d  —  u_y',i));  xlabel('t [s]');  ylabel('u_y');
% end
% sgtitle('[Loop Shaping]  Physical inputs u_i','FontSize',12);
% 
% % --- Figure: Formation error norms ---
% figure('Name','[LS] Error Norms','Position',[50 50 700 350]);
% hold on; grid on;
% plot(t,vecnorm(delta2),'Color',colors{2},'LineWidth',lw,'DisplayName','||\delta_2||');
% plot(t,vecnorm(delta3),'Color',colors{3},'LineWidth',lw,'DisplayName','||\delta_3||');
% xlabel('t [s]');  ylabel('[m]');
% title('[Loop Shaping]  Formation error norms');
% legend('Location','best');

% Errore relativo agente 2 vs agente 1
fe_rel_21 = squeeze(sqrt( ...
    (Y(1,2,:) - Y(1,1,:) - (h{2}(1)-h1(1))).^2 + ...
    (Y(2,2,:) - Y(2,1,:) - (h{2}(2)-h1(2))).^2 ));

% Errore relativo agente 3 vs agente 1  
fe_rel_31 = squeeze(sqrt( ...
    (Y(1,3,:) - Y(1,1,:) - (h{3}(1)-h1(1))).^2 + ...
    (Y(2,3,:) - Y(2,1,:) - (h{3}(2)-h1(2))).^2 ));

figure('Name','Errore relativo');
subplot(2,1,1); plot(t, fe_rel_21, 'b', 'LineWidth', 1.2);
ylabel('||y_2 - y_1 - (h_2-h_1)||'); xlabel('t [s]'); grid on;
title('Errore relativo agente 2 vs 1');

subplot(2,1,2); plot(t, fe_rel_31, 'r', 'LineWidth', 1.2);
ylabel('||y_3 - y_1 - (h_3-h_1)||'); xlabel('t [s]'); grid on;
title('Errore relativo agente 3 vs 1');


%% DIAGNOSI CONVERGENZA
fprintf('\n=== DIAGNOSI ===\n');
figure;
% 1. Formation error a regime
k_reg = round(0.8*Nt);
fprintf('\n1. Formation error a regime:\n');
for i = 1:3
    fe = squeeze(sqrt(FE(1,i,:).^2 + FE(2,i,:).^2));
    fprintf('   Agente %d: max=%.4f | rms=%.4f | finale=%.4f\n', ...
        i, max(fe(k_reg:end)), rms(fe(k_reg:end)), fe(end));
end

% 2. Errore relativo (formazione)
fprintf('\n2. Errore relativo a regime:\n');
fe_rel_21 = squeeze(sqrt( ...
    (Y(1,2,:)-Y(1,1,:)-(h{2}(1)-h1(1))).^2 + ...
    (Y(2,2,:)-Y(2,1,:)-(h{2}(2)-h1(2))).^2 ));
fe_rel_31 = squeeze(sqrt( ...
    (Y(1,3,:)-Y(1,1,:)-(h{3}(1)-h1(1))).^2 + ...
    (Y(2,3,:)-Y(2,1,:)-(h{3}(2)-h1(2))).^2 ));
fprintf('   A2-A1: max=%.4f | finale=%.4f\n', ...
    max(fe_rel_21(k_reg:end)), fe_rel_21(end));
fprintf('   A3-A1: max=%.4f | finale=%.4f\n', ...
    max(fe_rel_31(k_reg:end)), fe_rel_31(end));

% 3. Tracking assoluto nodo 1
fprintf('\n3. Tracking assoluto nodo 1:\n');
track_1 = zeros(Nt,1);
for k = 1:Nt
    track_1(k) = sqrt((Y(1,1,k)-r(1,k))^2 + (Y(2,1,k)-r(2,k))^2);
end
fprintf('\n3. Tracking assoluto nodo 1:\n');
fprintf('   max=%.4f | rms=%.4f | finale=%.4f\n', ...
    max(track_1(k_reg:end)), rms(track_1(k_reg:end)), track_1(end));

% 4. Velocita' degli agenti a regime (devono seguire rdot)
fprintf('\n4. Velocita'' agente 1 a fine sim vs rdot atteso:\n');
fprintf('   vx reale=%.4f | vx atteso=%.4f\n', X(2,1,end), rdot(1,end));
fprintf('   vy reale=%.4f | vy atteso=%.4f\n', X(4,1,end), rdot(2,end));

% 5. Virtual input a regime (non deve esplodere)
fprintf('\n5. Virtual input a regime:\n');
for i = 1:3
    nu_norm = squeeze(sqrt(NU(1,i,:).^2 + NU(2,i,:).^2));
    fprintf('   Agente %d: max=%.4f | finale=%.4f\n', ...
        i, max(nu_norm(k_reg:end)), nu_norm(end));
end


%% Animazione formazione
% Mostra i 3 agenti, il triangolo, e il reference r(t)

figure('Name', 'Animazione Formazione', 'Position', [100 100 700 700]);

% Calcola i limiti degli assi una volta sola
all_x = [squeeze(Y(1,1,:)); squeeze(Y(1,2,:)); squeeze(Y(1,3,:)); r(1,:)'];
all_y = [squeeze(Y(2,1,:)); squeeze(Y(2,2,:)); squeeze(Y(2,3,:)); r(2,:)'];
ax_lim = [min(all_x)-0.5, max(all_x)+0.5, min(all_y)-0.5, max(all_y)+0.5];

colors_ag = {'b','r','g'};

% Decimazione per velocizzare l'animazione
% ogni step_anim passi di simulazione = un frame
step_anim = 100;

for k = 1:step_anim:Nt

    clf;
    hold on;
    axis(ax_lim);
    axis equal;
    grid on;
    xlabel('x [m]'); ylabel('y [m]');
    title(sprintf('Formazione — t = %.2f s', t(k)));

    %--- Traiettoria percorsa fino a k (scia) ---
    for i = 1:3
        plot(squeeze(Y(1,i,1:k)), squeeze(Y(2,i,1:k)), ...
             '--', 'Color', colors_ag{i}, ...
             'LineWidth', 0.8);
    end

    % --- Traiettoria reference completa (grigio sfondo) ---
    plot(r(1,:), r(2,:), 'k:', 'LineWidth', 1.0);

    % --- Posizione reference al passo k ---
    plot(r(1,k), r(2,k), 'k+', ...
         'MarkerSize', 12, 'LineWidth', 2);

    % --- Triangolo desiderato al passo k ---
    x_des = arrayfun(@(i) r(1,k) + h{i}(1) - h1(1), 1:3);
    y_des = arrayfun(@(i) r(2,k) + h{i}(2) - h1(2), 1:3);
    fill([x_des x_des(1)], [y_des y_des(1)], ...
         'k', 'FaceAlpha', 0.05, 'EdgeColor', 'k', ...
         'LineStyle', '--', 'LineWidth', 1.0);

    % --- Triangolo reale al passo k ---
    x_real = squeeze(Y(1,:,k));
    y_real = squeeze(Y(2,:,k));
    fill([x_real x_real(1)], [y_real y_real(1)], ...
         'c', 'FaceAlpha', 0.15, 'EdgeColor', 'c', 'LineWidth', 1.2);

    % --- Agenti come cerchi colorati ---
    for i = 1:3
        plot(Y(1,i,k), Y(2,i,k), 'o', ...
             'MarkerSize', 12, ...
             'MarkerFaceColor', colors_ag{i}, ...
             'MarkerEdgeColor', 'k', ...
             'LineWidth', 1.5);
        % etichetta agente
        text(Y(1,i,k)+0.05, Y(2,i,k)+0.05, ...
             sprintf('A%d', i), 'FontSize', 9, ...
             'Color', colors_ag{i}, 'FontWeight', 'bold');
    end

    % --- Formation error norm nel titolo ---
    fe_now = 0;
    for i = 1:3
        fe_now = fe_now + norm(FE(:,i,k));
    end
    title(sprintf('t = %.2f s  |  ||FE|| totale = %.4f m', t(k), fe_now));

    legend('Scia A1','Scia A2','Scia A3', ...
           'r(t) completo','r(t) attuale', ...
           'Triangolo desiderato','Triangolo reale', ...
           'A1','A2','A3', ...
           'Location','northeast', 'FontSize', 7);

    drawnow;
end


