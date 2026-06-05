%% =========================================================
%  FEEDBACK LINEARIZATION - 3 Nonlinear Planar Robots
%  =========================================================
%  For each agent:
%   1. Define the nonlinear model:  v_dot = phi(p,v) + B(p,v)*u
%   2. Compute the feedback-linearizing input:  u = B^{-1}*(nu - phi)
%   3. Verify that the closed-loop output satisfies:  y_ddot = nu
%  =========================================================

clear; clc;

%% =========================================================
%  AGENT 1 - Anisotropic nonlinear drag
%
%  px_dot = vx
%  vx_dot = -0.35*vx - 0.08*vx*|vx| + (1/1.2)*ux
%  py_dot = vy
%  vy_dot = -0.45*vy - 0.10*vy*|vy| + (1/1.2)*uy
%
%  Output: y1 = [px; py]
% ==========================================================

% --- Define a test state ---
p1 = [1.0; 2.0];    % [px; py]
v1 = [0.5; -0.3];  % [vx; vy]

% --- Identify phi1 and B1 ---
% v_dot = phi1 + B1*u  =>  read directly from the equations above
phi1 = [ -0.35*v1(1) - 0.08*v1(1)*abs(v1(1));   % drift x
         -0.45*v1(2) - 0.10*v1(2)*abs(v1(2)) ];  % drift y

B1 = (1/1.2) * eye(2);   % constant diagonal matrix

% --- Check nonsingularity ---
det_B1 = det(B1);
fprintf('=== AGENT 1 ===\n');
fprintf('phi1      = [%.4f; %.4f]\n', phi1(1), phi1(2));
fprintf('B1        = (1/1.2)*I  =>  det(B1) = %.6f  (always > 0)\n', det_B1);

% --- Choose a virtual input nu1 (what we WANT y_ddot to be) ---
nu1 = [0.2; -0.1];   % desired acceleration

% --- Compute the physical input u1 via feedback linearization ---
u1 = B1 \ (nu1 - phi1);   % equivalent to inv(B1)*(nu1 - phi1)

% --- Verify:  phi1 + B1*u1  should equal  nu1 ---
y_ddot_1 = phi1 + B1*u1;
fprintf('nu1       = [%.4f; %.4f]  (desired)\n',   nu1(1), nu1(2));
fprintf('y_ddot_1  = [%.4f; %.4f]  (achieved)\n',  y_ddot_1(1), y_ddot_1(2));
fprintf('error     = %.2e\n\n', norm(y_ddot_1 - nu1));


%% =========================================================
%  AGENT 2 - Position-dependent apparent mass
%
%  px_dot = vx
%  vx_dot = -0.40*vx + 1/(1 + 0.25*sin^2(px)) * ux
%  py_dot = vy
%  vy_dot = -0.55*vy + 1/(1 + 0.30*cos^2(py)) * uy
%
%  Output: y2 = [px; py]
% ==========================================================

% --- Define a test state ---
p2 = [0.8; 1.5];
v2 = [0.3;  0.4];

% --- Identify phi2 and B2 ---
phi2 = [ -0.40*v2(1);
         -0.55*v2(2) ];

b11 = 1 / (1 + 0.25*sin(p2(1))^2);   % always in (0, 1]
b22 = 1 / (1 + 0.30*cos(p2(2))^2);   % always in (0, 1]
B2  = diag([b11, b22]);

% --- Check nonsingularity ---
% b11 = 1/(1 + 0.25*sin^2(px)) >= 1/(1+0.25) = 0.8  > 0
% b22 = 1/(1 + 0.30*cos^2(py)) >= 1/(1+0.30) = 0.77 > 0
% => det(B2) = b11*b22 > 0  always
det_B2 = det(B2);
fprintf('=== AGENT 2 ===\n');
fprintf('phi2      = [%.4f; %.4f]\n', phi2(1), phi2(2));
fprintf('B2        = diag([%.4f, %.4f])  =>  det(B2) = %.6f  (always > 0)\n', b11, b22, det_B2);

% --- Choose a virtual input nu2 ---
nu2 = [0.15; 0.05];

% --- Compute the physical input u2 ---
u2 = B2 \ (nu2 - phi2);

% --- Verify ---
y_ddot_2 = phi2 + B2*u2;
fprintf('nu2       = [%.4f; %.4f]  (desired)\n',   nu2(1), nu2(2));
fprintf('y_ddot_2  = [%.4f; %.4f]  (achieved)\n',  y_ddot_2(1), y_ddot_2(2));
fprintf('error     = %.2e\n\n', norm(y_ddot_2 - nu2));


%% =========================================================
%  AGENT 3 - Nonlinear coupling
%
%  px_dot = vx
%  vx_dot = -0.30*vx + 0.15*sin(py)*vy
%           + (1 + 0.20*cos^2(py))*ux  +  0.08*sin(px)*uy
%  py_dot = vy
%  vy_dot = -0.35*vy + 0.15*cos(px)*vx
%           +  0.08*sin(py)*ux  +  (1 + 0.20*sin^2(px))*uy
%
%  Output: y3 = [px; py]
% ==========================================================

% --- Define a test state ---
p3 = [0.5; 1.0];
v3 = [0.2; -0.5];

% --- Identify phi3 and B3 ---
% phi3 = terms that do NOT multiply u
phi3 = [ -0.30*v3(1) + 0.15*sin(p3(2))*v3(2);
         -0.35*v3(2) + 0.15*cos(p3(1))*v3(1) ];

% B3 = matrix that multiplies u  (read column by column from equations)
%   vx_dot:  B3(1,1)*ux + B3(1,2)*uy
%   vy_dot:  B3(2,1)*ux + B3(2,2)*uy
B3 = [ 1 + 0.20*cos(p3(2))^2,   0.08*sin(p3(1));
       0.08*sin(p3(2)),           1 + 0.20*sin(p3(1))^2 ];

% --- Check nonsingularity ---
% Let:  a = 1 + 0.20*cos^2(py) in [1.00, 1.20]
%       d = 1 + 0.20*sin^2(px) in [1.00, 1.20]
%       b = 0.08*sin(px)        in [-0.08, 0.08]
%       c = 0.08*sin(py)        in [-0.08, 0.08]
%
% det(B3) = a*d - b*c
%         >= 1*1 - 0.08*0.08 = 1 - 0.0064 = 0.9936 > 0  always
det_B3 = det(B3);
fprintf('=== AGENT 3 ===\n');
fprintf('phi3      = [%.4f; %.4f]\n', phi3(1), phi3(2));
fprintf('B3        =\n'); disp(B3);
fprintf('det(B3)   = %.6f  (lower bound = %.4f  => always > 0)\n', det_B3, 1 - 0.08^2);

% --- Choose a virtual input nu3 ---
nu3 = [-0.05; 0.10];

% --- Compute the physical input u3 ---
u3 = B3 \ (nu3 - phi3);

% --- Verify ---
y_ddot_3 = phi3 + B3*u3;
fprintf('nu3       = [%.4f; %.4f]  (desired)\n',   nu3(1), nu3(2));
fprintf('y_ddot_3  = [%.4f; %.4f]  (achieved)\n',  y_ddot_3(1), y_ddot_3(2));
fprintf('error     = %.2e\n\n', norm(y_ddot_3 - nu3));


%% =========================================================
%  RELATIVE DEGREE SUMMARY
% ==========================================================
fprintf('=== RELATIVE DEGREE (same for all agents) ===\n');
fprintf('y   = p          =>  output = position\n');
fprintf('y_dot  = v       =>  input u does NOT appear\n');
fprintf('y_ddot = phi+B*u =>  input u appears here\n');
fprintf('Vector relative degree: (r_x, r_y) = (2, 2)\n');
fprintf('Sum = 4 = dim(state)  =>  no zero dynamics\n');