%% =========================================================
%  FEEDBACK_LINEARIZATION.m
%  =========================================================
%  For each of the 3 agents:
%   1. Defines the nonlinear model  v_dot = phi(p,v) + B(p,v)*u
%   2. Shows the FL derivation step by step (demonstrative)
%   3. Verifies that  phi + B*u_fl = nu  (numerical check)
%
%  OUTPUT: struct array  agents(i)  with fields
%    .phi   — function handle  phi(p,v)     -> 2x1 drift vector
%    .B     — function handle  B(p,v)       -> 2x2 decoupling matrix
%    .u_fl  — function handle  u_fl(p,v,nu) -> 2x1 physical input
%    .name  — string description of the agent
%
%  These variables stay in the workspace and are used directly
%  by sim_state_feedback.m and sim_loop_shaping.m
%  (run this file first, then run either simulation)
%  =========================================================

clear; clc;

%% =========================================================
%  TEST POINT (same for all agents)
% ==========================================================
p_test  = [1.0;  2.0];
v_test  = [0.5; -0.3];
nu_test = [0.2; -0.1];


%% =========================================================
%  AGENT 1 — Anisotropic nonlinear drag
%
%  px_dot = vx
%  vx_dot = -0.35*vx - 0.08*vx*|vx|  +  (1/1.2)*ux
%  py_dot = vy
%  vy_dot = -0.45*vy - 0.10*vy*|vy|  +  (1/1.2)*uy
% ==========================================================
fprintf('============================================================\n');
fprintf('  AGENT 1 — Anisotropic nonlinear drag\n');
fprintf('============================================================\n');

% Step 1: identify phi1 and B1
% phi1 = everything that does NOT multiply u
% B1   = matrix that multiplies u
%
%   phi1(p,v) = [-0.35*vx - 0.08*vx*|vx|]
%               [-0.45*vy - 0.10*vy*|vy|]
%
%   B1(p,v) = (1/1.2) * I_2   (constant, no state dependence)

agents(1).name = 'Anisotropic nonlinear drag';

agents(1).phi = @(p,v) [ -0.35*v(1) - 0.08*v(1)*abs(v(1));
                          -0.45*v(2) - 0.10*v(2)*abs(v(2)) ];

agents(1).B   = @(p,v)  (1/1.2) * eye(2);

% Step 2: nonsingularity
% det(B1) = (1/1.2)^2 = 0.6944  > 0  always  (constant matrix)
det_B1 = (1/1.2)^2;
fprintf('phi1(p_test, v_test) = [%.4f; %.4f]\n', ...
    agents(1).phi(p_test,v_test));
fprintf('B1 = (1/1.2)*I  =>  det(B1) = %.4f  (constant, always > 0)\n', det_B1);

% Step 3: feedback-linearizing law
%   u1 = B1^{-1} * (nu1 - phi1)  =>  y_ddot = nu1
agents(1).u_fl = @(p,v,nu) agents(1).B(p,v) \ (nu - agents(1).phi(p,v));

% Step 4: numerical verification
y_ddot_1 = agents(1).phi(p_test,v_test) ...
          + agents(1).B(p_test,v_test) * agents(1).u_fl(p_test,v_test,nu_test);
fprintf('Verification:  ||y_ddot - nu|| = %.2e  (should be ~0)\n\n', ...
    norm(y_ddot_1 - nu_test));


%% =========================================================
%  AGENT 2 — Position-dependent apparent mass
%
%  px_dot = vx
%  vx_dot = -0.40*vx  +  1/(1+0.25*sin^2(px)) * ux
%  py_dot = vy
%  vy_dot = -0.55*vy  +  1/(1+0.30*cos^2(py)) * uy
% ==========================================================
fprintf('============================================================\n');
fprintf('  AGENT 2 — Position-dependent apparent mass\n');
fprintf('============================================================\n');

% Step 1: identify phi2 and B2
%
%   phi2(p,v) = [-0.40*vx]
%               [-0.55*vy]
%
%   B2(p,v) = diag( 1/(1+0.25*sin^2(px)),  1/(1+0.30*cos^2(py)) )
%             diagonal, both entries always in (0, 1]

agents(2).name = 'Position-dependent apparent mass';

agents(2).phi = @(p,v) [ -0.40*v(1);
                          -0.55*v(2) ];

agents(2).B   = @(p,v)  diag([ 1 / (1 + 0.25*sin(p(1))^2);
                                1 / (1 + 0.30*cos(p(2))^2) ]);

% Step 2: nonsingularity
% b11 = 1/(1+0.25*sin^2(px)) >= 1/1.25 = 0.80 > 0
% b22 = 1/(1+0.30*cos^2(py)) >= 1/1.30 = 0.77 > 0
% det(B2) = b11*b22 >= 0.80*0.77 = 0.616 > 0  always
B2_test = agents(2).B(p_test, v_test);
fprintf('phi2(p_test, v_test) = [%.4f; %.4f]\n', ...
    agents(2).phi(p_test,v_test));
fprintf('B2 = diag([%.4f, %.4f])\n', B2_test(1,1), B2_test(2,2));
fprintf('b11 >= 0.80,  b22 >= 0.77  =>  det(B2) >= 0.616 > 0  always\n');
fprintf('det(B2) at test point = %.6f\n', det(B2_test));

% Step 3: feedback-linearizing law
agents(2).u_fl = @(p,v,nu) agents(2).B(p,v) \ (nu - agents(2).phi(p,v));

% Step 4: numerical verification
y_ddot_2 = agents(2).phi(p_test,v_test) ...
          + agents(2).B(p_test,v_test) * agents(2).u_fl(p_test,v_test,nu_test);
fprintf('Verification:  ||y_ddot - nu|| = %.2e  (should be ~0)\n\n', ...
    norm(y_ddot_2 - nu_test));


%% =========================================================
%  AGENT 3 — Nonlinear coupling
%
%  px_dot = vx
%  vx_dot = -0.30*vx + 0.15*sin(py)*vy
%           + (1+0.20*cos^2(py))*ux  +  0.08*sin(px)*uy
%  py_dot = vy
%  vy_dot = -0.35*vy + 0.15*cos(px)*vx
%           +  0.08*sin(py)*ux  +  (1+0.20*sin^2(px))*uy
% ==========================================================
fprintf('============================================================\n');
fprintf('  AGENT 3 — Nonlinear coupling\n');
fprintf('============================================================\n');

% Step 1: identify phi3 and B3
%
%   phi3(p,v) = [-0.30*vx + 0.15*sin(py)*vy]
%               [-0.35*vy + 0.15*cos(px)*vx]
%
%   B3(p,v) = [1+0.20*cos^2(py),   0.08*sin(px)]
%             [0.08*sin(py),         1+0.20*sin^2(px)]

agents(3).name = 'Nonlinear coupling';

agents(3).phi = @(p,v) [ -0.30*v(1) + 0.15*sin(p(2))*v(2);
                          -0.35*v(2) + 0.15*cos(p(1))*v(1) ];

agents(3).B   = @(p,v)  [ 1 + 0.20*cos(p(2))^2,   0.08*sin(p(1));
                           0.08*sin(p(2)),            1 + 0.20*sin(p(1))^2 ];

% Step 2: nonsingularity
% Let: a = 1+0.20*cos^2(py) in [1.00, 1.20]
%      d = 1+0.20*sin^2(px) in [1.00, 1.20]
%      b = 0.08*sin(px)     in [-0.08, 0.08]
%      c = 0.08*sin(py)     in [-0.08, 0.08]
%
% det(B3) = a*d - b*c >= 1.0 - 0.08^2 = 0.9936 > 0  always
B3_test = agents(3).B(p_test, v_test);
fprintf('phi3(p_test, v_test) = [%.4f; %.4f]\n', ...
    agents(3).phi(p_test,v_test));
fprintf('B3 at test point =\n');
fprintf('  [%.4f   %.4f]\n', B3_test(1,1), B3_test(1,2));
fprintf('  [%.4f   %.4f]\n', B3_test(2,1), B3_test(2,2));
fprintf('det(B3) = a*d - b*c >= 1 - 0.08^2 = %.4f > 0  always\n', 1-0.08^2);
fprintf('det(B3) at test point = %.6f\n', det(B3_test));

% Step 3: feedback-linearizing law
agents(3).u_fl = @(p,v,nu) agents(3).B(p,v) \ (nu - agents(3).phi(p,v));

% Step 4: numerical verification
y_ddot_3 = agents(3).phi(p_test,v_test) ...
          + agents(3).B(p_test,v_test) * agents(3).u_fl(p_test,v_test,nu_test);
fprintf('Verification:  ||y_ddot - nu|| = %.2e  (should be ~0)\n\n', ...
    norm(y_ddot_3 - nu_test));


%% =========================================================
%  RELATIVE DEGREE — same for all agents
% ==========================================================
fprintf('============================================================\n');
fprintf('  RELATIVE DEGREE (common to all agents)\n');
fprintf('============================================================\n');
fprintf('y      = p            =>  output is position\n');
fprintf('y_dot  = v            =>  input u does NOT appear\n');
fprintf('y_ddot = phi(p,v) + B(p,v)*u  =>  input appears here\n');
fprintf('Vector relative degree: (r_x, r_y) = (2, 2)\n');
fprintf('Sum = 4 = dim(state)  =>  complete linearization, no zero dynamics\n\n');
fprintf('agents struct is ready in the workspace.\n');
fprintf('Run sim_state_feedback.m or sim_loop_shaping.m\n');