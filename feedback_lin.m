clear all; clc; close all;

p_test  = [1.0;  2.0];
v_test  = [0.5; -0.3];
nu_test = [0.2; -0.1];

% agent 1
phi1 = @(p,v) [ -0.35*v(1) - 0.08*v(1)*abs(v(1));
                -0.45*v(2) - 0.10*v(2)*abs(v(2)) ];
B1   = @(p,v)  (1/1.2) * eye(2);
det_B1 = (1/1.2)^2;
assert(abs(det_B1) > 1e-10, 'Agent 1: B singular!');
u1_fl = @(p,v,nu) B1(p,v) \ (nu - phi1(p,v));

% agent 2
phi2 = @(p,v) [ -0.40*v(1); -0.55*v(2) ];
B2   = @(p,v)  diag([ 1/(1 + 0.25*sin(p(1))^2); 1/(1 + 0.30*cos(p(2))^2) ]);
det_B2 = det(B2(p_test, v_test));
assert(abs(det_B2) > 1e-10, 'Agent 2: B singular!');
u2_fl = @(p,v,nu) B2(p,v) \ (nu - phi2(p,v));

% agent 3
phi3 = @(p,v) [ -0.30*v(1) + 0.15*sin(p(2))*v(2);
                -0.35*v(2) + 0.15*cos(p(1))*v(1) ];
B3   = @(p,v) [ 1 + 0.20*cos(p(2))^2,  0.08*sin(p(1));
                0.08*sin(p(2)),          1 + 0.20*sin(p(1))^2 ];
det_B3 = det(B3(p_test, v_test));
assert(abs(det_B3) > 1e-10, 'Agent 3: B singular!');
u3_fl = @(p,v,nu) B3(p,v) \ (nu - phi3(p,v));

% summary
Phi = {phi1, phi2, phi3};
B   = {B1,   B2,   B3  };
U_fl = {u1_fl, u2_fl, u3_fl};