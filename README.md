# Formation Control of Nonlinear Mobile Robots

Distributed **formation control** for a team of **three heterogeneous nonlinear planar robots**, built in **MATLAB**. Each robot is first made linear with **input-output feedback linearization**. A distributed controller then drives the team into a **triangular formation** that tracks a moving reference. Two controller designs are implemented and compared: **state feedback** and **frequency-domain loop shaping**.

> Final project (Part II) for the course *Modeling and Control of Cyber-Physical Systems*, Politecnico di Torino, A.Y. 2025–2026.
> 📄 Project specification: [`Project_PartII.pdf`](Project_PartII.pdf)

---

## Problem

Each agent *i* has state $x_i = (p_{x,i}, v_{x,i}, p_{y,i}, v_{y,i})$, measured output $y_i = p_i$ (its planar position), and nonlinear control-affine dynamics

$$\dot p_i = v_i, \qquad \dot v_i = \varphi_i(p_i, v_i) + B_i(p_i, v_i)\,u_i$$

The three robots differ on purpose:

| Agent | Nonlinearity |
|---|---|
| 1 | Anisotropic nonlinear (quadratic) drag |
| 2 | Position-dependent apparent mass |
| 3 | Nonlinear cross-coupling between the x and y axes |

## Approach

1. **Input-output feedback linearization**: for each agent, the vector relative degree is (2, 2), so there are no zero dynamics. The decoupling matrix $B_i$ is shown to be nonsingular, and the control law

   $$u_i = B_i^{-1}\left[\nu_i - \varphi_i(p_i, v_i)\right]$$

   reduces every robot to a double integrator, $\ddot y_i = \nu_i$.
2. **Distributed formation controller**, designed on the virtual inputs $\nu_i$ using the communication graph (Laplacian eigenvalues) to reach a triangular formation around a moving leader/reference:
   - **State feedback** (`sim_state_feedback.m`): constant or circular reference trajectories.
   - **Loop shaping** (`sim_loop_shaping.m`, `controllore_con_due_reti.m`): a compensator built from two compensating networks, tuned on the modal loops $\lambda_k P(s)$ with $P(s) = 1/s^2$. The design is then tested for robustness against:
     - common vs differential **output disturbances**
     - **sensor noise** (GPS-like, σ ≈ 1 cm)
     - optional **feedforward** of the reference
3. **Evaluation**: trajectories, position tracking, formation errors $\delta_i = (y_i - y_1) - (h_i - h_1)$, error norms in transient and steady state, and the physical control effort $u_i$.

## Tech stack

**MATLAB** with the Control System Toolbox (`tf`, Bode-based loop shaping), plus custom nonlinear simulation loops.

## Project structure

```
├── main.m                      # Agent definitions + feedback linearization
├── utils.m                     # Helpers: agent builder, IOFL law, ...
├── sim_state_feedback.m        # Formation control via state feedback
├── sim_loop_shaping.m          # Formation control via loop shaping (+ disturbances/noise)
├── controllore_con_due_reti.m  # Two-network compensator design on the modal loops
└── Project_PartII.pdf          # Project specification
```

## How to run

1. Open the folder in MATLAB.
2. Run `main.m` to define the agents and compute the linearizing laws.
3. Run `sim_state_feedback.m` or `sim_loop_shaping.m`. At the top of `sim_loop_shaping.m` you can change the reference type, disturbance type (`none` / `common` / `differential`), sensor noise and feedforward.
