# Robust MPC for Omnidirectional Mobile Robot Navigation

This project implements a robust Model Predictive Control (MPC) framework for an
omnidirectional mecanum-wheel robot using a spatial-domain (s-domain) formulation.
The controller performs stable trajectory tracking and smooth obstacle avoidance
under disturbances.

> This work follows the spatial-domain MPC formulation proposed in [1], which resolves prediction inconsistency issues of time-domain models for omnidirectional robots.

Final project for **ME 550 – Nonlinear Optimal Control**, University of Washington.

---

## System Overview
![System architecture for s-domain MPC pipeline](/img/system%20structure.png)

---

## My Contributions
My contributions include implementing:
- nonlinear dynamics model and s-domain formulation  
- linearization & discretization  
- QP-based MPC solver  
- Gauss–Markov disturbance model  
- reference-shaping obstacle avoidance  

---

## Features
- Nonlinear mecanum-wheel robot model  
- s-domain error-state formulation  
- Robust MPC with jerk inputs  
- Smooth reference-based obstacle avoidance  
- Full MATLAB simulation + visualization  

---

## How to Run

1. Clone the repository and open the folder in MATLAB.
2. Make sure all `.m` files are on the MATLAB path.
3. Run the main simulation script:
    ```matlab
    main
    ```
4. (Optional) Enable obstacle avoidance:
    ```matlab
    build_ref_with_obstacle
    ```
5. (Optional) Visualize the robot trajectory:
    ```matlab
    animate_robot_path
    ```

---

## Results
### MPC vs PID Simulation Results 
<p align="center">
  <img src="/img/MPC_vs_PID.gif" width="40%">
  <img src="/img/MPC__PID.jpg" width="55%">
</p>

<p align="center">
  <em>Left: Path tracking result. Right: Tracking error vs time.</em>
</p>

### Obstacle Avoidance
<p align="center">
    <img src="/img/MPC_obstacle_avoidance.gif" width="50%">
</p>
<p align="center">
  <em>Obstacle Avoidance by using Smooth Lateral Reference Shaping.</em>
</p>

---

## Reference
[1] Y. Han and Q. Zhu, “Robust Optimal Control of Omni-directional Mobile Robot using Model Predictive Control Method,” in 2019 Chinese Control Conference (CCC), July 2019, pp. 4679–4684. doi: 10.23919/ChiCC.2019.8865344.

