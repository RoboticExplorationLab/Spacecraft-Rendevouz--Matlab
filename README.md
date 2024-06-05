# spacecraft-rendevouz-matlab

This repo is an example of a chaser spacecraft performing a fuel-optimal rendevouz manuever toward a target spacecraft. The implementation is done in MATLAB R2024a. It is currently solving with the SeDuMi solver (free solver included in CVX).  

## Package Requirements
- Download the CVX package for MATLAB: https://github.com/cvxr/cvx
- Follow the installation instructions in the following page: https://cvxr.com/cvx/download/

## Implementation 
- Model: Clohessy Wiltshire (CW) Equations
- Objective: Minimize the L1 norm of the controls. (assuming impulsive burns)
- Constraints: 
    - Dynamics Constraint (Discrete CW Equations)
    - Thrust limit
    - Initial Chaser Position
    - Goal Chaser Position

## Variables that can tuned for a specific appliation: 
- Target Satellite Orbit altitude
- Chaser initial position 
- Mass
- Timestep, amount of knot points, total time to execute the maneuver
- Solver settings


