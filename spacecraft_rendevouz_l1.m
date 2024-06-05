%updated code 
clc;
clear; 
cvx_clear;

%mu - gravitational parameter (m3/s2)
mu = 3.986e14

%radius of Earth (m)
r_e = 6.371e6

%altitude of the target orbit (m)
%altitude = 422137

%a - radius of the target body's circular orbit (m)
%a = r_e + altitude

%sma from ppt
a = 6.85635e6

%constant mean motion of target n = sqrt(mu/a^3)
%units: rad/s
n = sqrt(mu/a^3)

%number of states
nx = 6

%number of control inputs
nu = 3

%number of knot points
%this is to plan 2 minitues worth of controls
N = 121

%timestep (seconds)
%need to wait at least 60 seconds between 
%successive thruster actuations. Encoded with 
%this timestep
dt = 1

%time history
thist = linspace(1, N*dt, N)

%Target Spacecraft State
x0_target = zeros(1, nx)'

%Chaser Spacecraft Initial State
%position - m
%velocity - m/s

%arbitrary state
%x0_chaser = [15.0, 20.0, 15.0, 10.0, 10.0, 0.0]'

%using the true data from the ppt
x0_chaser = [-0.2906951578333974,
             -258.73424997857484,
             1984.0093925290266,
             -2.2132577166418366,
             4.1627960172263556e-5,
             -0.00031920885703584645]

%Clohessy Wiltshire Equations in statespace form

A = zeros(nx, nx)
A(1:3, 4:6) = eye(3)
A(4, 1) = 3*n^2
A(6, 3) = -n^2
A(4, 5) = 2*n
A(5, 4) = -2*n

%mass of the satellite (kg) from data
m = 5.22

B = 1/m*[zeros(3,3); eye(3)]

%Spacecraft Linear Dynamics
%CW equations
function xdot = spacecraft_dynamics(x, u)

    xdot = A*x + B*u

end

%Discretize the dynamics model 
%with matrix exponential

H = expm(dt*[A B; zeros(nu, nx+nu)])

%Discrete Dynamics Matrices
Ad  = H(1:nx, 1:nx)
Bd = H(1:nx, (nx+1):end)

%initial position in target center
%coordinate frame
x_initial = x0_chaser

%Goal is the target position
x_goal = x0_target

cvx_begin

   %solver 
   cvx_solver sedumi

   %precision for solver tolerance
   cvx_precision('low')

   %maximum iterations
   %cvx_solver_settings('maxit', 1000);

    %state trajectory 
   variable X(nx, N)

   %Controls
   variable U(nu, N-1)

   %minimize the L1 norm to get bang bang controls
   minimize(norm(U(:),1) + norm(X(:), 2))

   subject to

       %initial state
       X(:,1) == x_initial
    
       %Goal constraint
       %X(:, N) == x_goal
    
       %Thrust Limit Constraint
       %currently infeasible. removing this 
       %it solves when removing this constraint and 
       %we get impulsive controls
       for k=1:(N-1)
           %maximum 4.6 mm/s for 1 second burn
           norm(U(:,k), 1) <= 4.6e-3

       end
    
       %Dynamics Constraints
    
       for k=1:N-1
           X(:,k+1) == Ad*X(:,k) + Bd*U(:,k)
       end

       %successive controls must be 60 seconds apart
       for k = 2:60
            U(:, k) == zeros(3,1)
       end

       for k = 62:N-1
            U(:, k) == zeros(3,1)
       end
       
cvx_end






figure('Name', 'My Plot', 'NumberTitle', 'off', 'Visible', 'on');

% Plot the data
plot3(X(1,:), X(2,:), X(3,:), '-b', 'LineWidth', 1.5);
hold on
scatter3(X(1,1), X(2,1), X(3,1), "green", "filled")
scatter3(X(1,N), X(2,N), X(3,N), "red", "filled")
hold off

% Add title and labels
title('Chaser Trajectory');
xlabel('X (m)');
ylabel('Y (m)');
zlabel("Z (m)")

legend('Trajectory', 'Start Point', 'End Point', 'Location', 'best');

% Create a new figure window
figure('Name', 'Multiple Subplots', 'NumberTitle', 'off', 'Visible', 'on', 'WindowStyle', 'normal');

% First subplot
subplot(3, 1, 1); 
plot(thist(1, 1:N-1), U(1,:), '-b', 'LineWidth', 1.5); 
title('X-Control');
xlabel('Time (s)');
ylabel('Delta v (m/s)');

% Second subplot
subplot(3, 1, 2);
plot(thist(1, 1:N-1), U(2,:), '-r', 'LineWidth', 1.5); 
title('Y-Control');
xlabel('Time (s)');
ylabel('Delta v (m/s)');

% Third subplot
subplot(3, 1, 3); 
plot(thist(1, 1:N-1), U(3,:), '-g', 'LineWidth', 1.5); 
title('Z-Control');
xlabel('Time (s)');
ylabel('Delta v (m/s)');