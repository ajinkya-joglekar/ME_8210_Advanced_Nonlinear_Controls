%%
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%%
x = sym('x', 'real');
u = sym('u','real');
% Define system dynamics
f = @(x, u) x - x^3 + u;  % Dynamics: x_dot = x - x^3 + u
g = @(x) 1;  % Control input matrix, assuming a single control input

% Define cost function
q = @(x,u) x^2 +u^2;  % Cost function: q(x) = x^2
psi = @(x) 0;  % Terminal cost, assuming psi(x(T)) = 0


% Define constants
C = 1;  % Control cost, assuming a single control input

% Define state and co-state variables
syms x p real;

% Define Hamiltonian
% H = p*f(x, -inv(C)*g(x)) + 0.5*(p*g(x)*inv(C)*g(x)*p + q(x));
H = q(x,u) + p*f(x,u);

% Solve for u*
du = diff(H,u);
solu_u = solve(du,u);


% Compute partial derivatives of the Hamiltonian
dH_dx = diff(H, x);
dH_dp = diff(H, p);

% Define matrices A, B, and Q
A = -3*x^2 + 1;
B = g(x);
Q = 0.5;

% Compute matrix E
E = A - B*inv(C)*B' - Q - A';

% Define the Hamiltonian dynamical system
Hamiltonian_system = [dH_dx -dH_dp; -E -dH_dx];

% Time parameters
T = 1;

% Initial conditions
x0 = 1;  % Initial state
pT = 0;  % Terminal co-state

% Solve the Hamiltonian system
[t, sol] = ode45(@(t, z) Hamiltonian_dynamics(t, z, A, B, Q, C, g), [0 T], [x0; pT]);

% Display results
disp('Optimal State and Co-state Trajectories:');
disp([t, sol(:, 1), sol(:, 2)]);



%% Redo


syms x1 x2 p1 p2 u;
Dx1 = x2;
Dx2 = -x2 + u;
% Cost function inside the integral
syms g;
g = 0.5*u^2;
% Hamiltonian
syms p1 p2 H;
H = g + p1*Dx1 + p2*Dx2;
% Costate equations
Dp1 = -diff(H,x1);
Dp2 = -diff(H,x2);
% solve for control u
du = diff(H,u);
sol_u = solve(du, u);

% Substitute u to state equations
Dx2 = subs(Dx2, u, sol_u);
% convert symbolic objects to strings for using ’dsolve’
eq1 = strcat('Dx1=',char(Dx1));
eq2 = strcat('Dx2=',char(Dx2));
eq3 = strcat('Dp1=',char(Dp1));
eq4 = strcat('Dp2=',char(Dp2));
sol_h = dsolve(eq1,eq2,eq3,eq4);

% case a: (a) x1(0)=x2(0)=0; x1(2) = 5; x2(2) = 2;
conA1 = 'x1(0) = 0';
conA2 = 'x2(0) = 0';
conA3 = 'x1(2) = 5';
conA4 = 'x2(2) = 2';
sol_a = dsolve(eq1,eq2,eq3,eq4,conA1,conA2,conA3,conA4);

% Extract symbolic solutions for x1 and x2
sol_x1_a = sol_a.x1;
sol_x2_a = sol_a.x2;

% Create a time vector
t = linspace(0, 2, 100);

% Evaluate symbolic solutions for x1 and x2 at different time points
x1_values_a = double(subs(sol_x1_a, 't', t));
x2_values_a = double(subs(sol_x2_a, 't', t));

% Plot x1 and x2
figure;
subplot(2,1,1);
plot(t, x1_values_a);
xlabel('Time');
ylabel('x1');
title('Plot of x1 over Time');

subplot(2,1,2);
plot(t, x2_values_a);
xlabel('Time');
ylabel('x2');
title('Plot of x2 over Time');


%% Functions

% function dzdt = Hamiltonian_dynamics(t, z, A, B, Q, C, g)
%     x = z(1);
%     p = z(2);
% 
%     % Compute matrix E
%     E = A - B*inv(C)*B' - Q - A';
% 
%     % Compute matrix Fn
%     Fn = [diff_H_diff_p(x, p, C, g); -E*x - diff_H_diff_x(x, p, A, B, Q, C, g)];
% 
%     % Hamiltonian system equations
%     dzdt = Fn;
% 
% end
% 
% function result = diff_H_diff_p(x, p, C, g)
%     % Compute partial derivatives of the Hamiltonian with respect to p
%     result = g(x)'*inv(C)*g(x)*p;
% end
% 
% function result = diff_H_diff_x(x, p, A, B, Q, C, g)
%     % Compute partial derivatives of the Hamiltonian with respect to x
%     result = -p*(-3*x^2 + 1) - 0.5*g(x)'*inv(C)*g(x)*p;
% end


%%
% Define system dynamics
f = @(t, z) [z(2); -z(1) - z(2) * (1 - (1 + 2 * sin(z(1)))^2) + (1 + 2 * sin(z(1))) * z(3)];

% Define cost function
L = @(u) u^2;

% Define Hamiltonian
H = @(x, p, u) L(u) + p(1) * f(0, [x; p(2); 0]) + p(2) * f(0, [x; p(2); 0]);

% Optimal control condition
du_optimal = @(p) -p(3) / (1 + 2 * sin(p(1)));

% Set up the ODE system for Hamiltonian dynamics
ode_system = @(t, z) [f(t, z(1:2)); -jacobian(H(z(1), z(3:4), du_optimal(z(3:4))))];

% Time vector
tspan = [0 10]; % adjust as needed

% Initial conditions
z0 = [0; 0; 0; 0]; % [x1; x2; p1; p2]

% Solve the ODE system
[t_values, z_values] = ode45(ode_system, tspan, z0);

% Extract solutions
x1_values_optimal = z_values(:, 1);
x2_values_optimal = z_values(:, 2);
u_values_optimal = arrayfun(du_optimal, z_values(:, 3:4));

% Plot the results
figure;

subplot(3, 1, 1);
plot(t_values, x1_values_optimal);
xlabel('Time');
ylabel('x1');
title('State x1');

subplot(3, 1, 2);
plot(t_values, x2_values_optimal);
xlabel('Time');
ylabel('x2');
title('State x2');

subplot(3, 1, 3);
plot(t_values, u_values_optimal);
xlabel('Time');
ylabel('u');
title('Optimal Control input u');
