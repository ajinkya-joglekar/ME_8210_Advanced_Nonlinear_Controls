% Hartman-Grobman for Hamiltonian dynamics
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% Test code 
% x = sym('x', 'real');
% u = sym('u','real');
% % Define system dynamics
% f = @(x, u) x - x^3 + u;  % Dynamics: x_dot = x - x^3 + u
% g = @(x) 1;  % Control input matrix, assuming a single control input
% 
% % Define cost function
% q = @(x,u) x^2 +u^2;  % Cost function: q(x) = x^2
% psi = @(x) 0;  % Terminal cost, assuming psi(x(T)) = 0
% 
% 
% % Define constants
% C = 1;  % Control cost, assuming a single control input
% 
% % Define state and co-state variables
% syms x p real;
% 
% % Define Hamiltonian
% % H = p*f(x, -inv(C)*g(x)) + 0.5*(p*g(x)*inv(C)*g(x)*p + q(x));
% H = q(x,u) + p*f(x,u);
% 
% % Solve for u*
% du = diff(H,u);
% solu_u = solve(du,u);
% 
% 
% % Compute partial derivatives of the Hamiltonian
% dH_dx = diff(H, x);
% dH_dp = diff(H, p);
% 
% % Define matrices A, B, and Q
% A = -3*x^2 + 1;
% B = g(x);
% Q = 0.5;
% 
% % Compute matrix E
% E = A - B*inv(C)*B' - Q - A';
% 
% % Define the Hamiltonian dynamical system
% Hamiltonian_system = [dH_dx -dH_dp; -E -dH_dx];
% 
% % Time parameters
% T = 1;
% 
% % Initial conditions
% x0 = 1;  % Initial state
% pT = 0;  % Terminal co-state
% 
% % Solve the Hamiltonian system
% [t, sol] = ode45(@(t, z) Hamiltonian_dynamics(t, z, A, B, Q, C, g), [0 T], [x0; pT]);
% 
% % Display results
% disp('Optimal State and Co-state Trajectories:');
% disp([t, sol(:, 1), sol(:, 2)]);



%% Redo
syms x  p u T;
Dx = x - x^3 + u;
g = 1;
% Cost function inside the integral
syms q J;
q = x^2;
J = x^2 +u^2;
% Hamiltonian
syms  H;
H = J + p*Dx;
% Costate equations
Dp = -diff(H,x);

% solve for control u
du = diff(H,u);
sol_u = solve(du, u);

% Substitute u into the Hamiltonian 

H = subs(H,u,sol_u);

%Find Dh/dx and Dh/du
% This is same as the state equation. i.e -dH/dp which gives us the dynamics
x_dot = diff(H,p);
p_dot = -diff(H,x);

% Define matrices A, B, and Q
A = -3*x^2 + 1;
B = eye(1);
C= eye(1);
Q = eye(1);

% According to equation 7 and 8, lets divide the Hamiltonian into linear
% and nonlinear parts
df_dx = diff(Dx,x);
dq_dx = diff(q,x);
d2q_dx2 = diff(dq_dx,x);

E = [subs(df_dx,x,0),-1;-0.5*subs(d2q_dx2,x,0),-subs(df_dx,x,0)']; %Linear part
F_n = [x_dot;p_dot]-[A,-B*inv(C)*B';-Q,-A']*[x;p];

% Convert the non-linear hamiltonian control into linear form
syms z;
z = [x',p']';
H_z = polynomial_basis(x,p);
disp(H_z)


%%


%% Traidtional approach for hamiltonian based on pontryagin principle
% convert symbolic objects to strings for using ’dsolve’
eq1 = strcat('Dx1=',char(x_dot));
% eq2 = strcat('Dx2=',char(del));
eq2 = strcat('Dp1=',char(p_dot));
% eq4 = strcat('Dp2=',char(Dp2));
sol_h = dsolve(eq1,eq2);


% Testing out the solution when x(5) is free
eq1b = subs(sol_h.x1,'t',0);
eq2b = subs(sol_h.p1,'t',2) == subs(sol_h.x1,'t',2);

sol_b = solve(eq1b,eq2b);

% case a: (a) x1(0)=x2(0)=0; x1(2) = 5; x2(2) = 2;
conA1 = 'x(0) = 0';
% conA2 = 'x2(0) = 0';
conA3 = 'x(2) = 5';
% conA4 = 'x2(2) = 2';
sol_a = dsolve(eq1,eq2,conA1,conA3);

% Extract symbolic solutions for x1 and x2
sol_x1 = sol_a.x1;
% sol_x2_a = sol_a.x2;

% Create a time vector
t = linspace(0, 2, 100);

% Evaluate symbolic solutions for x1 and x2 at different time points
x_values_a = double(subs(sol_x1, 't', t));
% x2_values_a = double(subs(sol_x2_a, 't', t));

% Plot x1 and x2
figure;
subplot(2,1,1);
plot(t, x_values_a);
xlabel('Time');
ylabel('x1');
title('Plot of x1 over Time');

% subplot(2,1,2);
% plot(t, x2_values_a);
% xlabel('Time');
% ylabel('x2');
% title('Plot of x2 over Time');


%% Functions

function dzdt = Hamiltonian_dynamics(t, z, A, B, Q, C, g)
    x = z(1);
    p = z(2);

    % Compute matrix E
    E = A - B*inv(C)*B' - Q - A';

    % Compute matrix Fn
    Fn = [diff_H_diff_p(x, p, C, g); -E*x - diff_H_diff_x(x, p, A, B, Q, C, g)];

    % Hamiltonian system equations
    dzdt = Fn;

end

function result = diff_H_diff_p(x, p, C, g)
    % Compute partial derivatives of the Hamiltonian with respect to p
    result = g(x)'*inv(C)*g(x)*p;
end

function result = diff_H_diff_x(x, p, A, B, Q, C, g)
    % Compute partial derivatives of the Hamiltonian with respect to x
    result = -p*(-3*x^2 + 1) - 0.5*g(x)'*inv(C)*g(x)*p;
end

function basis = polynomial_basis(x,p)
% syms x p;

% Define the maximum degree for the polynomial basis
max_degree = 2; % Adjust as needed

% Generate all possible combinations of monomials up to the specified degree
basis = [];
for degree_x = 0:max_degree
    for degree_y = 0:max_degree
        basis = [basis; x^degree_x * p^degree_y];
    end
end

% Display the generated polynomial basis
disp('Polynomial Basis:');
disp(basis);
end
