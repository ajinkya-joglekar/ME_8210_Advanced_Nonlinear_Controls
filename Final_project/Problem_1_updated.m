% Hartman-Grobman for Hamiltonian dynamics
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% Hamiltonian optimal control problem 1
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
H_x = polynomial_basis(x,p);
H_x = subs(H_x,p,0);
H_z = polynomial_basis(x,p);
disp(H_z)


%% Functions

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
