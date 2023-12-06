% Hartman-Grobman for Hamiltonian dynamics
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])


%%

hartman_grobman()
%%
function hartman_grobman()
    % System dynamics
    dx = @(x, u) x - x^3 + u;

    % Cost function
    J = @(x, u) x^2 + u^2; % Since C is identity matrix in this case

    % Hamiltonian
    hamiltonian = @(x, lambda, u) lambda * dx(x, u) + 0.5*J(x, u);

    % Adjoint equation
    adjoint_eq = @(t, lambda, x, u) -2*x - (1 - 3*x^2) * lambda;

    % Time parameters
    T = 1;  % Final time
    tspan = [T, 0];  % Integrate backward in time

    % Initial condition for the adjoint variable
    lambda_terminal = 0;

    % Integrate the adjoint equation backward in time
    [t_adjoint, lambda_adjoint] = ode45(@(t, lambda) adjoint_eq(t, lambda, x(t), u(t)), tspan, lambda_terminal);

    % Display or use Lagrange multiplier
    disp('Lagrange Multiplier at Initial Time:');
    disp(lambda_adjoint(end));
end

% Define the state and control trajectories using a simple example
function x = x(t)
    % Define the state trajectory
    x = 0;
end

function u = u(t)
    % Define the control trajectory
    u =1;
end

