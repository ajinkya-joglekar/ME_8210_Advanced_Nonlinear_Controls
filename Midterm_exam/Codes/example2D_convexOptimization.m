%% Code Files 
% This code uses two function files 
% 1) PolyEigenfunc1: This solves the convex optimization 
% 2) monomial_basis: This generates monomial basis of degree deg...
                     %...We can add also add some other polynomial basis at the end of this file eg. sin(x)    
% Toolbox required: YALMIP, Gurobi (or SDP3) 

clear all; close all; clc;
set(0,'DefaultLineLineWidth',2.5) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% Duffing oscillator Dynamics
lambda1=-2.5; lambda2=-1.0;
f = @(t, x) [lambda1*x(1,:); lambda2*(x(2,:)-x(1,:).^2)];
n =2;% # of state variable
Dom = [-2 2];  % Domain for eigenfunction 
g = @(t,x) [1; 0]; 
fN_u= @(t, x,u) f(t, x)+ g(t, x)*u; 

%% linear system at origin
x = sym('x',[2;1]); % symbolic state variable
A = double(subs(jacobian(f(0,x),x),x,[0;0])); % system matrix jacobian
u = sym('u',[1;1]);
B = double(subs(jacobian(fN_u(0,x,u),u),[x;u],[0;0;0])); % input matrix
[~,D1,W1] = eig(A);
[W, D] = cdf2rdf(W1, D1); % required in case of complex eigenvalues for realification

%% Eigenfunciton computation
% call PolyEigenfunc1 file where deg: order of monomial basis 
N = 1e3; % # of data points in domain for eigenfunction computation
deg =3; % order poly
[U, Phi, DPhi] = PolyEigenfunc1(deg , n, N, Dom,  f);
% U: parameterization coefficient 
% Phi: Nonlinear eigenfunction 
% DPhi: Jacobian of Phi

%% plots Eigenfunction
% create a grid of data points for ploting in the Domain
xvec = linspace(Dom(1),Dom(2),300); 
yvec = linspace(Dom(1),Dom(2),300);
[XX,YY] = ndgrid(xvec,yvec);
% zero initialization for phi1, phi2
phi1 = zeros(size(XX)); %NL parts of eigenfunctions
phi2 = zeros(size(XX));
% Obtain the eigenfunction value for the grid points
for i = 1:length(xvec)
    for j = 1:length(yvec)
        PHI = Phi([xvec(i);yvec(j)]);
        phi1(i,j) = PHI(1);
        phi2(i,j) = PHI(2);
    end
end
%% Plots
% Eigenfunction 1
figure(1)
p1 = pcolor(XX,YY,phi1);
colormap jet
colorbar
set(gca,'ColorScale','linear') % For log scale , comment for normal scale 
set(p1,'Edgecolor','none')
title('$\phi_1(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')

% Eigenfunction 2
figure(2)
p2 = pcolor(XX,YY,phi2);
colormap jet
colorbar
set(gca,'ColorScale','linear')
set(p2,'Edgecolor','none')
title('$\phi_2(x)$','Interpreter','latex')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'fontsize',20,'FontName','Times New Roman')
