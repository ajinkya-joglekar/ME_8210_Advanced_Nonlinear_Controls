function [U, Phi,Phi_lin, DPhi] = PolyEigenfunc1(deg , n, N, Dom,  f)
% PolyEigenfunc  returns paramter and basis function and its derivative
%   to approximate eigenfunctions.
[Psi, ~] = monomial_basis(deg, n); % Call basis function file 
Psi = Psi(n+1:end); %remove linear part
Nbs = length(Psi); % length of only nonlinear monomials
x = sym('x',[n,1]);
DPsi = jacobian(Psi,x);
Psi = matlabFunction(Psi,'Vars',{x});
DPsi = matlabFunction(DPsi,'Vars',{x});
%
X = rand(n,N)*2*Dom(2) - Dom(2); % generate N random data points in Domain
E = double(subs(jacobian(f(0,x),x),x,zeros(n,1))); % Linearization at origin
Ez = E*X; %Linearization --> A in class notation
G_z = Psi(X); % Basis functions
F_z = f(0, X); % 
Df_z = zeros(Nbs,1);
for i=1:size(X,2)
    Df_z(:,i) = DPsi(X(:,i))*f(0, X(:,i));
end
% Use YALMIP solver
Uy = sdpvar(n,Nbs); % optimization variable
Objective =  norm(Uy*Df_z + F_z - E*Uy*G_z - Ez, 'fro'); % convex optimization objective
Constraints = [];
% opt = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
opt = sdpsettings('solver','sdpt3','verbose',0,'cachesolvers',1);
optimize(Constraints,Objective,opt);
U = value(Uy);
[~,D1,W1] = eig(E);
[W, ~] = cdf2rdf(W1, D1);
Phi =  @(x) W'*x +  W'*(U*Psi(x)); % Nonlinear eigenfunction in function handle form
Phi_lin =  W'*x ; % Linear part of eigenfunction
DPhi = @(x) W'  +  W'*(U*DPsi(x)); % Jacobian of nonlinear eigenfunction
end


