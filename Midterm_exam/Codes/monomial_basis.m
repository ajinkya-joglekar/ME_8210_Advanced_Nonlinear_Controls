function [Psi, DPsi] = monomial_basis(deg, dim)
%  [Psi, DPsi] = monomial_basis(deg, dim) returns a monomial basis function and its derivative 
%   There is not a 1 included in Psi, only the linear and higher order terms
%   deg = degree of monomial
%   dim = number of states
k = linspace(2, deg, deg-1);
d = dim;
x=sym('x',[d,1]);
assume(x,'real')

Psi = [x.'];
% Following line 13-20 generates combinations of states of degree deg
for i=1:size(k,2)
    m = nchoosek(k(i)+d-1,d-1); 
    dividers = [zeros(m,1),nchoosek((1:(k(i)+d-1))',d-1),ones(m,1)*(k(i)+d)]; 
    a = diff(dividers,1,2)-1;
    for i = 1:size(a,1)
        Psi = [Psi prod(x.' .^ a(i,:))];
    end
end
DPsi = jacobian(Psi,x);
Psi = Psi';
end

