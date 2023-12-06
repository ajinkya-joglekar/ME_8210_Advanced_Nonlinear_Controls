 function xdot=examplesystem(t,x)

 


v1=-x(2);
v2=-x(1)+x(2)*(1 - x(1)^2+0.1*(x(1)^4));


xdot=[v1;v2];