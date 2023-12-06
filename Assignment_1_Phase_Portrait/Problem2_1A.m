 function xdot=examplesystem(t,x)

 


v1=-x(1)+2*x(1)^3+x(2);
v2=-x(1)-x(2);


xdot=[v1;v2];