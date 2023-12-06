 function xdot=examplesystem(t,x)

 


v1=(x(1)-x(2))*(1-x(1)^2-x(2)^2);
v2=-(x(1)+x(2))*(1-x(1)^2-x(2)^2);


xdot=[v1;v2];