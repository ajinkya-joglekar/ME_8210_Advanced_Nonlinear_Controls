 function xdot=examplesystem(t,x)

 

%Duffing oscillator

a=0.1;

%v1=-x(2);
%v2=x(1)-(1-x(1)^2)*x(2);%(x(1)-x(1)^3)-0.5*x(2);
%v1= 2*x(1)^3-2*x(1)-x(2);
%v2= -2*x(1)-x(2)*(6*x(1)^2-2);


%v1=1./(cos(x(2))+2)*(-cos(x(2))*(x(1)-2*x(2))+4*(x(1)+sin(x(2))));
%v2=1./(cos(x(2))+2)*(x(1)-2*x(2)+2*(x(1)+sin(x(2))));


v1=x(2);
v2=-x(2)+x(1)-x(1)^3
%v2=-x(1)-0.5*x(2)*(1-(1+sin(x(1)))^2);

% v1=x(1)-x(1)^3+x(2);
% v2=2*x(1)-x(2);

%v1=x(2);
%v2=-x(1)-x(2)-x(2)*abs(x(2));


%(1-x(1)^2)*x(2)-x(1);%-x(2)-x(1)^3;

%v1=-x(1)*(1-x(1))-10*x(2);
%v2=-2*x(2)-10*x(1);%x(1)-x(2)-x(2)*abs(x(2));%-x(2)-x(1)^3;


% v1=-2*x(1)+x(1)^2+0.0*x(2)^2;
% v2=-6*x(2)+0.0*2*x(1)*x(2);
%x(1)-x(2)-x(2)*abs(x(2));%-x(2)-x(1)^3;

%%%%convex function+bump function
%a=100;

%v1=-(2*x(1)-a*exp(-(x(1)-5)^2-(x(2)-5)^2)*(x(1)-5));
%v2=-(2*x(2)-a*exp(-(x(1)-5)^2-(x(2)-5)^2)*(x(2)-5));


%Motzkin example x4y2 + x2y4 ? 3x2y2

%v1=-(4*x(1)^3*x(2)^2+2*x(1)*x(2)^4-6*x(1)*x(2));
%v2=-(2*x(1)^4*x(2)+x(1)^2*4*x(2)^3-6*x(1)^2*x(2));


% v1=-(2*x(1)-100*exp(-x(1)^2-x(2)^2)*x(1));
% 
% v2=-(2*x(2)-100*exp(-x(1)^2-x(2)^2)*x(2));
xdot=[v1;v2];