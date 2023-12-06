% Path integral Eigenfunctions
% Anaytical Example
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% set up the system for positive lambda2
% system dynamics
lambda1=-0.5+1.32i; lambda2=-0.5-1.32i; %Linearization eigen values
f = @(t, x) [lambda1*(x(2,:)); lambda2*(x(1,:) -x(1,:).^3 - x(2,:))];

% convergence check
if(lambda2>0)
    if(2*lambda1-lambda2 < 0)
        disp('Convergece Condition is satified')
    else
        disp('Convergece Condition failed!!')
    end
end

if(lambda2<0)
    if(2*lambda1+lambda2<0)
        disp('Convergece Condition is satified')
    else
        disp('Convergece Condition failed!!')
    end
end

% setup for phase portrait
alpha=1;
Dom = [-2 2]; ds = 0.05;
[X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
u = alpha.*(lambda1.*Y);
v = alpha.*(lambda2.*(X-X.^3-Y));

%% Path integral setup
% determine the eigenvalues and left-eigenvectors
x = sym('x',[2;1], 'real');
A = double(subs(jacobian(f(0,x),x),x,[1;0])); %Point corresponding to the mentioned eigenvalues
fn = @(x) f(0,x)-A*x;
[~,D,W] = eig(A); 
l1 = D(1,1); l2 = D(2,2);
w1 = W(:,1); w2 = W(:,2);

% define the function for path integral
g1 = @(x) w1'*fn(x);
g2 = @(x) w2'*fn(x);

%% Eigenfunction computation using the path Integral
w_bar = waitbar(0,'1','Name','Calcualting path integral...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%define grid where eigenfunction is to be calculated
grid = Dom(1):ds:Dom(2); 
[q1,q2] = meshgrid(grid);
QQ = [q1(:)';q2(:)'];
U = f(0,QQ);

% setup the initial conditions and eigenfunction variables
x_0 = [q1(:),q2(:)]; 
phi1_est=[];phi2_est=[];

%% Find the Lyaponav function
Q = -1*eye(2); % Identity matrix as an example

% Solve for P using Lyapunov equation: A^T * P + P * A = -Q
LAMBDA = eig(A,"matrix");
P = lyap(LAMBDA, -Q);
V = [];
dV = [0];

% 2D Distance Function
distance_2d = @(x1, y1, x2, y2) sqrt((x2 - x1)^2 + (y2 - y1)^2);

%%

% ode options
options = odeset('RelTol',1e-9,'AbsTol',1e-10);

% path integral calculation
for i = 1:length(x_0)
    % update the progress bar
    waitbar(i/length(x_0),w_bar,sprintf(string(i)+'/'+string(length(x_0))))
    % simulate trajectory from initial condition 'x_0(i,:)'
    [t,x] = ode45(@(t,x)f(t,x),[0 3],x_0(i,:), options);
    % numerically approximate eigenfunction values
    phi1_est = [phi1_est, w1'*x_0(i,:)' + trapz(t,exp(-l1*t).*g1(x')')];
    phi2_est = [phi2_est, w2'*x_0(i,:)' + trapz(t,exp(-l2*t).*g2(x')')];
    V_t = [phi1_est(i),phi2_est(i)]*P*[phi1_est(i);phi2_est(i)];
    V = [V,V_t];
    if i > 1 
        dist = distance_2d(phi1_est(i), phi2_est(i), phi1_est(i-1), phi2_est(i-1));
        gradient = (V(i) - V(i-1))/dist;
        dV = [dV,gradient];
    end
end
% close the progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

% reshape the eigenfunction estimates into x1-x2 grid
phi1_est = reshape((phi1_est),size(q2));
phi1_est = real(phi1_est);
phi2_est = reshape((phi2_est),size(q2));
phi2_est = real(phi2_est);
V = reshape((V),size(q2));
V = real(V);
dV = reshape((dV),size(q2));
dV = real(dV);

%% plots the eigenfunctions
figure(1)
subplot(1,2,1)
% visualize the approximation of the 1st eigenfunction
p3 = pcolor(q1,q2,phi1_est); hold on;
set(p3,'Edgecolor','none')
colormap jet
% plot the phase portrait on top of the eigenfunction
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-2,2]); ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('$\phi_1$: $\lambda_1=$'+string(lambda1), 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar

subplot(1,2,2)
% visualize the approximation of the 2nd eigenfunction
p4 = pcolor(q1,q2,phi2_est); hold on;
set(p4,'Edgecolor','none')
colormap jet
% plot the phase portrait on top of the eigenfunction
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-2,2]); ylim([-2,2])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('$\phi_2$: $\lambda_2=$'+string(lambda2), 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar
%% Calculate gradient and ROI
% V_dot = gradient(V);
ROA = zeros(size(dV));
size_V_dot = size(dV);
for i = 1:size_V_dot(1)
    for j = 1:size_V_dot(1)
        currentElement = dV(i, j);
        if currentElement <= 0 
            ROA(i,j) = 1;
        else
            ROA(i,j) = 0;
        end
    end
end
%%

% Create a custom colormap with only green and red
customColormap = [1, 0, 0;  % Red
                  0, 1, 0]; % Green

figure(2)
subplot(1,2,1)
p1 = pcolor(q1,q2,V);
hold on
set(p1,'Edgecolor','none')
xlim([-2,2]); ylim([-2,2])
ax1 = gca;
colormap(ax1,"jet"); 
axis square
set(ax1,'FontSize',15);
l = streamslice(X,Y,u,v); 
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-2,2]); ylim([-2,2])
hold off;
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title(['Lyaponav function value over new basis'],['$\phi_1$: $\lambda_1=$'+string(lambda1),'$\phi_2$: $\lambda_2=$'+string(lambda2)], 'Interpreter','latex')
box on
ax1.LineWidth=2;
colorbar
hold off

subplot(1,2,2)

p2 = pcolor(q1,q2,ROA); hold on;
shading flat;  % Ensure the color is flat (no interpolation)
set(p2, 'Edgecolor', 'none');
xlim([-2,2]); ylim([-2,2])
ax2 = gca;
colormap(ax2,customColormap);
axis square
set(ax2,'FontSize',15);
l = streamslice(X,Y,u,v); 
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-2,2]); ylim([-2,2])
hold off;
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title(['Region of attraction over new basis'],['$\phi_1$: $\lambda_1=$'+string(lambda1),'$\phi_2$: $\lambda_2=$'+string(lambda2)], 'Interpreter','latex')
box on
ax2.LineWidth=2;
colorbar

%%
threshold = 1e-2;
stable_region = phi2_est;
stable_region(-threshold < stable_region < threshold) = 0;
stable_region(stable_region ~= 0) = 1;

figure(3)
p1 = pcolor(q1,q2,stable_region);
hold on
set(p1,'Edgecolor','none')
xlim([-2,2]); ylim([-2,2])
ax1 = gca;
colormap(ax1,"jet"); 
axis square
set(ax1,'FontSize',15);
l = streamslice(X,Y,u,v); 
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-2,2]); ylim([-2,2])
hold off;
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title(['Stable manifold of $\phi_2$ coordinates'],['corresponding to $\lambda_2=$'+string(lambda2)],'Interpreter','latex')
box on
ax1.LineWidth=2;
colorbar
hold off
