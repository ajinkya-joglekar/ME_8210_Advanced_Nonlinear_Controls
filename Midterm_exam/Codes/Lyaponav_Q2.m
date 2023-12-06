% Define the nonlinear system
dx1 = @(x1, x2) x2;
dx2 = @(x1, x2) -x1 -x1.^3 - x2;

% Define the Lyapunov function
V = @(x1, x2) 7/4 * x1.^2 + 2/4 *x1*x2 + 3/4*x2.^2;

% Define the gradient of the Lyapunov function
dV = @(x1, x2) ((7/2 * x1 + 1/2 * x2)*dx1(x1, x2) + (3/2 * x2 + 1/2 * x1)*dx2(x1, x2));

% Set up a grid of initial conditions
[x1, x2] = meshgrid(-3:0.1:3, -3:0.1:3);

[xtraj1,xtraj2] = meshgrid(-3:0.25:3,-3:0.25:3);

dV_grid = []; % Gradient of Lyaponav function across the grid


for i = 1:numel(x1)
    x1_0 = x1(i);
    x2_0 = x2(i);

    % Simulate the system using ode45 or another ODE solver
    [~, X] = ode45(@(t, x) [dx1(x(1), x(2)); dx2(x(1), x(2))], [0, 0.25], [x1_0; x2_0]);

    % Check if the Lyapunov function decreases or stays the same over time
    dV_val = dV(x1_0, x2_0);
    dV_grid = [dV_grid,dV_val];
end

dV_grid = reshape((dV_grid),size(x1));


%% Phase Portrait
% setup for phase portrait
alpha=1;
lambda1=-0.5; lambda2=-0.5; %Linearization eigen values
Dom = [-3 3]; ds = 0.05;
[X,Y] = meshgrid(Dom(1):0.5:Dom(2),Dom(1):0.5:Dom(2));
u = alpha.*(lambda1.*Y);
v = alpha.*(lambda2.*(X-X.^3-Y));

%% Find a new ROA
ROA = zeros(size(dV_grid));
size_dV_grid = size(dV_grid);
for i = 1:size_dV_grid(1)
    for j = 1:size_dV_grid(1)
        currentElement = dV_grid(i, j);
        if currentElement <= 0 
            ROA(i,j) = 1;
        else
            ROA(i,j) = 0;
        end
    end
end

%% Plot results

% Create a custom colormap with only green and red
customColormap = [1, 0, 0;  % Red
                  0, 1, 0]; % Gren


% plot the phase portrait on top of the lyaponav function value
figure(1)
subplot(1,2,1)
% visualize the approximation of the 1st eigenfunction
p3 = pcolor(x1,x2,dV_grid); hold on;
set(p3,'Edgecolor','none')
colormap jet
% plot the phase portrait on top of the eigenfunction
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1); set(l,'Color','k');
xlim([-3,3]); ylim([-3,3])
axes = gca;
axis square
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('Lyaponav function value across grid', 'Interpreter','latex')
box on
axes.LineWidth=2;
colorbar


subplot(1,2,2)

p2 = pcolor(x1,x2,ROA); hold on;
shading flat;  % Ensure the color is flat (no interpolation)
set(p2, 'Edgecolor', 'none');
xlim([-3,3]); ylim([-3,3])
ax2 = gca;
colormap(ax2,customColormap);
axis square
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1); set(l,'Color','k');
set(ax2,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
title('Region of attraction using Lyaponav function', 'Interpreter','latex')
box on
ax2.LineWidth=2;
colorbar

