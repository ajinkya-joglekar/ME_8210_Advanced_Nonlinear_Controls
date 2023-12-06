% Define the transfer function
%% Question 11A
numerator1 = [1];
denominator1 = [1, 1, 1];

% Create the transfer function object
sys1 = tf(numerator1, denominator1);


% Plot the Nyquist plot
nyquist(sys1);
grid on;
title('Nyquist Plot for 11A');
% Add the transfer function to the title
% title(['Nyquist Plot of Transfer Function: ', num2str(numerator1), '/', num2str(denominator2)]);

%% Question 11B
numerator2 = [1,-1];
denominator2 = [1, 0.1, 0.9];

% Create the transfer function object
sys2 = tf(numerator2, denominator2);

% B = sys2/(1+sys2);
% minreal(B);

% Plot the Nyquist plot
nyquist(sys2);
grid on;
title('Nyquist Plot for 11B');

%%
% Define the transfer function
numerator = [1, 0, -1];
denominator = conv([1, 1], conv([1, 0, 1], [1, 0, -1]));

% Create the transfer function object
sys = tf(numerator, denominator);

% Plot the Nyquist plot
nyquist(sys);
grid on;

% Add the critical point (-1, j0) to the plot
hold on;
plot(-1, 0, 'ro'); % Critical point
hold off;

% Add labels and title
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Nyquist Plot');

% Determine stability using the circle criterion
critical_point = -1 + 0i;
is_stable = is_stable_circle_criterion(sys, critical_point);

% Display stability result
if is_stable
    disp('The system is stable according to the circle criterion.');
else
    disp('The system is not stable according to the circle criterion.');
end


%%
num = [1,0,-1];
den = [1,1,1,1];
w = 3; % for example 3 rad/s
val = polyval(num,j*w)/polyval(den,j*w);


%%
% Function to check stability using the circle criterion
function stability = is_stable_circle_criterion(sys, critical_point)
    % Find the frequency response of the system
    [magnitude, phase] = bode(sys);

    % Convert magnitude and phase to complex numbers
    Gjw = magnitude .* exp(1i * deg2rad(phase));

    % Check if the Nyquist plot encircles the critical point
    stability = any(encircle(Gjw, critical_point));
end
