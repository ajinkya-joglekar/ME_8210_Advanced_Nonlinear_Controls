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
numerator2 = [-1,1];
denominator2 = [1, -0.9, 1];

% Create the transfer function object
sys2 = tf(numerator2, denominator2);

% Plot the Nyquist plot
nyquist(sys2);
grid on;
title('Nyquist Plot for 11B');
