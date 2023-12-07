fig1 = openfig('C:\Users\Manny\MATLAB Drive\2.74FinalProject\best_dynamics\theta_1.fig');
axes1 = get(fig1, 'Children');
data1 = get(axes1, 'Children');
x1 = get(data1, 'XData');
y1 = get(data1, 'YData');

fig2 = openfig('C:\Users\Manny\MATLAB Drive\2.74FinalProject\best_dynamics\theta_2.fig');
axes2 = get(fig2, 'Children');
data2 = get(axes2, 'Children');
x2 = get(data2, 'XData');
y2 = get(data2, 'YData');

newFig = figure(13);
plot(x1, y1, 'b-', 'DisplayName', 'Theta 1'); % blue line for fig 1 data
hold on;
plot(x2, y2, 'r-', 'DisplayName', 'Theta 2'); % red dashed line for fig 2 data
hold off;

legend('show');

xlabel('Time (s)');
ylabel('Angle (rad)');
title('Angle over time');