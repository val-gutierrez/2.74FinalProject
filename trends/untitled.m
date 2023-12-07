

fig1 = openfig('C:\Users\Manny\MATLAB Drive\2.74FinalProject\trends\momentumOverMassWK01.fig');
axes1 = get(fig1, 'Children');
data1 = get(axes1, 'Children');
x1 = get(data1, 'XData');
y1 = get(data1, 'YData');

fig2 = openfig('C:\Users\Manny\MATLAB Drive\2.74FinalProject\trends\momentumOverMassWK05.fig');
axes2 = get(fig2, 'Children');
data2 = get(axes2, 'Children');
x2 = get(data2, 'XData');
y2 = get(data2, 'YData');

fig3 = openfig('C:\Users\Manny\MATLAB Drive\2.74FinalProject\trends\momentumOverMassWK1.fig');
axes3 = get(fig3, 'Children');
data3 = get(axes3, 'Children');
x3 = get(data3, 'XData');
y3 = get(data3, 'YData');

newFig = figure;
plot(x1, y1, 'b-', 'DisplayName', 'k = .01'); % blue line for fig 1 data
hold on;
plot(x2, y2, 'r-', 'DisplayName', 'k = .05'); % red dashed line for fig 2 data
plot(x3, y3, 'g-', 'DisplayName', 'k = .1'); % red dashed line for fig 2 data
hold off;

legend('show');

xlabel('Mass (kg)');
ylabel('Momentum (kg·m/s)');
title('Momentum of End-Effector vs Mass');
