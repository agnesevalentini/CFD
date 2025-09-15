%comparison between u profile in analytical solution and simulation



data = load("u_poiseuille.dat");
Ux = data(:,1);
y = data(:,10);



%y = 0:0.01:1;
%parameters of the function

nu = 0.05;
b = 0.5;
dpdx = -4.98;

plot(y, Ux , 'b-', 'LineWidth', 40, 'DisplayName', 'Simulation');
hold on;
plot(y, u(y, nu, b, dpdx), 'y-', 'LineWidth', 2, 'DisplayName', 'Analytical solution');
hold off;
xlim([0 1]);
%ylim([0 1]);
xlabel("y");
ylabel("u");
legend('Location', 'southeast');
%title("Velocity profile");

