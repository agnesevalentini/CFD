% script to plot the velocity and make comparison with Ghia et al (must be in the same folder as the data to load)


%data = load("horizontal_1000.dat");
data = load("vertical_1000.dat");
u = data(:,1);
y = 0:0.01:1;
%v = data(:,2);
%x = data(:,7);
%plot(x,u)
%hold on;

% data for Re = 100 for both horizontal and vertical velocity
%% plot for v
%x_ref = [1.000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000];
%v_ref = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000];

% plot for u
%y_ref = [1.00000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000];   % y positions
%u_ref = [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, -0.00000]; % u(y) values

% -------------------------------------

% data for Re = 1000 for both horizontal and vertical velocity

y_ref = [1.00000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000];   % y positions
u_ref = [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000];

%x_ref = [1.000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000];
%v_ref = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000];

% Plot simulation vs paper
plot(y, u, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation');
%plot(x, v, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation');
hold on;
plot(y_ref, u_ref, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Data from Ghia et al.');
%plot(x_ref, v_ref, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Data from Ghia et al.');
hold off;

%xlabel('x');
%ylabel('v');
xlabel('y');
ylabel('u');
%title('Lid-driven cavity: v at y = 0.5');
legend('Location', 'southeast');



