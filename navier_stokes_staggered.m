addpath('./functions/')
%% Parameters
Re = 1000; % Reynolds number  (also tried with Re = 1000)
nt = 25000; % max time steps (20000 for Re = 1000)
Lx = 1; Ly = 1; % domain size
Nx = 128; Ny = 128; % Number of grids
dt = 0.001; % time step;

%% Grid
% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;

%% Data arrays
u = zeros(Nx+1,Ny+2); % velocity in x direction (u)
v = zeros(Nx+2,Ny+1); % velocity in y direction (v)
p = zeros(Nx,Ny); % pressure (lagurange multiplier)
uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center

%% Time simulation
for ii = 1:nt
    % Setting BCs
    bctop = 1; % Top lid u
    u(:,1) = -u(:,2); v(:,1) = 0;             %bottom
    u(:,end) = 2*bctop-u(:,end-1);  v(:,end) = 0;  %top
    u(1,:) = 0;    v(1,:) = -v(2,:);             %left
    u(end,:) = 0;  v(end,:) = -v(end-1,:);    %right

    % Diffusive term (operator L)
    Lux = (u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1))/dx^2; % nx-1 * ny
    Luy = (u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end))/dy^2; % nx-1 * ny
    Lvx = (v(1:end-2,2:end-1)-2*v(2:end-1,2:end-1)+v(3:end,2:end-1))/dx^2; % nx * ny-1
    Lvy = (v(2:end-1,1:end-2)-2*v(2:end-1,2:end-1)+v(2:end-1,3:end))/dy^2; % nx * ny-1

    % Advective term (operator N)
    % 1. interpolate velocity at cell center/cell corner
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    uco = (u(:,1:end-1)+u(:,2:end))/2;
    vco = (v(1:end-1,:)+v(2:end,:))/2;
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
    % 2. multiply
    uuce = uce.*uce;
    uvco = uco.*vco;
    vvce = vce.*vce;
    % 3-1. get derivative for u
    Nu = (uuce(2:end,:) - uuce(1:end-1,:))/dx;
    Nu = Nu + (uvco(2:end-1,2:end) - uvco(2:end-1,1:end-1))/dy;
    % 3-2. get derivative for v
    Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
    Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;

    % Get intermidiate velocity (called u* before)
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + dt*(-Nu + (Lux+Luy)/Re);
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) + dt*(-Nv + (Lvx+Lvy)/Re);

    % update velocity with Poisson Equation for pressure
    % RHS of pressure Poisson eq.
    b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
        + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);
    % Solve for p
    p = solvePoissonEquation_2dDCT(b,Nx,Ny,dx,dy); % (using cosine transform, faster)
    % p = solvePoissonEquation_direct(b,Nx,Ny,dx,dy); % Direct method
    % The new divergent free velocity field
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (p(2:end,:)-p(1:end-1,:))/dx;
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (p(:,2:end)-p(:,1:end-1))/dy;

    % get velocity at the cell center (for visualization)
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
end

%% Setting Visualization for the contourplot
[Xce,Yce] = meshgrid(xce,yce); % cell center
[~,h_abs] = contourf(Xce',Yce',uce, -0.4:0.1:1.1);
xlim([0 Lx]); ylim([0 Ly]);
colorbar ();
% colorbar("ticks", levels);


%% COmparison plot between simulation and paper for u component of velocity
% Extract the line of uce corresponding to x=0.5
u_centerline = uce(64,:);
% y coordinate of the centers of the grid squares
y_plot = (dy/2):dy:(1-(dy/2));

% reference data
% Reference data from paper (y,u)
y_ref = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000];   % y positions
%u_ref Reynolds 100
%u_ref = [1.0000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, -0.00000]; % u(y) values

%u_ref Reynolds 1000
u_ref = [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000];


% Plot simulation vs paper
plot(y_plot, u_centerline, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation');
hold on;
plot(y_ref, u_ref, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', 'Data from Ghia et al.');
hold off;

xlabel('y position');
ylabel('u(y) velocity');
title('Lid-driven cavity: u(y) at x = 0.5 with Re = 1000');
legend('Location', 'southeast');



%% COmparison plot between simulation and paper for v component of velocity
% Extract the line of uce corresponding to y = 0.5
v_centerline = vce(:,64);
% x coordinate of the centers of the grid squares
x_plot = (dx/2):dx:(1-(dx/2));


%% plot for v
x_ref = [1.000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000];

% v_ref Reynolds 100
%v_ref = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000];

%v_ref Reynolds 1000
v_ref = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000];


% Plot simulation vs paper
plot(x_plot, v_centerline, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulation');
hold on;
plot(x_ref, v_ref, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Data from Ghia et al.');
hold off;

%xlabel('x position');
%ylabel('v velocity');
%title('Lid-driven cavity: v at y = 0.5');
%legend('Location', 'southeast');



