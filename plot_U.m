% function to plot the U to check the steady state in the Poiseuille case (must be in the same folder as the data to load)
% first idea didn't work
%data = load("U");
%t = data(:,1);
%U = data(:,2);
%plot(t,U)

%xlim([0 25]);
%xlabel('t');
%ylabel('U');


fid = fopen('U','r');

% Skip lines starting with '#'
raw = textscan(fid, '%f (%f %f %f)', 'CommentStyle', '#');

fclose(fid);

% Extract columns
time = raw{1};
Ux   = raw{2};
Uy   = raw{3};
Uz   = raw{4};

% Plot
figure;
plot(time, Ux, 'b-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('U');
%title('Probe velocity component U_x vs Time');
grid on;

