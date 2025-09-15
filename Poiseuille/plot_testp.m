% function to plot the p to check the steady state (must be in the same folder as the data to load)
data = load("p");
t = data(:,1);
p = data(:,2);
plot(t,p)

xlim([0 25]);
xlabel('t');
ylabel('p');
