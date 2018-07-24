%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
%% Basilisk Animated Post Processing
close all
clear
clc

%% Auto Import from Basilisk
g = 50; %m/s^2
Lo = 10;
rho = 1; %kg/m^3
timestep = 0.005;
maximumtimestep = 5.095;
sigma = 1;
a_eps = [];

for i = 0:timestep:maximumtimestep                                          %range = total time; step = time step
     myfilename = sprintf('data_16_07/Static/sigma0/wavepos-%g.txt', i*1000);
     mydata{fix(1+i/timestep)} = importdata(myfilename);
end

myfilename = sprintf('data_16_07/Static/sigma0/wave-50-1e-09.txt');
amplitude{1} = importdata(myfilename);

%Data Cleaning and Acquiring Maximum Depth
for i = 1:length(mydata)
    datax = mydata{i};
    try
        a_eps(i) = -min(datax(:,2));
    catch
    end
end

lc = sqrt(sigma/(rho*g)); %Numerical Capillary Length

%% Plotting
%Surface Waves Plot
figure
title('Basilisk Output for Interface Deformation')
for i = 1:length(mydata)
    data = mydata{i};
    try
        plot(data(:,1),data(:,2))
    catch
       %warning('Data is Lost in Time Step Number %d', i)
    end
    axis([0 20 -2 2])
    pause(0.01);     
end
hold on
%Pressure Plot (Comment out for Transient Plotting)
x = (0:0.05:20);
y = 0.3+exp(-((x-10).^2)/0.5^2);
plot(x,y)
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
%hold off

%% Verification
%Convergence Check
figure
amp = amplitude{1};
plot(amp(1:2000,1),amp(1:2000,2))
xlabel('Time')
ylabel('Maximum Amplitude')

%Static Pressure Equilibrium Check (Integral Method)
fun = @(x)(-exp(-((x-5).^2)/0.5^2));
pressure_int = integral(fun, 0,Lo) %#ok<NOPTS>
xval = data(:,1);
yval = data(:,2);
f = fit(xval,yval,'gauss1') %#ok<NOPTS>
figure
plot(f,x,y)
coeff = coeffvalues(f);
fit_fun = @(x)(coeff(1)*0.95*exp(-((x-coeff(2)).^2)/(0.95*coeff(3))^2));
interface_int = integral(fit_fun,0,Lo) %#ok<NOPTS>