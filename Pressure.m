%% Basilisk Animated Post Processing
close all
clear
clc

%% Auto Import from Basilisk
timestep = 0.005; %0.05
maximumtimestep = 5.095;%9.95
sigma = 0.0001;

for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_02_07/Static/wavepos-%g.txt', i*1000);
     mydata{fix(1+i/timestep)} = importdata(myfilename);
end

myfilename = sprintf('data_02_07/Static/wave-128-%g.txt', sigma);
amplitude{1} = importdata(myfilename);
%% Plotting
%Surface Waves Plot
figure
for i = 1:length(mydata)
    data = mydata{i};
    try
    plot(data(:,1),data(:,2))
    catch
       %warning('Data is Lost in Time Step Number %d', i)
    end
    axis([0 10 -2 1])
    pause(0.01);     
end
hold on
%Pressure Plot (Comment out for Transient Plotting)
x = (0:0.05:10);
y = -exp(-((x-5).^2)/0.5^2);
plot(x,y)
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
%hold off
%% Verification
%Convergence Checkfit_fun = @(x)(coeff(1)*exp(-((x-coeff(2)).^2)/coeff(3)^2));
figure
amp = amplitude{1};
plot(amp(:,1),amp(:,2))
xlabel('Time')
ylabel('Maximum Amplitude')

%Static Pressure Equilibrium Check (Integral Method)
fun = @(x)(-exp(-((x-5).^2)/0.5^2));
pressure_int = integral(fun, 0,5);
xval = data(:,1);
yval = data(:,2);
f = fit(xval,yval,'gauss1') %#ok<NOPTS>
figure
plot(f,x,y)
coeff = coeffvalues(f);
fit_fun = @(x)(coeff(1)*exp(-((x-coeff(2)).^2)/coeff(3)^2));
interface_int = integral(fit_fun,0,5);


