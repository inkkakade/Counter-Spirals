%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
clc
clear
close all

timestep = 0.005;
maximumtimestep = 5.095;

for i = 0:timestep:maximumtimestep                                          %range = total time; step = time step
     myfilename = sprintf('data_16_07/Static/sigma0/wavepos-%g.txt', i*1000);
     mydata1{fix(1+i/timestep)} = importdata(myfilename);
end

figure
for i = 1:length(mydata1)
    data = mydata1{i};
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
title('Basilisk Output for Interface Deformation at sigma = 0')
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
hold off

for i = 0:timestep:maximumtimestep                                          %range = total time; step = time step
     myfilename = sprintf('data_16_07/Static/sigmaalc/wavepos-%g.txt', i*1000);
     mydata2{fix(1+i/timestep)} = importdata(myfilename);
end

figure
for i = 1:length(mydata2)
    data = mydata2{i};
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
title('Basilisk Output for Interface Deformation at sigma such that a = lc')
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
hold off

for i = 0:timestep:maximumtimestep                                          %range = total time; step = time step
     myfilename = sprintf('data_16_07/Static/sigmaalc10/wavepos-%g.txt', i*1000);
     mydata3{fix(1+i/timestep)} = importdata(myfilename);
end

figure
for i = 1:length(mydata3)
    data = mydata3{i};
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
title('Basilisk Output for Interface Deformation at sigma such that a/lc = 10')
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
hold off

for i = 0:timestep:maximumtimestep                                          %range = total time; step = time step
     myfilename = sprintf('data_16_07/Static/sigmaalc75/wavepos-%g.txt', i*1000);
     mydata3{fix(1+i/timestep)} = importdata(myfilename);
end

figure
for i = 1:length(mydata3)
    data = mydata3{i};
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
title('Basilisk Output for Interface Deformation at very high sigma ')
legend('Basilisk Solution','Applied Pressure Shape')
xlabel('Length of Box')
ylabel('Height')
hold off
%% Amplitude Convergence
myfilename = sprintf('data_16_07/Static/sigma0/wave-50-1e-09.txt');
amplitude{1} = importdata(myfilename);

myfilename = sprintf('data_16_07/Static/sigmaalc/wave-50-1.25.txt');
amplitude{2} = importdata(myfilename);

myfilename = sprintf('data_16_07/Static/sigmaalc10/wave-50-12.5.txt');
amplitude{3} = importdata(myfilename);

myfilename = sprintf('data_16_07/Static/sigmaalc75/wave-50-75.txt');
amplitude{4} = importdata(myfilename);

figure
amp = amplitude{1};
plot(amp(:,1),amp(:,2))
xlabel('Time')
ylabel('Maximum Amplitude')
axis([0 20 0 15])
title('Basilisk Output for Interface Deformation at sigma = 0')

figure
amp = amplitude{2};
plot(amp(:,1),amp(:,2))
xlabel('Time')
ylabel('Maximum Amplitude')
title('Basilisk Output for Interface Deformation at sigma such that a/lc = 1')

figure
amp = amplitude{3};
plot(amp(:,1),amp(:,2))
xlabel('Time')
ylabel('Maximum Amplitude')
title('Basilisk Output for Interface Deformation at sigma such that a/lc = 0.1')

figure
amp = amplitude{4};
plot(amp(:,1),amp(:,2))
xlabel('Time')
ylabel('Maximum Amplitude')
title('Basilisk Output for Interface Deformation at very high sigma')

