%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
close all
clc
clear

%% Constant Definition (keeping Basilisk values & units)
g =10;                             %Acceleration due to gravity
rho = 1;                            %rho_l
Lo = 2;                             %Square Box
H = Lo/2;                           %Depth
lambda = 1;                      %Wavelength on Basilisk
%fo = sqrt(g/(2*pi*lambda));         
k = 2*pi/lambda;                    %Wave Number

max_sigma = 4.8;                      %Maximum value for Surface Tension (sigma) used in Basilisk)
delta_sigma = 0.2;                 %Step size of sigma
sigma  = [0:delta_sigma:max_sigma]; %from range of values considered in Basilisk
mydata = cell(1, length(sigma));
lc = sqrt(sigma./(rho*g));          %Numerical Capillary Length
L_by_lc = lambda./lc;
klc = 2*pi./L_by_lc;

%Automated file importing
for i = 0:delta_sigma:max_sigma
  myfilename = sprintf('data_17_07/Dispersion/One/wave-%d-%g.txt', g, i);
  mydata{fix(5*i+1)} = importdata(myfilename);
end

%% Finding the Period of the Wave Using Local Maximums 
period = zeros(1,length(mydata));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydata{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freq = 1./period; 
freq_ana = sqrt(klc.*(1+(klc.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Co = sqrt(2*g./lc);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_num = freq.*sqrt(g./lc)*lambda;
Cphase_ana = freq_ana./klc*35.56;%.*sqrt(g./lc)*lambda;

%% Plotting
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydata)
    data = mydata{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klc,freq,'ok')
hold on
freq_th = 1/pi*sqrt((1+klc.^2));
loglog(klc,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale')

%Plotting Analytical Dispersion Relation
figure
plot(klc,Cphase_num,'ok')
hold on
plot(klc,Cphase_ana, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydata{1};
cdata = mydata{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klc,Cphase_num./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klc/2.4,Cphase_num/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water')