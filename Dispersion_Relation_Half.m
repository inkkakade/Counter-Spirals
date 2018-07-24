%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
close all
clc
clear

%% Constant Definition (keeping Basilisk values & units)
g1 = 10; g2 = 20; g3 = 30; g4 = 40; g5 = 50;                           %Acceleration due to gravity
rho = 1;                              %rho_l
Lo = 2;                               %Square Box
H = Lo/2;                             %Depth
lambda = 0.5;                         %Wavelength on Basilisk
%fo = sqrt(g/(2*pi*lambda));         
k = 2*pi/lambda;                      %Wave Number

max_sigma = 4.8;                      %Maximum value for Surface Tension (sigma) used in Basilisk)
delta_sigma = 0.2;                    %Step size of sigma
sigma  = [0:delta_sigma:max_sigma];   %from range of values considered in Basilisk
mydata = cell(1, length(sigma));

lcg1 = sqrt(sigma./(rho*g1));         %Numerical Capillary Length
lcg2 = sqrt(sigma./(rho*g2));          
lcg3 = sqrt(sigma./(rho*g3));         
lcg4 = sqrt(sigma./(rho*g4));          
lcg5 = sqrt(sigma./(rho*g5));        
                                      
L_by_lcg1 = lambda./lcg1;
L_by_lcg2 = lambda./lcg2;
L_by_lcg3 = lambda./lcg3;
L_by_lcg4 = lambda./lcg4;
L_by_lcg5 = lambda./lcg5;

klcg1 = 2*pi./L_by_lcg1;
klcg2 = 2*pi./L_by_lcg2;
klcg3 = 2*pi./L_by_lcg3;
klcg4 = 2*pi./L_by_lcg4;
klcg5 = 2*pi./L_by_lcg5;

%Automated file importing
for i = 0:delta_sigma:max_sigma
  myfilename = sprintf('data_17_07/Dispersion/Half/wave-%d-%g.txt', g1, i);
  mydatag1{fix(5*i+1.5)} = importdata(myfilename);
  myfilename = sprintf('data_17_07/Dispersion/Half/wave-%d-%g.txt', g2, i);
  mydatag2{fix(5*i+1.5)} = importdata(myfilename);
  myfilename = sprintf('data_17_07/Dispersion/Half/wave-%d-%g.txt', g3, i);
  mydatag3{fix(5*i+1.5)} = importdata(myfilename);
  myfilename = sprintf('data_17_07/Dispersion/Half/wave-%d-%g.txt', g4, i);
  mydatag4{fix(5*i+1.5)} = importdata(myfilename);
  myfilename = sprintf('data_17_07/Dispersion/Half/wave-%d-%g.txt', g5, i);
  mydatag5{fix(5*i+1.5)} = importdata(myfilename);
end

%% Finding the Period of the Wave Using Local Maximums G1
period = zeros(1,length(mydatag1));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydatag1{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freqg1 = 1./period; 
freq_anag1 = sqrt(klcg1.*(1+(klcg1.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Cog1 = sqrt(2*g1./lcg1);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_numg1 = freqg1.*sqrt(g1./lcg1)*lambda;
Cphase_anag1 = freq_anag1./klcg1;%.*sqrt(g./lc)*lambda;

%% Finding the Period of the Wave Using Local Maximums G2
period = zeros(1,length(mydatag2));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydatag2{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freqg2 = 1./period; 
freq_anag2 = sqrt(klcg2.*(1+(klcg2.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Cog2 = sqrt(2*g2./lcg2);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_numg2 = freqg2.*sqrt(g2./lcg2)*lambda;
Cphase_anag2 = freq_anag2./klcg2;%.*sqrt(g./lc)*lambda;

%% Finding the Period of the Wave Using Local Maximums G3
period = zeros(1,length(mydatag3));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydatag3{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freqg3 = 1./period; 
freq_anag3 = sqrt(klcg3.*(1+(klcg3.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Cog3 = sqrt(2*g3./lcg3);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_numg3 = freqg3.*sqrt(g3./lcg3)*lambda;
Cphase_anag3 = freq_anag3./klcg3;%.*sqrt(g./lc)*lambda;

%% Finding the Period of the Wave Using Local Maximums G4
period = zeros(1,length(mydatag4));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydatag4{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freqg4 = 1./period; 
freq_anag4 = sqrt(klcg4.*(1+(klcg4.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Cog4 = sqrt(2*g4./lcg4);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_numg4 = freqg4.*sqrt(g4./lcg4)*lambda;
Cphase_anag4 = freq_anag4./klcg4;%.*sqrt(g./lc)*lambda;

%% Finding the Period of the Wave Using Local Maximums G5
period = zeros(1,length(mydatag5));
for i = 1:length(sigma)
   total_timeperiods = 0;
   indices = [];
   values = [];
   r = mydatag5{i}; 
   [values,indices,f]= max_local(r(:,2),5,1);                              %Check Resolution Condition
   for j = 1:length(indices)-1
       total_timeperiods = total_timeperiods+(r(indices(j+1),1)-r(indices(j),1)); %Finding total time period of a single wave
   end 
   period(i) = period(i) + total_timeperiods/length(indices);
end

%Reference from Kundu book on Fluid Mechanics
%c = sqrt(((g*period/2*pi) + (2*pi*sigma)./(rho*period)).*tanh(2*pi*H./period));
freqg5 = 1./period; 
freq_anag5 = sqrt(klcg5.*(1+(klcg5.^2)))/(2*pi);  %Analytical relation 
% Dimensionalizing numerical variables to compare with actually dispersion
% relation
Cog5 = sqrt(2*g5./lcg5);
%Co = 50;
%C =  Co.*sqrt(1./klc + klc);
Cphase_numg5 = freqg5.*sqrt(g5./lcg5)*lambda;
Cphase_anag5 = freq_anag5./klcg5;%.*sqrt(g./lc)*lambda;

%% Plotting G1
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag1)
    data = mydatag1{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G1')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg1,freqg1,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg1.^2));
loglog(klcg1,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G1')

%Plotting Analytical Dispersion Relation
figure
plot(klcg1,Cphase_numg1,'ok')
hold on
plot(klcg1,Cphase_anag1, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on%% Plotting
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag1)
    data = mydatag1{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G1')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg1,freqg1,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg1.^2));
loglog(klcg1,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G1')

%Plotting Analytical Dispersion Relation
figure
plot(klcg1,Cphase_numg1,'ok')
hold on
plot(klcg1,Cphase_anag1, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1 G1')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag1{1};
cdata = mydatag1{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G1')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg1,Cphase_numg1,'ok') %./8.4515
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G1')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg1/2.4,Cphase_numg1/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G1')
xlabel('klc')
ylabel('Phase Velocity')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag1{1};
cdata = mydatag1{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G1')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg1,Cphase_numg1./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G1')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg1/2.4,Cphase_numg1/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G1')

%% Plotting G2
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag2)
    data = mydatag2{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G2')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg2,freqg2,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg2.^2));
loglog(klcg2,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G2')

%Plotting Analytical Dispersion Relation
figure
plot(klcg2,Cphase_numg2,'ok')
hold on
plot(klcg2,Cphase_anag2, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1 G2')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag2{1};
cdata = mydatag2{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G2')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg2,Cphase_numg2./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G2')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg2/2.4,Cphase_numg2/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G2')

%% Plotting G3
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag3)
    data = mydatag3{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G3')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg3,freqg3,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg3.^2));
loglog(klcg3,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G3')

%Plotting Analytical Dispersion Relation
figure
plot(klcg3,Cphase_numg3,'ok')
hold on
plot(klcg3,Cphase_anag3, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1 G3')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag3{1};
cdata = mydatag3{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G3')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg3,Cphase_numg3./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G3')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg3/2.4,Cphase_numg3/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G3')

%% Plotting G4
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag4)
    data = mydatag4{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G4')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg4,freqg4,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg4.^2));
loglog(klcg4,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G4')

%Plotting Analytical Dispersion Relation
figure
plot(klcg4,Cphase_numg4,'ok')
hold on
plot(klcg4,Cphase_anag4, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1 G4')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag4{1};
cdata = mydatag4{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G4')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg4,Cphase_numg4./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G4')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg4/2.4,Cphase_numg4/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G4')

%% Plotting G5
%Plotting all Basilisk Data
figure
hold on
for i = 1:length(mydatag5)
    data = mydatag5{i};
    plot(data(:,1),data(:,2))
end 
xlabel('Time')
ylabel('Amplitude')
title('Wave Forms for Increasing Surface Tensions G5')
hold off

%Plotting Numerical Dispersion Relation
figure
loglog(klcg5,freqg5,'ok')
hold on
freq_th = 1/pi*sqrt((1+klcg5.^2));
loglog(klcg5,freq_th,'r')
legend('Numerics','Analytical')
%loglog(klc,fo*(klc/ko).^(1))
grid on
xlabel('klc')
ylabel('Frequency')
title('Numerical Dispersion Relation in log-log scale G5')

%Plotting Analytical Dispersion Relation
figure
plot(klcg5,Cphase_numg5,'ok')
hold on
plot(klcg5,Cphase_anag5, '-r')
%plot(klc,C,'r--')
legend('Numerics','Analytical')
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Analytical Dispersion around klc = 1 G5')

%Comparison of Surface Tension Effects on Waves 
figure
hold on
gdata = mydatag5{1};
cdata = mydatag5{end};
plot(gdata(:,1),gdata(:,2),40+cdata(:,1),cdata(:,2))
legend('Gravity Waves','Capillary Waves')
xlabel('Time')
ylabel('Amplitude')
title('Effect of Surface Tension G5')
hold off

%Fully Normalised Plotting Analytical Dispersion Relation (with Co and lc)
figure
plot(klcg5,Cphase_numg5./8.4515,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion G5')
axis([0 2 0.75 1.5])

%Fully Normalised Plotting Analytical Dispersion Relation (Water)
figure
plot(klcg5/2.4,Cphase_numg5/23,'ok')
hold on
grid on
xlabel('klc')
ylabel('Phase Velocity')
title('Fully Normalised Analytical Dispersion for Water G5')