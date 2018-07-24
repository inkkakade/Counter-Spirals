%%
% *Post Processing of Dispersion Relation from Basilisk Data by Kakadé*
close all
clear
clc

%% Auto Import from Basilisk
g = 50; %m/s^2
rho = 1; %kg/m^3
Lo = 20;
a_charac = 0.5;
sigmastep = 0.1;
max_sigma = 5;

sigma = [0:0.1:5];
%sigma_low = [0:0.01:0.1];
sigma_high = [6:1:7];

lc = sqrt(sigma/(rho*g)); 
lc_comp = sqrt(0.01/(rho*g)); 

%lc_low = sqrt(sigma_low/(rho*g));
lc_high = sqrt(sigma_high/(rho*g));

for i = 0:sigmastep:max_sigma                                          %range = total time; step = time step
     myfilename = sprintf('data_04_07/Bond/wave-256-%g.txt', i);
     mydata{fix(10*i+1)} = importdata(myfilename);
end

for i = 6:1:7                                        %range = total time; step = time step
     myfilename = sprintf('data_11_07/Bond/wavefinalpos-%g.txt', i);
     mydatahigh{fix(i-5)} = importdata(myfilename);
end

%% Bond Number Plot
a = [];
for i = 1:length(mydata)
    data = mydata{i};
    maxi = max(data(:,2));
    a(i) = maxi;
end

for i = 1:length(mydatahigh)
    data = mydatahigh{i};
    maxi = min(data(:,2));
    a = cat(2,a,maxi);
end

Bo = (a_charac./lc).^2;
Bo(1) = (a_charac./lc_comp).^2;                                              %Inserting high Bond number for exteremely low surface tensions of the order of 10e-4
Bo = cat(2,Bo,(a_charac./lc_high).^2);

figure
plot(Bo(2:end),a(2:end),'ok')
grid on
title('Depth of Penetration versus Bond Number log-log Plot')
axis([0 140 -0.2 0.2])
ylabel('Depth of Penetration')
xlabel('Bond Number')

sigmastepbond = 0.5;
max_sigmabond = 5;
sigma_bond = [0:0.5:5];
lc_bond = sqrt(sigma_bond./rho*g);

for i = 0:sigmastepbond:max_sigmabond
    for j = 1:1:max_sigmabond/sigmastepbond
     mybondfilename = sprintf('data_05_07/Bond/wavefinalpos-%g.txt', i);
     mybonddata{j} = importdata(mybondfilename);
    end
end

int_vector = [];
for i = 1:length(mybonddata)
    bonddata = mybonddata{i};
    xvals = bonddata(:,1);
    yvals = bonddata(:,2);
    fitting = fittype('a*exp(-(b*(x - c)^2))+d');
    f = fit(xvals, yvals,fitting,'StartPoint',[1,0.7,10,0.05]);
    coeff = coeffvalues(f);
    fun = @(x)coeff(1).*exp(-coeff(2)*(x - coeff(3)).^2) + coeff(4);
    int_vol2D = integral(fun,0,Lo);
    int_vector(i) = int_vol2D;
end


bondvol = (a_charac./lc_bond).^2;    
%bondvol(1) = (a_charac/lc_comp).^2;
figure
plot(bondvol(2:end),int_vector)
grid on
title('Volume of Fluid Displaced versus Bond Number Plot For Mass Conservation')
ylabel('Volume of Water Integrated Across Interface')
xlabel('Bond Number')