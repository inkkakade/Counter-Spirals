%%
% *Post Processing of 3D Simulation from Basilisk Data by Ishaan Kakadé*
clc
clear
close all 

%% Auto Import from Basilisk
total_t = 5.095;                     
delta_t = 0.005;
 
for i = 0:delta_t:total_t
  try
  myfilename = sprintf('data_23_07/R5S2/datafacqcc-%g.txt', i);
  mydata{fix(1+i/delta_t)} = importdata(myfilename);
  catch
     if isempty(mydata{fix(1+i/delta_t)})
      warning('Data loss on %g time step', i)
     end 
  end
end

%% Fast Plot
az = 90;
el = 90;
figure
for i = 1:length(mydata)
    try
    data = mydata{i};
    xdata = data(:,1);
    ydata = data(:,2);
    zdata = data(:,3);
    sf = fit([xdata, ydata],zdata,'linearinterp');
    plot(sf)
    shading interp;
    view([az,el])
    axis([-100 100 -100 100 -8  8])
    pause(0.001);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    catch 
        warning('Data missing at %d time step', i)
    end
end      

%% Plotting with GIF image export (Slow)

az = 90;
el = 90;
h = figure;
axis([-100  100 -100 100 -8 8])
 % this ensures that getframe() returns a consistent size
filename = 'R5S2.gif';
for i = 1:length(mydata)
    data = mydata{i};
    try
    xdata = data(:,1);
    ydata = data(:,2);
    zdata = data(:,3);
    sf = fit([xdata, ydata],zdata,'linearinterp');
    plot(sf)
    shading interp;
    view(az,el)
    axis([-100  100 -100 100 -8  8])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
    catch
        warning('boom')
    end
end

%% Temporal Analysis

plotdata = [];
figure
for i = 1:length(mydata)
    try
    data = mydata{i};
    dplot = data(1000,3);
    plotdata(i) = dplot;
    catch
        warning('Data missing at %d time step', i)
    end
end
plot(plotdata)
