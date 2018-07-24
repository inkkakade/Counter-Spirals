%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
close all
clc
clear

%% Auto Import from Basilisk Data
timestep = 0.005; %0.05
maximumtimestep = 5.000;%9.95

for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_10_07/new/Speed3/wavepos-%g.txt', i*1000);
     mydata{fix(1+i/timestep)} = importdata(myfilename);
end

%% Plotting with GIF image export
h = figure;
axis([0 20 -2 1]) % this ensures that getframe() returns a consistent size
filename = 'S3.gif';
for i = 1:length(mydata)-600
    data = mydata{i};
    try
    x = data(:,1);
    y = data(:,2);
    plot(x,y,'ob')
    axis([0 20 -2 1])
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
    end
end

%% Quick Plot
figure
for i = 1:length(mydata)
    data = mydata{i};
    try
    plot(data(:,1),data(:,2), 'ob')
    catch
    end
    axis([0 20 -3 1])
    pause(0.001); 
    xlabel('Length')
    ylabel('Displacement')
    title('Basilisk Output')
end
