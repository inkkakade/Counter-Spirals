%%
% *Post Processing of Dispersion Relation from Basilisk Data by Ishaan Kakadé*
%% Basilisk Animated Post Processing
close all
clear
clc

%% Auto Import from Basilisk
timestep = 0.005; %0.05
maximumtimestep = 5.095;%9.95

for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_04_07/Transient/S1/wavepos-%g.txt', i*1000);
     mydata{fix(1+i/timestep)} = importdata(myfilename);
end

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
    axis([0 20 -2 1])
    pause(0.001);     
end
hold on

%% Data Import: Varying Speed of Wave
%S = 1
for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_04_07/Transient/S1/wavepos-%g.txt', i*1000);
     mydataS1{fix(1+i/timestep)} = importdata(myfilename);
end
%%
%S = 3
for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_04_07/Transient/S3/wavepos-%g.txt', i*1000);
     mydataS3{fix(1+i/timestep)} = importdata(myfilename);
end
%%
%S = 6
for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_04_07/Transient/S6/wavepos-%g.txt', i*1000);
     mydataS6{fix(1+i/timestep)} = importdata(myfilename);
end
%%
%S = 11
for i = 0:timestep:maximumtimestep %range = total time; step = time step
     myfilename = sprintf('data_04_07/Transient/S11/wavepos-%g.txt', i*1000);
     mydataS11{fix(1+i/timestep)} = importdata(myfilename);
end        
%% Surface Waves Plot
% figure
% for i = 1:length(mydataS1)
%     data = mydataS1{i};
%     try
%     plot(data(:,1),data(:,2), 'ob')
%     catch
%        %warning('Data is Lost in Time Step Number %d', i)
%     end
%     axis([0 20 -3 1])
%     pause(0.001); 
%     xlabel('Length')
%     ylabel('Displacement')
%     title('Basilisk Output')
% end
% 

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for i = 1:length(mydataS1)
    data = mydataS1{i};
    try
    x = data(:,1);
    y = data(:,2);
    plot(x,y) 
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


% 
% position = 10;
% c = 0;
% for i=1:length(mydata)
%     plotdata = mydataS1{i};   
%     size(plotdata);
%     if ~isempty(plotdata)
%         c=c+1;
%         [val,indice] = min(abs(plotdata(:,1)-position));
%         Y(c) = plotdata(indice,2);
%     end
% end
% plot(Y)