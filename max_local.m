function [values,indices,f]=max_local(r,resolution,pas)

counter=0;

% resolution = the mean resolution the curve on n points
% pas = the step is used to perform the sliding average on n points, 
%       to compare with the current value

%performs a rolling average on the curve (smoothing to avoid unwanted local maxima)
for i=1+resolution:length(r)-resolution
    f(i)=mean(r(i-resolution:i+resolution));
end

%Search for local maxima, sliding average between index = i-pas and index = i + pas
for i=pas+1:(length(r)-pas-2*resolution)
    [val,indice]=max(f(i-pas:i+pas));
    if ((pas+1) == indice)
        counter=counter+1;
        values(counter)=r(i);
        indices(counter)=i;
    end
end

