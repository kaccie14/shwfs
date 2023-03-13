function [x, y, dx, dy] = loadWavefrontDifferences(filename)
%LOAD2WAVEFRONTDIFFERENCES positions and difference output are in pixels

% extract raw numbers
fid_temp = fopen(filename,'r');
fgets(fid_temp);
fgets(fid_temp);
data = fscanf(fid_temp, '%f\t');
fclose(fid_temp);

% order data so each row represents single subaperture
data = reshape(data, [5, length(data) / 5])';

% sampled positions (reference centroids)
x = data(:,2);
y = data(:,3);

% wavefront differences
dx = (data(:,4) - data(:,2)); 
dy = (data(:,5) - data(:,3));
dx = dx - mean(dx);
dy = dy - mean(dy);

end