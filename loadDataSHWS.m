function [Mx, My, A] = loadDataSHWS(filename,pitchPixel,focal)
%LOADDATASHWS Numbers read in are raw centroid and reference locations

% wavefront sensor parameters
pitchPixel = 7.4e-6;        % 7.4-um pixels
focal = 10e-3;              % 10-mm focal lenslets

% extract raw numbers
fid_temp = fopen(filename,'r');
fgets(fid_temp);
fgets(fid_temp);
data = fscanf(fid_temp, '%f\t');
fclose(fid_temp);

% order data so each row represents single subaperture
data = reshape(data, [5, length(data) / 5])';

% Wavefront slopes vector
wx = -(data(:,4) - data(:,2)) * pitchPixel / focal;
wy = (data(:,5) - data(:,3)) * pitchPixel / focal;
wx = wx - mean(wx);
wy = wy - mean(wy);

% reference locations
x = round(data(:,2), 5);    % rounding for safety
%x = (x - mean(x)) * pitchPixel;

% reference location meshgrids
u = unique(x);
N = length(u);
A = false(N);
for idx = 1:N
    len = length(find(x == u(idx)));
    a = round(0.5 * (N - len)) + 1;
    b = a + len - 1;
    A(a:b,idx) = true;
end

% 'measured' wavefront slope meshgrids
Mx = zeros(N);
My = zeros(N);
k = 1;
for idx = 1:numel(A)
    if A(idx)
        Mx(idx) = wx(k);
        My(idx) = wy(k);
        k = k + 1;
    end
end

% data is index left-to-right then top-to-bottom
Mx = padarray(Mx', [1 1]);
My = padarray(My', [1 1]);
A = padarray(A, [1 1]);





% % wavefront sensor parameters
% pitchPixel = 6.4e-6;        % 7.4-um pixels
% % focal = 5.9041e-3;              % 10-mm focal lenslets
% padding_size=0;
% 
% % extract raw numbers
% fid_temp = fopen(filename,'r');
% fgets(fid_temp);
% fgets(fid_temp);
% data = fscanf(fid_temp, '%f\t');
% fclose(fid_temp);
% 
% % order data so each row represents single subaperture
% data = reshape(data, [5, length(data) / 5])';
% 
% % Wavefront slopes vector
% wx = -(data(:,4) - data(:,2)) * pitchPixel / focal;
% wy = -(data(:,5) - data(:,3)) * pitchPixel / focal;
% wx = wx - mean(wx);
% wy = wy - mean(wy);
% %pv = [min([wx;wy]) max([wx;wy])];
% 
% % reference locations
% x = round(data(:,2), 5);    % rounding for safety
% %x = data(:,2);
% x = (x - mean(x)) * pitchPixel;
% %             y = round(data(:,3), 5);
% %             y = (y - mean(y)) * pitchPixel;
% 
% % reference location meshgrids
% u = unique(x);
% N = length(u);
% A = false(N);
% for idx = 1:N
%     len = length(find(x == u(idx)));
%     a = round(0.5 * (N - len)) + 1;
%     b = a + len - 1;
%     A(a:b,idx) = true;
% end
% 
% % 'measured' wavefront slope meshgrids
% Mx = zeros(N);
% My = zeros(N);
% k = 1;
% for idx = 1:numel(A)
%     if A(idx)
%         Mx(idx) = wx(k);
%         My(idx) = wy(k);
%         k = k + 1;
%     end
% end
% 
% 
% 
% % data is index left-to-right then top-to-bottom
% Mx = padarray(Mx', [padding_size padding_size]);
% My = padarray(My', [padding_size padding_size]);
% A = padarray(A, [padding_size padding_size]);
