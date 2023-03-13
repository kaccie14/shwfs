function [W, rms] = zonalIterativeFourier(Mx, My, A, h, n)
%ITERATIVEFOURIER reconstruct wavefront with iterative ffts
%   matrices Mx, My and A must be same size
%
%   Mx - measured x-slopes
%   My - measured y-slopes
%   A - aperture
%   h - sampling spacing (i.e., subaperture pitch)

% spatial frequency meshgrid
[u, v] = freqspace(size(A), 'meshgrid');
u2v2 = u.^2 + v.^2;

% Initial condition
Wx = Mx; % x-slopes in zero-padded meshgrid
Wy = My; % y-slopes in zero-padded meshgrid

% initialize reconstruction error tracking variables
rms = zeros(n + 1,1);
E = [Mx My];
e = E([A A]);
N = numel(e) / 2;
rms(1) = sqrt(e'*e) / N;

for idx = 1:n
    FWx = fftshift(fft2(Wx));
    FWy = fftshift(fft2(Wy));
    FW = (-1i * h / pi) * (u.*FWx + v.*FWy) ./ u2v2;
    o = floor(0.5 * size(FW)) + 1;          % origin
    FW(o(1), o(2)) = 0.0;                   % origin being 0 forces zero mean
    W = real(ifft2(ifftshift(FW)));
    
    % Compute gradient of estimated wavefront
    FWx = (1i * pi / h) * u .* FW;
    FWy = (1i * pi / h) * v .* FW;
    Wx = ifft2(ifftshift(FWx));
    Wy = ifft2(ifftshift(FWy));
    
    % Keep track of residual RMS
    E = [Wx - Mx, Wy - My];
    e = E([A A]);
    rms(idx + 1) = sqrt(e'*e) / N;
    
    % Put measurement (slope) signals back inside signal boundaries
    Wx(A) = Mx(A);
    Wy(A) = My(A);
end% for

W = A .* W;
end% iterativeFourier