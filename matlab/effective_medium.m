% Compute the elastic effective medium of a polymonomineral

% Wipe the slate clean
clc;

% Define some necessary variables
phi = 0*pi/180;    % z1 rotation
theta = 0*pi/180;  % x-rotation
psi = 0*pi/180;    % z2-rotation

% standard deviations of the rotation angles 
sdphi = 20*pi/180;     
sdtheta = 20*pi/180; 
sdpsi = 20*pi/180;

npts = 1000;
volume_fraction = ones(2*npts, 1)./(2*npts);

% Define the stiffness tensor 
C = [   148.3183 81.2892 68.1977 0 0 0;... 
        81.2892 148.3183 68.1977 0 0 0;...
        68.1977 68.1977 159.5873 0 0 0;...
        0 0 0 31.5959 0 0;...
        0 0 0 0 31.5959 0;...
        0 0 0 0 0 33.5145];


%% 

phi = [phi + rand(npts, 1).*sdphi; phi - rand(npts, 1).*sdphi];
theta = [theta + rand(npts, 1).*sdtheta; theta - rand(npts, 1).*sdtheta];
psi = [psi + rand(npts, 1).*sdpsi; psi - rand(npts, 1).*sdpsi];


euler_set = zeros(2*npts, 3);
euler_set(:,1) = phi;
euler_set(:,2) = theta;
euler_set(:,3) = psi;



[V, R, H] = vrh_homogenization(C, euler_set, volume_fraction);

V = V.*(C ~=0);
R = R.*(C ~=0);
H = H.*(C ~=0);

% Confess your sins then go do it all over again