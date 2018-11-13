% Compute the elastic effective medium of a polymonomineral
tic
% Wipe the slate clean
clc;

%% Define some necessary variables to create a synthetic crystal fabric
phi = 0*pi/180;    % z1 rotation
theta = 0*pi/180;  % x-rotation
psi = 0*pi/180;    % z2-rotation

% standard deviations of the rotation angles 
sdphi = 40*pi/180;     
sdtheta = 25*pi/180; 
sdpsi = 40*pi/180;

% define how many points to use for the homogenization and each crystals contribution to the volume
npts = 100000;
volume_fraction = ones(npts, 1)./(npts);

% Define the stiffness tensor 
C = [   148.3183 81.2892 68.1977 0 0 0;... 
        81.2892 148.3183 68.1977 0 0 0;...
        68.1977 68.1977 159.5873 0 0 0;...
        0 0 0 31.5959 0 0;...
        0 0 0 0 31.5959 0;...
        0 0 0 0 0 33.5145];


%% Synthesize the set of euler angles for a z-x-z
rnd1 = rand(npts,1);
rnd1 = -(1-rnd1) + rnd1;

rnd2 = rand(npts,1);
rnd2 = -(1-rnd2) + rnd2;

rnd3 = rand(npts,1);
rnd3 = -(1-rnd3) + rnd3;

phi = phi + rnd1.*sdphi;
theta = theta + rnd2.*sdtheta;
psi = psi + rnd3.*sdpsi;


euler_set = zeros(npts, 3);
euler_set(:,1) = phi;
euler_set(:,2) = theta;
euler_set(:,3) = psi;


%% Homogenize

[V, R, H] = vrh_homogenization(C, euler_set, volume_fraction);

% The material is now orthorhombic so we only need the upper left 3-by-3 and the lower right diagonal
V = V.*(C ~=0);
R = R.*(C ~=0);
H = H.*(C ~=0);
toc

