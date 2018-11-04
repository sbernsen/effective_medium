function [V, R, H] = vrh_homogenization(C, euler_set, fraction)
% VRH_HOMOGENIZATION compute the voigt, reuss and hill average from a set
% of crystal c-axis orientations 
%
% INPUT
%   C - the 6-by-6 stiffness matrix 
%   euler_set - the m-by-3 set of euler angles (phi, theta, psi)
%
% OUTPUT
%
%--------------------------------------------------------------------------

m = length(euler_set(:,1) );

% You must comply 
S = inv(C);

% Allocate space
V = zeros(6, 6);
R = V;
H = V;

for i = 1:m
   rot = euler_rotation(euler_set(i,1), euler_set(i,2), euler_set(i,3) );
   M = bond(rot);
   N = M';
   V = V + fraction(i).*( M*C*M');
   R = R + ( fraction(i) ).*(N*S*N');
   
end

R = inv(R);
H = (V + R)./2;


end


function [R] = euler_rotation(phi, theta, psi)
% EULER_ROTATION computes the z-x-z rotation matrix from the three euler 
% angles.
%
% INPUT
%   phi, theta, psi - the euler angles for each subsequent rotation input
%                     as radians
%--------------------------------------------------------------------------

D = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
C = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
B = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];

R = D*B*C;


end

function [M] = bond(R)
%BOND computes the 6-by-6 bond rotation matrix to be applied to the  


M = [  R(1,1) R(1,2) R(1,3) ...
                2*R(1,2)*R(1,3)	2*R(1,3)*R(1,1) 2*R(1,1)*R(1,2);...
            R(2,1) R(2,2) R(2,3) ...
                2*R(2,2)*R(2,3) 2*R(2,3)*R(2,1) 2*R(2,1)*R(2,2);...
            R(3,1) R(3,2) R(3,3)... 
                2*R(3,2)*R(3,3)	2*R(3,3)*R(3,1) 2*R(3,1)*R(3,2);...
            R(2,1)*R(3,1) R(2,2)*R(3,2) R(2,3)*R(3,3)...
                R(2,2)*R(3,3) + R(2,3)*R(3,2)...
                R(2,1)*R(3,3) + R(2,3)*R(3,1)...
                R(2,2)*R(3,1) + R(2,1)*R(3,2);...
            R(3,1)*R(1,1) R(3,2)*R(1,2) R(3,3)*R(1,3)...
                R(1,2)*R(3,3) + R(1,3)*R(3,2)... 
                R(1,3)*R(3,1) +	R(1,1)*R(3,3)... 
                R(1,1)*R(3,2) + R(1,2)*R(3,1);...
            R(1,1)*R(2,1) R(1,2)*R(2,2) R(1,3)*R(2,3)...
                R(1,2)*R(2,3) + R(1,3)*R(2,2)...
                R(1,3)*R(2,1) + R(1,1)*R(2,3)...
                R(1,1)*R(2,2) + R(1,2)*R(2,1) ];


end


