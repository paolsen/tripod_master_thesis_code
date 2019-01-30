function [M] = GenRot(phi,theta,psi)%z,y,x
M =      [[ cos(theta)*cos(phi)     -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)    sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)              0]
          [ sin(phi)*cos(theta)       cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)     -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)           0]
          [ -sin(theta)                            cos(theta)*sin(psi)                                 cos(theta)*cos(phi)                         0]
          [ 0                                             0                                                         0                              1]];
end

