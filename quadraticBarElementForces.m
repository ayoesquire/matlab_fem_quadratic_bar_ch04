function y = quadraticBarElementForces(k,u)
% quadraticBarElementForces     This function returns the element nodal
%                               force vector given the elemenet stiffness
%                               matrix k and the element nodal displacement
%                               vector u.
y = k * u;                      