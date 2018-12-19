function y = quadraticBarElementStresses(k,u,A)
% quadraticBarElementStresses       This function return the element nodal
%                                   stress vector given the element
%                                   stiffness matrix k, the element nodal
%                                   displacement vector u, and the
%                                   cross-sectional area A.
y = k * u/A;