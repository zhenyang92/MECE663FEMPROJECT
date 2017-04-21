function [J,Jdet,N_quad]=Jacob(quad_order,N_element,element_N,p,B,h)
%This subscript is written to obtain Jacobian and the determinant of
%Jacobian at each quadrature point in a element-wise range when given the
%order of Legendre-Gauss Quadrature, the total number of elements, the
%element number and order of shape function
if B==0
N_quad=quad_order^2;
elseif B==1
N_quad=quad_order;
end
[Psi,dPsidz,dPside]=shape_function(quad_order,p,B,h);
% [quad_points_rhs,w]=quadrature(quad_order);
[global_coordinate]=global_coordinates(N_element,element_N,p);
%  evaluate Jacobian
 Jdet=zeros(1,N_quad);   
 for j=1:N_quad % loop over quadrature points
J(1:2,2*j-1:2*j)=[dPsidz(:,j)';dPside(:,j)']*global_coordinate; %Jacobian for one element, at different quadrature points are arranged in a 2X(2*N_quad) matrix 
Jdet(1,j)=det(J(1:2,2*j-1:2*j));
 end
     