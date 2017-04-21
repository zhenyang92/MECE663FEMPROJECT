function [q_avg]=heat_flux(N_element,p)
% This subscript is written to evaluate the heat flux at the left boundary. 
% for the bilinear quadrilateral element, the heat flux term is of order 1,
% so the the order of Gaussian quadrature is (1+1)/2=1. For biquadratic
% element, The heat flux term is of order 3, so the order of Gaussian
% quadrature is (3+1)/2=2

B=1; % for surface integral
N_shape=(p+1)^2;
[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);
L_flux=0.1/N_y; % size of heat flux element
 a=zeros(N_y,2);
% find out elements at boundary 3
flux_B_element=zeros(1,N_y);
for j=1:N_y 
flux_B_element(j)= 1+(j-1)*N_x;
end
size_flux=size(flux_B_element);
N_flux=size_flux(2); % number of dirichlet elements

if p==1 % order of Lagragian polynomials
    q=1; % order of Gaussian quadrature
    [quad_points_flux,w_flux,N_quad_flux] = quadrature(1,B,0);   
    % find out nodes at boundary 3
     for e=1:N_y
     a(e,1:2)=[c_matrix(flux_B_element(e),1),c_matrix(flux_B_element(e),4)];
     end
     flux_B_nodes=unique(a);
   
    
elseif p==2
    q=2;
    [quad_points_flux,w_flux,N_quad_flux] = quadrature(2,B,1); 
    % find out nodes at boundary 3
    for e=1:N_y
    a(e,1:3)=[c_matrix(flux_B_element(e),1),c_matrix(flux_B_element(e),4),c_matrix(flux_B_element(e),7)];
    end
    flux_B_nodes=unique(a);
    
    
end

% call Jacob.m to obatin Jacobian for the elements at boundary 3
J_flux=zeros(2,N_quad_flux*N_flux);
Jinv=zeros(2,N_quad_flux*N_flux);
if p==1
for e=1:N_flux
    [J_flux(:,(e-1)*2+1:e*2)]=Jacob(1,N_element,flux_B_element(e),1,B,0);
    Jinv(:,(e-1)*2+1:e*2)=inv(J_flux(:,(e-1)*2+1:e*2));
end
elseif p==2
    for e=1:N_flux
    [J_flux(:,(e-1)*4+1:e*4)]=Jacob(2,N_element,flux_B_element(e),2,B,1);
    for i=1:2 % loop over quadrature points
      Jinv(:,(e-1)*4+2*(i-1)+1:(e-1)*4+2*(i-1)+2)=inv(J_flux(:,(e-1)*4+2*(i-1)+1:(e-1)*4+2*(i-1)+2));
    end
    end
end
    
% call for shape_function.m to obtain the gradient of shape functions
if p==1
[Psi_flux,dPsidz_flux,dPside_flux]=shape_function(1,p,B,0);
elseif p==2
[Psi_flux,dPsidz_flux,dPside_flux]=shape_function(2,p,B,1);   
end
q=zeros(1,N_flux);
% 
[u_sys]=FEM_project(N_element,p,5e3);
if p==1
for e=1:N_flux % loop over heat flux elements
    for j=1:4 % loop over shape functions
        q(e)=q(e)+(-1)*(15*(L_flux/2)*Jinv(1,2*(e-1)+1)*dPsidz_flux(j)*u_sys(c_matrix(flux_B_element(e),j))*w_flux(2))+15*(L_flux/2)*Jinv(1,2*(e-1)+2)*dPside_flux(j)*u_sys(c_matrix(flux_B_element(e),j))*w_flux(2);
    end
end
elseif p==2
for e=1:N_flux % loop over heat flux elements
    for j=1:9 % loop over shape functions
        for i=1:2 % loop over quadrature points
        q(e)=q(e)+(-1)*(15)*(L_flux/2)*Jinv(1,(e-1)*4+2*(i-1)+1)*(dPsidz_flux(j,i))*u_sys(c_matrix(flux_B_element(e),j))*w_flux(1,2)+(-1)*(15)*(L_flux/2)*Jinv(1,(e-1)*4+2*(i-1)+2)*(dPside_flux(j,i))*u_sys(c_matrix(flux_B_element(e),j))*w_flux(1,2);
        end
    end
end
end
   q_total=sum(q); 
   q_avg=q_total/0.1;