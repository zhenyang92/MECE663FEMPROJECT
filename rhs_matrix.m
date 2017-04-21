function [f]=rhs_matrix(N_element,element_N,A,p,B,h)
% This subscribt is written to evaluate rhs using Legendre-Gauss Quadratur
% method in the element_wise when given given the order of Legendre-Gauss Quadrature, 
% the total number of elements and the element number
% for bilinear quadraliteral case, the order of rhs term is 6, so the order of
% gaussian quadrature is (6+1)/2=3.5=4
% for biquadratic quadraliteral case, the order of kij is 6, so the order
% of gaussian quadrature is (6+1)/2=3.5=7
if p==1 
    quad_order=4;
elseif p==2
    quad_order=7;
end

[Psi,dPsidz,dPside,N_quad,w,N_shape]=shape_function(quad_order,p,B,h);
[J,Jdet]=Jacob(quad_order,N_element,element_N,p,B,h); % call Jacob.m to obtain the Jacobian and the determinant of the Jacobian
[quad_points,w_rhs] = quadrature(quad_order,B,h); % call quadrature.m to obtain the quadrature points and weight functions

% mapping the coordinates from the master element to the real element for
% each real element
transform_coordinate=zeros(N_quad,2);
a=zeros(N_quad,1);
[global_coordinate]=global_coordinates(N_element,element_N,p);

% evaluate a(x,y)
    for j=1:N_quad % loop over quadrature points
        transform_coordinate(j,1:2)=[sum(global_coordinate(:,1).*Psi(:,j)),sum(global_coordinate(:,2).*Psi(:,j))];
         a(j)=A*transform_coordinate(j,1).^2;
    end

% evaluate F
    F=zeros(N_shape,N_quad);
    for j=1:N_quad %loop over quadrature point
        for k=1:N_shape %loop over shape functions
            F(k,j)=a(j)*Psi(k,j)*Jdet(j);
        end
    end
    
% evaluate rhs using numerical integration
f=zeros(N_shape,1);
 for k=1:N_shape % loop over shpae functions
 f(k)=F(k,:)*(w_rhs(:,1).*w_rhs(:,2));
 end


