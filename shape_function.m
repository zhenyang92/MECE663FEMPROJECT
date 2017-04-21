function [Psi,dPsidz,dPside,N_quad,w,N_shape]=shape_function(quad_order,p,B,h)
% This subsript is written to obtaian the values of quadraliteral shape functions and
% the gradient of the shape functions when given the Gaussian quadrature
% order and the shape function order
% quad_order=2;
      [quad_points,w,N_quad]=quadrature(quad_order,B,h); % quadrature points 
       N_shape=(p+1)^2;% number of shape functions
%        N_quad=quad_order^2;
       Psi=zeros(N_shape,N_quad);
       dPsidz=zeros(N_shape,N_quad);
       dPside=zeros(N_shape,N_quad);
if p==1 % Bilinear quadraliteral shape functions
    for j=1:N_quad % loop over quadrature points
        Psi(1,j)=1/4*(1-quad_points(j,1))*(1-quad_points(j,2));
        Psi(2,j)=1/4*(1+quad_points(j,1))*(1-quad_points(j,2));
        Psi(3,j)=1/4*(1+quad_points(j,1))*(1+quad_points(j,2));
        Psi(4,j)=1/4*(1-quad_points(j,1))*(1+quad_points(j,2));
        % gradient of shape functions
        dPsidz(1,j)=-1/4*(1-quad_points(j,2));
        dPsidz(2,j)=1/4*(1-quad_points(j,2));
        dPsidz(3,j)=1/4*(1+quad_points(j,2));
        dPsidz(4,j)=-1/4*(1+quad_points(j,2));
        dPside(1,j)=-1/4*(1-quad_points(j,1));
        dPside(2,j)=-1/4*(1+quad_points(j,1));
        dPside(3,j)=1/4*(1+quad_points(j,1));
        dPside(4,j)=1/4*(1-quad_points(j,1));
    end
    
elseif p==2 % BIquadratic quadraliteral shape functions
    for j=1:N_quad
        Psi(1,j)=1/4*(1-quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,1)*quad_points(j,2);
        Psi(3,j)=-1/4*(1+quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,1)*quad_points(j,2);
        Psi(9,j)=1/4*(1+quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,1)*quad_points(j,2);
        Psi(7,j)=-1/4*(1-quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,1)*quad_points(j,2);
        Psi(2,j)=-1/2*(1-quad_points(j,1)^2)*(1-quad_points(j,2))*quad_points(j,2);
        Psi(8,j)=1/2*(1-quad_points(j,1)^2)*(1+quad_points(j,2))*quad_points(j,2);
        Psi(6,j)=1/2*(1-quad_points(j,2)^2)*(1+quad_points(j,1))*quad_points(j,1);
        Psi(4,j)=-1/2*(1-quad_points(j,2)^2)*(1-quad_points(j,1))*quad_points(j,1);
        Psi(5,j)=(1-quad_points(j,1)^2)*(1-quad_points(j,2)^2);
        % gradient of shape functions
        dPsidz(1,j)=-1/4*(1-quad_points(j,2))*quad_points(j,1)*quad_points(j,2)+1/4*(1-quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,2);
        dPsidz(3,j)=-1/4*(1-quad_points(j,2))*quad_points(j,1)*quad_points(j,2)-1/4*(1+quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,2);                 
        dPsidz(9,j)=1/4*(1+quad_points(j,2))*quad_points(j,1)*quad_points(j,2)+1/4*(1+quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,2);                                       
        dPsidz(7,j)=1/4*(1+quad_points(j,2))*quad_points(j,1)*quad_points(j,2)-1/4*(1-quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,2);                            
        dPsidz(2,j)=quad_points(j,1)*(1-quad_points(j,2))*quad_points(j,2);
        dPsidz(8,j)=-quad_points(j,1)*(1+quad_points(j,2))*quad_points(j,2);
        dPsidz(6,j)=1/2*(1-quad_points(j,2)^2)*quad_points(j,1)+1/2*(1-quad_points(j,2)^2)*(1+quad_points(j,1));
        dPsidz(4,j)=1/2*(1-quad_points(j,2)^2)*quad_points(j,1)-1/2*(1-quad_points(j,2)^2)*(1-quad_points(j,1));
        dPsidz(5,j)=-2*quad_points(j,1)*(1-quad_points(j,2)^2);
        
        dPside(1,j)=-1/4*(1-quad_points(j,1))*quad_points(j,1)*quad_points(j,2)+1/4*(1-quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,1);
        dPside(3,j)=1/4*(1+quad_points(j,1))*quad_points(j,1)*quad_points(j,2)-1/4*(1+quad_points(j,1))*(1-quad_points(j,2))*quad_points(j,1);                 
        dPside(9,j)=1/4*(1+quad_points(j,1))*quad_points(j,1)*quad_points(j,2)+1/4*(1+quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,1);                                       
        dPside(7,j)=-1/4*(1-quad_points(j,1))*quad_points(j,1)*quad_points(j,2)-1/4*(1-quad_points(j,1))*(1+quad_points(j,2))*quad_points(j,1); 
        dPside(2,j)=1/2*(1-quad_points(j,1)^2)*quad_points(j,2)-1/2*(1-quad_points(j,1)^2)*(1-quad_points(j,2));
        dPside(8,j)=1/2*(1-quad_points(j,1)^2)*quad_points(j,2)+1/2*(1-quad_points(j,1)^2)*(1+quad_points(j,2));
        dPside(6,j)=-quad_points(j,2)*(1+quad_points(j,1))*quad_points(j,1);
        dPside(4,j)=quad_points(j,2)*(1-quad_points(j,1))*quad_points(j,1);       
        dPside(5,j)=-2*quad_points(j,2)*(1-quad_points(j,1)^2); 
    end
end
