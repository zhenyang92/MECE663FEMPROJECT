function [B_sys]=RobinBC_k(N_element,p,B,h)
% Robin boundary terms can be saperated to two parts, one added to the k_matrix and the other added to the rhs.
% For Robin boundary elements, eta is equal to 1 in the master
% element.
[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);
hh=20; %W/(m^2 K)
if p==1
for e=1:N_element
[global_coordinate((e-1)*4+1:e*4,1:2)]=global_coordinates(N_element,e,p);
end
elseif p==2
for e=1:N_element
[global_coordinate((e-1)*9+1:e*9,1:2)]=global_coordinates(N_element,e,p);
end  
end
% for bilinear quadraliteral case, the surface integral added to k_matrix
% where the robin boundary condition is at is of order 4, so the order if
% gaussian quadrature is (4+1)/2=2.5=3. Also because the convective heat
% transfer only happens at the upper surface, which corresponding to eta=1
% in the master element. Calculations will be performed on only three
% quadrature points which are(-sqrt(3/5),1),(0,1),(sqrt(3/5),1).

% for biquadratic quadraliteral case, the order of the surface integral is
% of order 8, so the order of gaussian quadrature is (8+1)/2=4.5=5; In the
% master elemetns, these quadrature points are located on the edge eta=1.
% Calculations will be performed on five points, which are
if p==1 
    quad_order=3;
elseif p==2
    quad_order=5;
end
[Psi,dPsidz,dPside,N_quad,w,N_shape]=shape_function(quad_order,p,B,h);

    
%the robin boundary terms added to the local k matrix    
B=zeros(N_shape,N_shape);   
L_Robin=sqrt(((0.25-0.1)^2+(0.1-0.05)^2))/(N_x/2); % element size of Robin boundary elements
if p==1
    for i=3:4 % loop over the last two shape funcations
         for k=1:N_quad % loop over quadrature points along zeta
            for j=3:4 % loop over the last two shape functions
                B(i,j)=B(i,j)+Psi(i,k)*Psi(j,k)*(L_Robin/2*hh)*w(k,1);
            end
         end     
    end
elseif p==2 
     for i=7:9 % loop over the last three shape funcations
         for k=1:N_quad % loop over quadrature points along zeta
            for j=7:9 % loop over the last three shape functions
                B(i,j)=B(i,j)+Psi(i,k)*Psi(j,k)*(L_Robin/2*hh)*w(k,1);
            end
         end     
    end
end
% find out Robin boundary elements
Robin_B_element=zeros(1,N_x/2);
for j=1:N_x/2
    Robin_B_element(j)=N_element-N_x/2+j;
end
Robin_size=size(Robin_B_element);
N_Robin=Robin_size(2);
% Assemble the local Robin boundary matrix to a gloabal Robin boundary matrix
B_sys=zeros(N_nodes,N_nodes);
for e=Robin_B_element(1):Robin_B_element(N_Robin) % loop over Robin boundary elements
    for i=1:N_shape % loop over rows
        for j=1:N_shape % loop over columns
            B_sys(c_matrix(e,i),c_matrix(e,j))=B_sys(c_matrix(e,i),c_matrix(e,j))+B(i,j);
        end

    end
end    
   
    


