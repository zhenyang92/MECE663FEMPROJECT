function[B_rhs_sys]=RobinBC_rhs(N_element,p,B,h)
% Robin boundary terms can be saperated to two parts, one added to the k_matrix and the other added to the rhs.
% For Robin boundary elements, eta is equal to 1 in the master
% element.Since the Robin boundary term added to the rhs is of order 2, the order
% of Gaussian quadrature is (2+1)/2=1.5=2. Calculations will be performed on
% only 2 quadrature points, which are(-sqrt(1/3),1),(sqrt(1/3),1)
[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);
hh=20; %W/(m^2 K)
T_inf=300;
% for bilinear quadraliteral case, the surface integral added to rhs matrix
% where the robin boundary condition is at is of order 2, so the order if
% gaussian quadrature is (2+1)/2=1.5=2. Also because the convective heat
% transfer only happens at the upper surface, which corresponding to eta=1
% in the master element. Calculations will be performed on only two
% quadrature points which are(-sqrt(1/3),1),(sqrt(1/3),1).

% for biquadratic quadraliteral case, the order of the surface integral is
% of order 4, so the order of gaussian quadrature is (4+1)/2=2.5=3; In the
% master elemetns, these quadrature points are located on the edge eta=1.
% Calculations will be performed on five points, which are
% (-sqrt(3/5),1),(0,1),(sqrt(3/5),1).
if p==1 
    quad_order=2;
elseif p==2
    quad_order=3;
end
[Psi,dPsidz,dPside,N_quad,w,N_shape]=shape_function(quad_order,p,B,h);
% The robin boundary terms added to the rhs matrix   
B_rhs=zeros(N_shape,1);   
L_Robin=sqrt(((0.25-0.1)^2+(0.1-0.05)^2))/(N_x/2); % element size of Robin boundary elements
if p==1
    for i=3:4 % loop over the last two shape funcations
         for k=1:N_quad % loop over quadrature points along zeta
            B_rhs(i)=B_rhs(i)+Psi(i,k)*T_inf*(L_Robin/2*hh)*w(k);
         end
    end
elseif p==2
     for i=7:9 % loop over the last three shape funcations
         for k=1:N_quad % loop over quadrature points along zeta
            B_rhs(i)=B_rhs(i)+Psi(i,k)*T_inf*(L_Robin/2*hh)*w(k);
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

% Assemble the local rhs Robin boundary matrix to a gloabal rhs Robin boundary matrix
B_rhs_sys=zeros(N_nodes,1);
for e=Robin_B_element(1):Robin_B_element(N_Robin) % loop over Robin boundary elements
    for i=1:N_shape % loop over rows
            B_rhs_sys(c_matrix(e,i))= B_rhs_sys(c_matrix(e,i))+B_rhs(i);
    end
end    
   
    

    
    