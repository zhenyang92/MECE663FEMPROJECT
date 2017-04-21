
function[k_m]=k_matrix(N_element,element_N,p,B,h)
% this subscript is wrriten to obtain the k matrix in a element-wise range
% when given the order of Legendre-Gauss Quadrature, 
% the total number of elements and the element number
% for bilinear quadraliteral case, the order of kij is 2, so the order of
% gaussian quadrature is (2+1)/2=1.5=2
% for biquadratic quadraliteral case, the order of kij is 6, so the order
% of gaussian quadrature is (6+1)/2=3.5=4
if p==1 
    quad_order=2;
elseif p==2
    quad_order=4;
end

[Psi,dPsidz,dPside,N_quad,w,N_shape]=shape_function(quad_order,p,B,h);
[J,Jdet]=Jacob(quad_order,N_element,element_N,p,B,h); % call Jacob.m to obtain the Jacobian and the determinant of the Jacobian
[quad_points,w_k] = quadrature(quad_order,B,h); % call quadrature.m to obtain the quadrature points and weight functions

kx=15; %W/(m k)
ky=5; %W/(m k)
 
dPsidxdy=zeros(2*N_quad,4);
dPsidx=zeros(N_quad,N_shape);
dPsidy=zeros(N_quad,N_shape);
for i=1:N_shape %loop over shape functions
    for j=1:N_quad % loop over quadrature points
    dPsidxdy(2*j-1:2*j,i)=(J(:,2*j-1:2*j))\[dPsidz(i,j);dPside(i,j)];
    dPsidx(j,i)=dPsidxdy(2*(j-1)+1,i);
    dPsidy(j,i)=dPsidxdy(2*j,i);
    end
end
% evaluating F
for k=1:N_quad
    for i=1:N_shape % loop over shape functions
       for j=1:N_shape % loop over shpae functions
           F(i,(k-1)*N_shape+j)=(kx*dPsidx(k,i)*dPsidx(k,j)+ky*dPsidy(k,i)*dPsidy(k,j))*Jdet(k);
        end
    end
    
end
% adding the elements F multiplied by weight functions in two directions 
k_m=zeros(N_shape,N_shape);
for i=1:N_shape %loop over shape functions
    for j = 1:N_shape % loop over shape functions
        for k = 1:N_quad % loop over quadrature points in one direction
        if p==1
    k_m(i,j)=[F(i,j),F(i,j+N_shape),F(i,j+N_shape*2),F(i,j+N_shape*3)]*(w_k(:,1).*w_k(:,2));
        elseif p==2
    k_m(i,j)=[F(i,j),F(i,j+N_shape),F(i,j+N_shape*2),F(i,j+N_shape*3),F(i,j+N_shape*4),F(i,j+N_shape*5),F(i,j+N_shape*6),F(i,j+N_shape*7),F(i,j+N_shape*8),F(i,j+N_shape*9),F(i,j+N_shape*10),F(i,j+N_shape*11),F(i,j+N_shape*12),F(i,j+N_shape*13),F(i,j+N_shape*14),F(i,j+N_shape*15)]*(w_k(:,1).*w_k(:,2));   
        end
        end
    end
end


     
        