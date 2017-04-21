
function[u_sys,u_sys_rearrange]=FEM_project(N_element,p,A)
% The main function to obtain the solution to my FEM project, consisting of
% several subscripts, when given total number of elements, order of shape
% function and value of A.

% call for mesh function to obatin conectivity matrix and nodal matrix.
N_shape=(p+1)^2;

[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);

% call rhs_matrix functions to obtain local rhs matrices
f_sys=zeros(N_nodes,1);
f=zeros(N_shape,N_element);
for e=1:N_element
[f(:,e)]=rhs_matrix(N_element,e,A,p,0,0);
end
% call k_matrix functions to obtain local k matrices
q_k=2;
% [quad_points_k,w_k,N_quad_k] = quadrature(q_k);
k_sys=zeros(N_nodes,N_nodes);
for e=1:N_element
    [k(1:N_shape,N_shape*(e-1)+1:N_shape*e)]=k_matrix(N_element,e,p,0,0);
%     [J,Jdet,Psi,dPsidz,dPside,N_quad,w]=Jacob(2,N_element,e,p);
end

% Assemble the local k_matrix and rhs_matrix
for e=1:N_element % loop over elements
    for i=1:N_shape % loop over rows
        for j=1:N_shape % loop over columns
            k_sys(c_matrix(e,i),c_matrix(e,j))=k_sys(c_matrix(e,i),c_matrix(e,j))+k(j,N_shape*(e-1)+i);
        end
            f_sys(c_matrix(e,i))= f_sys(c_matrix(e,i))+f(i,e);
    end
end

% apply boundary conditions 
% find out Dirichlet boundary elements
Dirichlet_B_element=zeros(1,N_y);
for j=1:N_y 
Dirichlet_B_element(j)= 1+(j-1)*N_x;
end
size_Dirichlet=size(Dirichlet_B_element);
N_Dirichlet=size_Dirichlet(2); % number of dirichlet elements
% size_Dirichlet=size_Dirichlet
% find out Dirichlet boundary nodes
if p==1
for e=1:N_y
a(e,1:2)=[c_matrix(Dirichlet_B_element(e),1),c_matrix(Dirichlet_B_element(e),4)];
end
elseif p==2
for e=1:N_y
a(e,1:3)=[c_matrix(Dirichlet_B_element(e),1),c_matrix(Dirichlet_B_element(e),4),c_matrix(Dirichlet_B_element(e),7)];
end
end
Dirichlet_B_nodes=unique(a);
% 
% %  put Dirichlet boundary conditions onto the Global matrices k_sys and f_sys
if p==1
for e=1:N_y+1% loop over Dirichlet boundary elements 
k_sys(Dirichlet_B_nodes(e),:)=0;
k_sys(Dirichlet_B_nodes(e),Dirichlet_B_nodes(e))=1;
f_sys(Dirichlet_B_nodes(e))=400;
end
elseif p==2
for e=1:N_y*2+1% loop over Dirichlet boundary elements 
k_sys(Dirichlet_B_nodes(e),:)=0;
k_sys(Dirichlet_B_nodes(e),Dirichlet_B_nodes(e))=1;
f_sys(Dirichlet_B_nodes(e))=400;
end
end
% 
% call function RobinBC_k to obtain a global Robin boundary matrix to be
% added to the global k matrix
[B_sys]=RobinBC_k(N_element,p,1,0);
% call function RobinBC_rhs to obtain a global Robin bounday matrix to be
% added to the global rhs f matrix;
[B_rhs_sys]=RobinBC_rhs(N_element,p,1,0);

k_sys=k_sys+B_sys;
f_sys=f_sys+B_rhs_sys;
u_sys=k_sys\f_sys;
% 
% rearrange the global coordinate and the solution for the purpose of
% plotting
if p==1
u_sys_rearrange=zeros(N_y+1,N_x+1);
for j=1:N_y+1 % loop over rows
    for i=1:N_x+1 % loop over columns
        u_sys_rearrange(j,i)=u_sys((j-1)*(N_x+1)+i);
    end
end
elseif p==2
u_sys_rearrange=zeros(N_y*2+1,N_x*2+1); 
for j=1:N_y*2+1 % loop over rows
    for i=1:N_x*2+1 % loop over columns
        u_sys_rearrange(j,i)=u_sys((j-1)*(N_x*2+1)+i);
    end
end
end

