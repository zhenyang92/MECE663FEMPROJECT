function [node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p)
% meshing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subscribt is written to obtain the nodal matrix, conectivity matrix
% as well as the global coordinate of the nodes when given total number of
% nodes used and the number of the node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=sqrt(N_element/2)*ones(1,5);%%%%%%%% % number of segments along each edge % N(3)=N(4)=N(5)=N(1)=N(2)
L=[0.1,0.15,0.1,0.1,0.05];%%%%%%%%%%%% % length of each edge
h=[0.1,0.15,0.1,0.1,0.05]./N;%%%%%%%%% % size of elements along each edge
N_x=N(1)+N(2);%%%%%%%%%%%%%%%%%%%%%%%% % number of elements in x-direction
N_y=N(3);%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % number of elements in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p==1 % Bilinear quadraliteral shape functions
N_nodes=(N_x+1)*(N_y+1);%%%%%%%%%%%%%% % number of total nodes
% node matrix 
x_location=zeros(N_y+1,N_x+1);
for i=1:N_x
    if i <=N(1)
    x_location(:,i+1)=x_location(:,i)+h(1);
    else
    x_location(:,i+1)=x_location(:,i)+h(2);
    end
end

y_location=zeros(N_y+1,N_x+1);
h_y=zeros(1,N_y+1);

for i=1:N_y %loop over y_direction
    for j=1:N_x+1 %loop over x_direction
       if j <=N(1)+1
       h_y(j)=h(3);
       else 
           L_y(j)=(L(2)+L(1)-x_location(i,j))/L(2)*(L(4)-L(5))+L(5);
           h_y(j)=L_y(j)/N(3);
       end
           y_location(i+1,j)=y_location(i,j)+h_y(j);
    end
end

node_matrix=zeros(N_nodes,2);

for i=1:N_y+1 %loop over y-direction
    for j=1:N_x+1 %loop over x-dirction
        node_matrix(j+(N_x+1)*(i-1),1)=x_location(i,j);
        node_matrix(j+(N_x+1)*(i-1),2)=y_location(i,j);
    end
end
% connectivitiy matrix 
c_matrix=zeros(N_element,4);
e=1;
for j=1:N_y
    while e<=j*N_x
        c_matrix(e,1:2)=[e+j-1,e+j];
        e=e+1;
    end
end
e=1;
for j=1:N_y
    while e<=j*N_x
        c_matrix(e,3:4)=[e+N_x+j+1,e+N_x+j];
        e=e+1;
    end
end
elseif p==2 %BIquadratic quadraliteral shape functions now the element size is half of that of Bilinear case
    N_shape=(p+1)^2;
    N_nodes=(N_x*2+1)*(N_y*2+1);
    x_location=zeros(N_y*2+1,N_x*2+1);
for i=1:N_x*2
    if i <=N(1)*2
    x_location(:,i+1)=x_location(:,i)+h(1)/2;
    else
    x_location(:,i+1)=x_location(:,i)+h(2)/2;
    end
end
y_location=zeros(N_y*2+1,N_x*2+1);
h_y=zeros(1,N_y*2+1);

for i=1:N_y*2 %loop over y_direction
    for j=1:N_x*2+1 %loop over x_direction
       if j <=N(1)*2+1
       h_y(j)=h(3)/2;
       else 
           L_y(j)=(L(2)+L(1)-x_location(i,j))/L(2)*(L(4)-L(5))+L(5);
           h_y(j)=L_y(j)/N(3)/2;
       end
           y_location(i+1,j)=y_location(i,j)+h_y(j);
    end
end

node_matrix=zeros(N_nodes,2);

for i=1:N_y*2+1 %loop over y-direction
    for j=1:N_x*2+1 %loop over x-dirction
        node_matrix(j+(N_x*2+1)*(i-1),1)=x_location(i,j);
        node_matrix(j+(N_x*2+1)*(i-1),2)=y_location(i,j);
    end
end
% connectivitiy matrix 
c_matrix=zeros(N_element,N_shape);

for i=1:N_y % loop over elements in y-direction
    for j=1:N_x % loop over elements in x-direction
        c_matrix((i-1)*N_x+j,1:3)=[1+(i-1)*2*(2*N_x+1)+(j-1)*2,2+(i-1)*2*(2*N_x+1)+(j-1)*2,3+(i-1)*2*(2*N_x+1)+(j-1)*2];
        c_matrix((i-1)*N_x+j,4:6)=[1+(i-1)*2*(2*N_x+1)+N_x*2+1+(j-1)*2,2+(i-1)*2*(2*N_x+1)+N_x*2+1+(j-1)*2,3+(i-1)*2*(2*N_x+1)+N_x*2+1+(j-1)*2];
        c_matrix((i-1)*N_x+j,7:9)=[1+(i-1)*2*(2*N_x+1)+(N_x*2+1)*2+(j-1)*2,2+(i-1)*2*(2*N_x+1)+(N_x*2+1)*2+(j-1)*2,3+(i-1)*2*(2*N_x+1)+(N_x*2+1)*2+(j-1)*2];
    end
end
end   



