
function plotting(N_element,p,A)
[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);
[u_sys,u_sys_rearrange]=FEM_project(N_element,p,A);
if p==1
mesh_plot=zeros(N_y+1,N_x+1);
elseif p==2
mesh_plot=zeros(N_y*2+1,N_x*2+1);
end
N_shape=(p+1)^2;
x=zeros(N_shape,N_element);
y=zeros(N_shape,N_element);
z=zeros(N_shape,N_element);
for e=1:N_element
    for i=1:N_shape
        x(i,e)=node_matrix(c_matrix(e,i),1);
        y(i,e)=node_matrix(c_matrix(e,i),2);
        z(i,e)=u_sys(c_matrix(e,i));
    end
end

figure
mesh(x_location,y_location,mesh_plot)
figure
meshc(x_location,y_location,u_sys_rearrange)
colorbar
figure
patch(x,y,z)
colorbar

q_avg_linear=zeros(1,5);
q_avg_quadratic=zeros(1,5);
    N_element_flux=[2,8,32,128,512];
for e=1:5
    [q_avg_linear(e)]=heat_flux(N_element_flux(e),1);
    [q_avg_quadratic(e)]=heat_flux(N_element_flux(e),2);
end

figure
plot(N_element_flux,q_avg_linear,N_element_flux,q_avg_quadratic); 
hlegend=legend('linear Lagrangian polynomials','quadratic Lagrangian polynomials');
set(hlegend,'fontsize',15);
hxlabel=xlabel('Number of element');
set(hxlabel,'fontsize',17);
hylabel=ylabel('average heat flux (W/m^2)');
set(hylabel,'fontsize',17);


