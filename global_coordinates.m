
function [global_coordinate]=global_coordinates(N_element,element_N,p)
[node_matrix,c_matrix,N_nodes,N_x,N_y,x_location,y_location]=meshing(N_element,p);
global_coordinate=[node_matrix(c_matrix(element_N,:),1),node_matrix(c_matrix(element_N,:),2)];