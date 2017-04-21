function [quad_points,w,N_quad] = quadrature(q,B,h)
% this subscript is written to obtain the quadrature points locations and
% the corresponding weight function values when given the order of Legendre-Gauss Quadrature
% B=0 stand for inner elements, B=1 stand for surface integral for Robin
% boundary conditions and for studying heat flux 
% B=1;
% q=22;
%% inner elements quadrature points
if B==0 % quadrature points for inner elements 
    if h==0
    N_quad=q^2;
if q==2
    zeta=[-1/sqrt(3);1/sqrt(3)];
    eta=[-1/sqrt(3);1/sqrt(3)];
    wzeta=[1;1];
    weta=[1;1];
elseif q==4
    zeta=[-0.861136;-0.339981;0.339981;0.861136];
    eta=[-0.861136;-0.339981;0.339981;0.861136];
    wzeta=[0.347855;0.652145;0.652145;0.347855];
    weta=[0.347855;0.652145;0.652145;0.347855];
elseif q==7
    zeta=[-0.9491079123427585;-0.741531185599;-0.40584515137739;0;0.40584515137739;0.741531185599;0.9491079123427585];
    eta=[-0.9491079123427585;-0.741531185599;-0.40584515137739;0;0.40584515137739;0.741531185599;0.9491079123427585];
    wzeta=[0.1294849661688697;0.2797053914892766;0.3818300505051189;0.4179591836734694;0.3818300505051189;0.2797053914892766;0.1294849661688697];
    weta=[0.1294849661688697;0.2797053914892766;0.3818300505051189;0.4179591836734694;0.3818300505051189;0.2797053914892766;0.1294849661688697];        
end
[A,B] = meshgrid(zeta,eta);
c=cat(2,A',B');
quad_points=reshape(c,[],2);
[C,D] = meshgrid(wzeta,weta);
d=cat(2,C',D');
w=reshape(d,[],2);
    end
%% boundary elements quadrature points for surface integral (Robin boundary condition, heat flux in convergence study)
elseif B==1 
          if q==1 % these quadurature points are for the heat flux at boundary 3 only, located on zeta=-1
          zeta=-1;
          eta=0;
          wzeta=0;
          weta=2;
          N_quad=q;
          elseif q==2
              if h ==0 % these quadurature points are for the robin boundary elements at boundary 2, added to rhs matrix
          zeta=[-sqrt(1/3),sqrt(1/3)]; % quadrature point locations along zeta
          eta=[1,1];
          wzeta=[1;1];
          weta=[0;0];
          N_quad=q;
              elseif h==1 % these quaduratre points are for heat flux at boundary 3, located on zeta=-1
          zeta=[-1,-1];    
          eta=[-sqrt(1/3),sqrt(1/3)];   
          wzeta=[0;0];    
          weta=[1;1];
          N_quad=2;
              end
          elseif q==3  % these quadurature points are for the robin boundary elements at boundary 2, added to k_matrix
          zeta=[-sqrt(3/5),0,sqrt(3/5)]; % quadrature point locations along zeta
          eta=[1,1,1];
          wzeta=[5/9;8/9;5/9];
          weta=[0;0;0];
          N_quad=q;
          elseif q==5
          zeta=[-0.9061798459386640,-0.5384693101056831,0,0.5384693101056831,0.9061798459386640]; % quadrature point locations along zeta
          eta=[1,1,1,1,1];
          wzeta=[0.2369268850561891;0.4786286704993665;0.5688888888888889;0.4786286704993665;0.2369268850561891];
          weta=[0;0;0;0;0];
          N_quad=q;
          end
          quad_points=[zeta',eta'];
          w=[wzeta,weta];
end


