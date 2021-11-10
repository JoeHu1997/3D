%% Phenomena considered include self-weight,linear internal damping,hydrodynamic force, internal tension and etc.

%% External Force
% A. Drag Force
% B. Froude-Krylov Force
% C. Hydrodynamic mass
% D. Internal damping force 
% E. Buoyancy

% Equivalent nodal force between two nodes
% F_node1   =   0.5*F2*Le+1/3*(F1-F2)*Le
% F_node2   =   0.5*F1*Le+1/3*(F2-F1)*Le

tR              =   zeros(node_num*3,1);
FD              =   zeros(node_num*3,1);
Ffk             =   zeros(node_num*3,1);
Fdamp           =   zeros(node_num*3,1);
[vT vN l rZ]    =   NodalVector(x,y,z,tU,node_num,Vorig,field_V);
[ele_cos u]     =   updateInfo(x,y,z,tU,ele_num,ele_info,order);
%% Force on node 
for     no =    1:node_num-1
%% A.Drag Force
    Fdn     =   -0.5*Cdn*rho_water*D*l(no)*norm(vN(no,:))*vN(no,:);
    Fdt     =   -0.5*Cdt*rho_water*D*l(no)*norm(vT(no,:))*vT(no,:);
    Fd      =   Fdn+Fdt;   
    
    FD(1+(no-1)*3:3+(no)*3)     =   FD(1+(no-1)*3:3+(no)*3)  +   0.5*[Fd Fd]';
    
%% B. Froude-Krylov Force
%    FK  =   rho_water*Area*l(no)*(field_A(1+(no)*3:3+(no)*3)+field_A(1+(no-1)*3:3+(no-1)*3))/2;
%    Ffk(1+(no-1)*3:3+(no)*3)    =   Ffk(1+(no-1)*3:3+(no)*3) +  0.5*[FK FK];
%% C. Hydrodynamic mass
%    ma   =   Ca*Area*l(no);    
%    Ae  =   [rX;rY;rZ;];
    
%% D. Internal damping force
    if (u(no,end)-ele_he(no))   >   10^-3
        Damp    =   Cc*dot(tV(1+(no)*3:3+(no)*3)-tV(1+(no-1)*3:3+(no-1)*3),rZ(no,:))*rZ(no,:);
        Fdamp(1+(no-1)*3:3+(no)*3)   =  Fdamp(1+(no-1)*3:3+(no)*3)  +   [-Damp Damp]';
    end
end
clear Fdn Fdt FK Damp Fd

%% E. Buoyancy
% Buoyancy has been calculated at Initialize.m.
% Fb: buoyancy of cable


%% F. Weight 
% Weight has been calculated at Initialize.m.
% Fw: buoyancy of cable
tR  =   tR  +   FD  +   Ffk +   Fdamp   +Fex;

