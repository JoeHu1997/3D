%OpenFOAM need to return no of catenary to MATLAB

close all
clear
clc  
root = pwd;

cd ([root '/matlab'])

cable_no =1;

%CableShape

load(['moordata'    num2str(cable_no) '.mat']);
% - read position and node of cable
[x y z order node_num ele_info ele_he ele_num ele_X]=data(cable_no);
% - read BCs of cable and assemble into global setting
BCs


%% Establish matrix such as [B_L], [B_NL], [M], [C]
% local stiffness without Jacabin
shape
% k=[BL]T[BL]=[BNL]T[BNL]
k  = local_k(shape_fuc,order); 
[a0    a1  a2  a3  a4  a5  a6  a7] =   paraNewmark(a,b,deltaT);

% establish mass matrix
if Ca==0
    M0 =     mass(rho_mass,0,Area,Ca,order,ele_num,ele_info,ele_he);
    Ma =     zeros(node_num*3,node_num*3);
else
    M0 =     mass(rho_mass,rho_water,Area,Ca,order,ele_num,ele_info,ele_he);
    Ma =     zeros(node_num*3,node_num*3);
end
M0  =   M0  +   lumpedM;

%% Initialize global variable
U{1}               =   tU;
V{1}               =   tV;
A{1}               =   tA;
F{1}               =   tF;     % internal force
R{1}               =   tR;     % external force 
N{1}               =   tN;
X{1}               =   X_tmp;
%N0                 =   tN;
cd (root)

call        =   0;
FOAMcall    =   1; 
MATLABcall  =   0;     % the time execute Matlab  MATLABcall
xfloat      =   float_pt_pos(1);
yfloat      =   float_pt_pos(2);
zfloat      =   float_pt_pos(3);
Fx          =   0;
Fy          =   0;
Fz          =   0;

if  pt(3,1)==1
    Ux=0;Uy=0;Uz=0;
else
    Ux          =   tV(float_pt_node(1));
    Uy          =   tV(float_pt_node(2));
    Uz          =   tV(float_pt_node(3));
end
First       =   true;
Uorig           =   tU;
Vorig           =   tV;
Aorig           =   tA;

%main MATLABcall-1
cd (root)
