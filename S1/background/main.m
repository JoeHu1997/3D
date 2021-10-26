%% This file is main program of the cable.
cd (root)
% call variable
%%
% Update the position of the floating end of the cable
% Obtain the motion of cable at last time step  
if  FOAMcall     <   MATLABcall 
    tU           =   U{MATLABcall}; 
    tV           =   V{MATLABcall}; 
    tA           =   A{MATLABcall};
    tF           =   F{MATLABcall};
    tN           =   N{MATLABcall};
    Uorig           =   tU;
    Aorig           =   tA;
end


float_pt_pos(1)	=	xfloat;
float_pt_pos(2)	=	yfloat;
float_pt_pos(3)	=	zfloat;
if  pt(3,1)==0
    tU(float_pt_node(1))    =   float_pt_pos(1)-x(float_pt_node(3)/3);
    tU(float_pt_node(2))    =   float_pt_pos(2)-y(float_pt_node(3)/3);
    tU(float_pt_node(3))    =   float_pt_pos(3)-z(float_pt_node(3)/3);
end



%%
% Calculate the external force actig on cable
cd (root)
cd ([root '/matlab'])
readField    
% Caculate force from fluid field
cd ([root '/matlab'])
CalculateExternalForce


%% Newton-Raphson iteration
NewtonRaphson


%%    
% Update the new position of the cable
X_tmp(:,1)=x+tU(1:3:end);
X_tmp(:,2)=y+tU(2:3:end);
X_tmp(:,3)=z+tU(3:3:end);
cd ([root '/matlab'])
probesLocation

%%
% Return to OpenFOAM
RF      =   tF-tR;
Fx      =   -1*RF(float_pt_node(1));
Fy      =   -1*RF(float_pt_node(2));
Fz      =   -1*RF(float_pt_node(3))+Adjust;


%%
% Save the data
cd (root) 
U{MATLABcall+1}   =   tU;
V{MATLABcall+1}   =   tV;
A{MATLABcall+1}   =   tA;
N{MATLABcall+1}   =   tN;
F{MATLABcall+1}   =   tF;
R{MATLABcall+1}   =   tR;
dis{MATLABcall}   =    u;
Force{MATLABcall} =   [Fx Fy Fz];    
X{MATLABcall+1}   =    X_tmp; 
save(['CableData' num2str(cable_no) '.mat'],'U','V','A','N','F','R','X','Force')
fclose('all');
