
load(['Data'  num2str(cable_no ) '.mat'])

% the position for 'lumpedM' need to adjust here
node_num       =   1   +   ele_num*order;
ele_info        =   zeros(ele_num,(order+1)+1);
I               =   ones(1,ele_num);
ele_info        =   I'.*[1:order+1 1];   
I               =   ones(1,(order+1)+1);
ele_info        =   ele_info    +   order*I.*((1:ele_num)-1)';

lumpedM   =     zeros(node_num*3,node_num*3);
C   =   zeros(node_num*3,node_num*3);
tU  =   zeros(node_num*3,1); 
tV  =   zeros(node_num*3,1); 
tA  =   zeros(node_num*3,1);   
tF  =   zeros(node_num*3,1);
tR  =   zeros(node_num*3,1);
tN  =   zeros(node_num,1);
Fw =    zeros(node_num*3,1); 
Fex =    zeros(node_num*3,1);
rho_mass =  w/g/Area;


if single==1 
    
    SingleShape
    
else
    
    ViscousDampMethod

end
T0  =   H
Fex =   Fw;
tU  =   zeros(node_num*3,1); 
tV  =   zeros(node_num*3,1); 
tA  =   zeros(node_num*3,1);   

CableOutput

clear I i no Va xL zL single L  fid fileID float_point fixed_point