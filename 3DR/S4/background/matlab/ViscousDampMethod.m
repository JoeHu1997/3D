

%% Predict PreTension
xL          =   ((float_point(1:2)-fixed_point(1:2))*(float_point(1:2)-fixed_point(1:2))')^0.5;
zL          =   float_point(3)-fixed_point(3);
Cos         =   (float_point(1:2)-fixed_point(1:2))/xL;
Cos         =   abs(Cos);

fix_pt_pos      =   [0.00  0.0     0.00];
float_pt_pos    =   [L     0.0     0.00];

syms    H  Va
f1  =   xL   ==  H/EA*L+H/w*(asinh((w*L+Va)/H)-asinh(Va/H));
f2  =   zL   ==  1/EA*(0.5*w*L^2+Va*L)+H/w*((1+((w*L+Va)/H)^2)^0.5-(1+(Va/H)^2)^0.5);
S   =   vpasolve([f1 f2]);
H   =   double(S.H);
Va  =   double(S.Va);
T0  =   abs(H);
clear f1 f2  S

%% Initialize for the shape data

x                =  linspace(0,float_pt_pos(1),node_num)';
y                =  zeros(node_num,1);
z                =  linspace(0,float_pt_pos(3),node_num)';

ele_he      =   zeros(ele_num,1);
ele_X       =   zeros(ele_num,order+1);
for     i   =   1:ele_num
    ele_he(i)   =   ((x(ele_info(i,1))-x(ele_info(i,order+1)))^2+(y(ele_info(i,1))-y(ele_info(i,order+1)))^2+(z(ele_info(i,1))-z(ele_info(i,order+1)))^2)^0.5;
    Pos             =   zeros(order+1,3);
    Pos(:,1)        =   x(ele_info(i,1:order+1));
    Pos(:,2)        =   y(ele_info(i,1:order+1));
    Pos(:,3)        =   z(ele_info(i,1:order+1));
    Pos             =   Pos-Pos(1,:);
    Pos             =   Pos'.^2;
    ele_X(i,:)      =   sum(Pos).^0.5;
end

   
fix_pt_dir  =   [1   1   1];
float_pt_dir  =   [1   1   1];
fix_pt_node=[];
float_pt_node=[];

dir1=find(fix_pt_dir(:));
loc1=[find(abs(x-fix_pt_pos(1)) < 0.001); find(abs(y-fix_pt_pos(2)) < 0.001); find(abs(z-fix_pt_pos(3)) < 0.001);];
fix_pt_node=[fix_pt_node 3*(mode(loc1)-1)+dir1;];
dir2=find(float_pt_dir(:));
loc2=[find(abs(x-float_pt_pos(1)) < 0.001); find(abs(y-float_pt_pos(2)) < 0.001); find(abs(z-float_pt_pos(3)) < 0.001);];
float_pt_node=[float_pt_node 3*(mode(loc2)-1)+dir2;];


clear I Pos loc1 dir1 loc2 dir2



%% ViscousDampMethod
shape
k  =    local_k(shape_fuc,order); 
[a0    a1  a2  a3  a4  a5  a6  a7] =   paraNewmark(a,b,deltaT);
for i=1:ele_num
        Fw([3+(i-1)*3 3+i*3],:)=Fw([3+(i-1)*3 3+i*3],:)-w*ele_he(i)/2;
end    
clear i
 
deltaA      =   zeros(node_num*3,1); 
deltaU      =   zeros(node_num*3,1);
deltaU1     =   zeros(node_num*3,1);
deltaR      =   zeros(node_num*3,1);


Uorig           =   tU;
Vorig           =   tV;
Aorig           =   tA;
tU(float_pt_node(1))=xL-L;%float_point(1)-float_pt_pos(1);
tU(float_pt_node(3))=zL;%float_point(3)-float_pt_pos(3);

tR  =   [];
tR  =   Fw;
tR(fix_pt_node(1))=tR(fix_pt_node(1))-1*T0;
tR(float_pt_node(1))=tR(float_pt_node(1))+1*T0;
T0=0;

first       =   1;
convergence =   0;

[ele_cos u] =   updateInfo(x,y,z,tU,ele_num,ele_info,order);
dir =   direction(deltaR,tR,node_num,fix_pt_node,float_pt_node,ele_cos);

RNORM       =   norm(tR);
visD        =   zeros(node_num*3,node_num*3);
nu0     =   nu;
N0=zeros(ele_num,1);
while convergence  ==  0
    % Nonlinear calculation 
    deltaR = tR     -   tF; 
    [K_L  Ma]  =   KL(Ca,k,EA,order,ele_num,ele_info,ele_he,ele_cos,rho_water,Area);
    C_nu    =    EA*2/min(ele_he)*k(1);
    K_NL            =   KNL(tN,k,ele_he,ele_num,ele_info,order);
    visD       =   nu0*C_nu*eye(node_num*3)/deltaT;
    K           =   K_L     +    K_NL   +visD;
    if first    ==  1
        deltaU1 =   displacement(K,node_num,fix_pt_node,deltaR,float_pt_node,dir);
        deltaU  =   deltaU1;
        R0      =   deltaR;
        first   =   0;
    else
        deltaU  =   displacement(K,node_num,fix_pt_node,deltaR,float_pt_node,dir);
    end
    tU          =   tU  +   deltaU;
    convergence = E_conv(R0,deltaR,deltaU,deltaU1,RNORM,ETOL,RTOL,DTOL,fix_pt_node,float_pt_node);
    nu0      =   nu0*(sum(deltaU.^2)/sum(tU.^2))^0.5;
    deltaU = tU-Uorig;
    [a_tmp  v]  =   solNewmark(a0,a2,a3,a6,a7,Aorig,Vorig,deltaU);
    tA=a_tmp;
    tV=v;
    [ele_cos u] =   updateInfo(x,y,z,tU,ele_num,ele_info,order);
	[tN tF]     = normalForce(k,EA,ele_X,u,ele_num,ele_he,ele_cos,ele_info,order,node_num,N0);
end

 clear a_tmp v nu0 
 
  X_tmp     =   zeros(node_num,3);
  if float_point(1)-fixed_point(1)>0
        X_tmp(:,1)=   (x+tU(1:3:node_num*3))*Cos(1)+fixed_point(1);
        X_tmp(:,2)=   (y+tU(2:3:node_num*3))*Cos(2)+fixed_point(2);
        X_tmp(:,3)=   z+tU(3:3:node_num*3)+fixed_point(3);
  else
       X_tmp(:,1)=   ((x+tU(1:3:node_num*3))*-1+xL)*Cos(1)+float_point(1);
       X_tmp(:,2)=   (y+tU(2:3:node_num*3))*Cos(2)+fixed_point(2);
       X_tmp(:,3)=   z+tU(3:3:node_num*3)+fixed_point(3);
  end

float_pt_pos    =   float_point;
fix_pt_pos      =   fixed_point;  
  
plot(X_tmp(:,1),X_tmp(:,3),'o-')
  