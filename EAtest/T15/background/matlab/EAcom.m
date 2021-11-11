clear all
close all
%%
timer = 0;
EA = linspace(500,200000,400);
for ii = 1:1:400

L               =   1.95;
w               =   1.243;
node_num        =   11;
order           =   1;                 
ele_num         =   10;
node_num       =   1   +   ele_num*order;
ele_info        =   zeros(ele_num,(order+1)+1);
I               =   ones(1,ele_num);
ele_info        =   I'.*[1:order+1 1];   
I               =   ones(1,(order+1)+1);
ele_info        =   ele_info    +   order*I.*((1:ele_num)-1)';
 
C   =   zeros(node_num*3,node_num*3);
tU  =   zeros(node_num*3,1); 
tV  =   zeros(node_num*3,1); 
tA  =   zeros(node_num*3,1);   
tF  =   zeros(node_num*3,1);
tR  =   zeros(node_num*3,1);
tN  =   zeros(node_num,1);
Fw =    zeros(node_num*3,1); 

fixed_point     =   [1.917 0 -0.9];  % Position of Fixed Point (anchor)
float_point     =   [0.2575 0 -0.076];  % Position of Connecting Point (refAttachmentPt)
xL          =   ((float_point(1:2)-fixed_point(1:2))*(float_point(1:2)-fixed_point(1:2))')^0.5;
zL          =   float_point(3)-fixed_point(3);
Cos         =   (float_point(1:2)-fixed_point(1:2))/xL;
Cos         =   abs(Cos);

fix_pt_pos      =   [0.00  0.0     0.00];
float_pt_pos    =   [L*Cos(1)     L*Cos(2)     0.00];
%%
%{
Va = linspace(-100,50,44);
EA =  linspace(100,2000,20);
A = zeros(44,1);
B = zeros(44,1);
figure(1)
hold off
for ii = 15:20
    for i = 1:44
    syms    H;
    eqn =   xL   ==  H/EA(ii)*L+H/w*(asinh((w*L+Va(i))/H)-asinh(Va(i)/H));
    S   =   vpasolve(eqn,H);
    A(i)  =   double(S);
    
    end
   
    
    
    for i = 1:44
    syms    H;
    eqn2 =   zL   == 1/EA(ii)*(0.5*w*L^2+Va(i)*L)+H/w*((1+((w*L+Va(i))/H)^2)^0.5-(1+(Va(i)/H)^2)^0.5);
    S1   =   vpasolve(eqn2,H);
    B(i)  =   double(S1);
    end
    plot(Va,A,Va,B)
    axis([-200,100,-200,50])
    pause(0.1)
end
 %}
%%
syms    H  Va
f1  =   xL   ==  H/EA(ii)*L+H/w*(asinh((w*L+Va)/H)-asinh(Va/H));
f2  =   zL   ==  1/EA(ii)*(0.5*w*L^2+Va*L)+H/w*((1+((w*L+Va)/H)^2)^0.5-(1+(Va/H)^2)^0.5);

S   =   vpasolve([f1 f2]);

H   =   double(S.H);
Va  =   double(S.Va);
T0(ii)  =   abs(H);
clear f1 f2  S
%%
x                =  linspace(0,float_pt_pos(1),node_num)';
y                =  linspace(0,float_pt_pos(2),node_num)';
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
tU(float_pt_node(1))=xL*Cos(1)-L*Cos(1);%float_point(1)-float_pt_pos(1);
tU(float_pt_node(2))=xL*Cos(2)-L*Cos(2);
tU(float_pt_node(3))=zL;%float_point(3)-float_pt_pos(3);

tR  =   [];
tR  =   Fw;
tR(fix_pt_node(1))=tR(fix_pt_node(1))-1*T0(ii);
tR(float_pt_node(1))=tR(float_pt_node(1))+1*T0(ii);

timer = timer + 1;
RNORM(ii)       =   norm(tR);
figure(1)

end
figure(1)
plot(EA,T0,EA,RNORM)



