 X_tmp  =   zeros(node_num,3);
 zL     =   float_point(3)-fixed_point(3);
 xL     =   ((float_point(1:2)-fixed_point(1:2))*(float_point(1:2)-fixed_point(1:2))')^0.5;
 Cos    =    (float_point(1:2)-fixed_point(1:2))/xL;
 
 if xL  ==  0
     Cos    =   [1 0];
 end
 
 X_tmp(:,3)=linspace(0,zL,node_num);
 X_tmp(:,2)=Cos(2)*linspace(0,xL,node_num);
 X_tmp(:,1)=Cos(1)*linspace(0,xL,node_num);
 
 
 % initialize 
 x  =   X_tmp(:,1);
 y  =   X_tmp(:,2);
 z  =   X_tmp(:,3);

 
tN(:)=T0;

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

fix_pt_pos   = fixed_point;
float_pt_pos = float_point;
fix_pt_dir   =   [1   1   1];
float_pt_dir =   [1   1   1];
fix_pt_node=[];
float_pt_node=[];

dir1=find(fix_pt_dir(:));
loc1=[find(abs(x-fix_pt_pos(1)) < 0.001); find(abs(y-fix_pt_pos(2)) < 0.001); find(abs(z-fix_pt_pos(3)) < 0.001);];
fix_pt_node=[fix_pt_node 3*(mode(loc1)-1)+dir1;];
dir2=find(float_pt_dir(:));
loc2=[find(abs(x-float_pt_pos(1)) < 0.001); find(abs(y-float_pt_pos(2)) < 0.001); find(abs(z-float_pt_pos(3)) < 0.001);];
float_pt_node=[float_pt_node 3*(mode(loc2)-1)+dir2;];


shape
k  =    local_k(shape_fuc,order); 
[a0    a1  a2  a3  a4  a5  a6  a7] =   paraNewmark(a,b,deltaT);
for i=1:ele_num
        Fw([3+(i-1)*3 3+i*3],:)=Fw([3+(i-1)*3 3+i*3],:)-w*ele_he(i)/2;
end    
clear i dir1 dir2 loc1 loc2