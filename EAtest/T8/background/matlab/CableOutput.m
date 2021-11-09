load(['Data'  num2str(cable_no ) '.mat'])
node_num       =   1   +   ele_num*order;

% node of each element
ele_info        =   zeros(ele_num,(order+1)+1);
I               =   ones(1,ele_num);
ele_info        =   I'.*[1:order+1 1];   
I               =   ones(1,(order+1)+1);
ele_info        =   ele_info    +   order*I.*((1:ele_num)-1)';    

fid=['node' num2str(cable_no) '.mat'];
save(fid,'ele_info')
fclose('all');

% position of each node
Pos     =   zeros(node_num,3);
x=[];
y=[];
z=[];
x   =   X_tmp(:,1);
y   =   X_tmp(:,2);
z   =   X_tmp(:,3);
fid=['ele' num2str(cable_no) '.mat'];
%save(fid,'Pos');
save(fid,'x','y','z','ele_he');
fclose('all');

fix_pt_pos  =   Pos(fix_pt_node(3)/3,:);
float_pt_pos=   Pos(float_pt_node(3)/3,:);

% file for sampling
fid     =   ['pointLocation' num2str(cable_no) '.txt'];
fileID  =   fopen(fid,'w');
for     i   =   1:node_num
    fprintf(fileID,'( %.3f %.3f %.3f )\n',x(i),y(i),z(i));
end

% file for BCs
fid     =   ['point' num2str(cable_no) '.txt'];
fileID  =   fopen(fid,'w');
fprintf(fileID,'%.3f %.3f %.3f \n',fixed_point);
fprintf(fileID,'%.3f %.3f %.3f \n',float_point);
fprintf(fileID,'%d \n',free);
fclose('all');

tF-tR;
Adjust=ans(float_pt_node(3));
N0=tN;

% file for parameter
fid=['moordata' num2str(cable_no) '.mat'];
save(fid,'EA','w','T0','g','ele_num','lumpedM','C','tU','tV','tA','tN','tF','tR','Fex','N0','a','b','ETOL','RTOL','DTOL','deltaT','Ca','Cc','Cdt','Cdn','D','nu','Area','rho_water','rho_mass','free','X_tmp','Adjust')
fclose('all');
cd ..