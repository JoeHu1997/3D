%% input BCs
FileBCs         =   ['point'    num2str(cable_no) '.txt'];
pt  =   dlmread(FileBCs);
clear FileBCs

% fixed point node and position
fix_pt_pos  =   pt(1,:);
fix_pt_dir  =   [1   1   1];
fix_pt_node=[];
dir=find(fix_pt_dir(:));
loc=[find(abs(x-fix_pt_pos(1)) < 0.001); find(abs(y-fix_pt_pos(2)) < 0.001); find(abs(z-fix_pt_pos(3)) < 0.001);];
fix_pt_node=[fix_pt_node 3*(mode(loc)-1)+dir;];
clear loc dir

% floating point and position
float_pt_pos  =   pt(2,:);
float_pt_dir  =   [1   1   1];
float_pt_node=[];
dir=find(float_pt_dir(:));
loc=[find(abs(x-float_pt_pos(1)) < 0.001); find(abs(y-float_pt_pos(2)) < 0.001); find(abs(z-float_pt_pos(3)) < 0.001);];
float_pt_node=[float_pt_node 3*(mode(loc)-1)+dir;];
clear loc dir


if pt(3,1)==1
    float_pt_node=[];
end