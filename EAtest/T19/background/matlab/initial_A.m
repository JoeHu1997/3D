function a  =   initial_A(M,node_num,fix_pt_node,tR,tF,float_pt_node)
%% displacement matrix
   if   sum(M)==0
       a       =   zeros(node_num*3,1);
   else
    a       =   zeros(node_num*3,1);
    a_fix   =   a;
    M_fix   =   M;
    R_fix   =   tR-tF;
    node    =   1:node_num*3;
    fix_pt_node    =   [fix_pt_node; float_pt_node;];
    
%   unknown force row Kij-->0
    M_fix(fix_pt_node,:)=[]; %unknown force row Kij-->0
    M_fix(:,fix_pt_node)=[]; %zero displacement column Kij-->0
    node(fix_pt_node)=[];
    a_fix(fix_pt_node)=[];
    R_fix(fix_pt_node)=[];
   
    % row Mij-->0 for 0
    tmp= sum(M_fix)==0;
    for i=1:length(tmp)
        if tmp(length(tmp)+1-i)==1
            M_fix(length(tmp)+1-i,:)=[];
            M_fix(:,length(tmp)+1-i)=[];
            node(length(tmp)+1-i)=[];
            a_fix(length(tmp)+1-i)=[];
            R_fix(length(tmp)+1-i)=[];
        end
    end
    a_fix=M_fix\R_fix;
    a(node)=a_fix;
   end
end
 
