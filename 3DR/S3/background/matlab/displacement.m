function u  =   displacement(K,node_num,fix_pt_node,R,float_pt_node,dir)
%% displacement matrix
    tmp     =  R==0;
    if sum(tmp)==node_num*3
        u       =   zeros(node_num*3,1);
        
    else
        u       =   zeros(node_num*3,1);
        u_fix   =   u;
        K_fix   =   K;
        R_fix   =   R;
        node    =   1:node_num*3;
        fix_pt_node    =   [fix_pt_node; float_pt_node;];
    
        %   unknown force row Kij-->0
        K_fix(fix_pt_node,:)=[]; %unknown force row Kij-->0
        K_fix(:,fix_pt_node)=[]; %zero displacement column Kij-->0
        node(fix_pt_node)=[];
        u_fix(fix_pt_node)=[];
        R_fix(fix_pt_node)=[];
        
        tmp     =   [];
        for i=1:3
            if dir(i)==1
                tmp=[tmp i:3:length(node)];
            end
        end
        K_fix(tmp,:)=[]; %unknown force row Kij-->0
        K_fix(:,tmp)=[]; %zero displacement column Kij-->0
        node(tmp)=[];
        u_fix(tmp)=[];
        R_fix(tmp)=[];
                
        
        if length(node)==1
            u_fix=R_fix/K_fix;
            u(node)=u_fix;
        else
            u_fix=K_fix\R_fix;
            u(node)=u_fix;
        end
    end
    
end
 


