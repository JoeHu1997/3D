function    dir =   direction(deltaR,tR,node_num,fix_pt_node,float_pt_node,ele_cos)
    dir         =   [0 0 0];
    node        =   zeros(node_num*3,1);
    node1       =   zeros(node_num*3,1);
    node2       =   zeros(node_num*3,1);
    tmp         =   find(deltaR==0);
    node1(tmp)  =   1;
    tmp         =   find(tR==0);
    node2(tmp)  =   1;
    for i=1:node_num*3
        if node1(i)==1 && node2(i)==1
           node(i)=1; 
        end
    end
    node([fix_pt_node float_pt_node])=[];
    
    if sum(node(2:3:end))==length(node)/3
        dir(2)   =   1; 
    end
    if sum(node(3:3:end))==length(node)/3
        dir(3)   =   1; 
    end
    
    for i=1:3
       tmp  =   abs(sum(ele_cos(:,i)));
       if tmp > 0
           dir(i)=0;
       end
    end
    
end