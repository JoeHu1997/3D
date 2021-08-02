function K_NL =   KNL(N,k,ele_he,ele_num,ele_info,order)
    
    K_NL=zeros(3*(order+ele_num),3*(order+ele_num));
    
    for no_ele=1:ele_num    %no of element : no_ele
        kij     =   N(no_ele)*2/ele_he(no_ele)*k;
        
        node=[];
        for i=1:order+1
           node=[node [1; 2; 3;]+ 3*(ele_info(no_ele,i)-1)];
        end       
        
        K_NL(node(1,:),node(1,:))=K_NL(node(1,:),node(1,:))+kij(:,:);
        K_NL(node(2,:),node(2,:))=K_NL(node(2,:),node(2,:))+kij(:,:);
        K_NL(node(3,:),node(3,:))=K_NL(node(3,:),node(3,:))+kij(:,:);
    end

end