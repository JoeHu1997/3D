function M = mass(rho_mass,rho_water,Area,Ca,order,ele_num,ele_info,ele_he)
    M   =   zeros(3*(order+ele_num),3*(order+ele_num));
    rho =   rho_mass    +   Ca*rho_water;
    for no_ele=1:ele_num    %no of element : no_ele
        mij     =   rho*Area*ele_he(no_ele)/2*eye(2);
        
%% Assemble local coordinate into global coordinate
        % node : the no of node assemble to global coordinate
        node=[];
        for i=1:order+1
           node=[node [1; 2; 3;]+ 3*(ele_info(no_ele,i)-1)];
        end  
        
        M(node(1,:),node(1,:))=M(node(1,:),node(1,:))+mij(:,:);
        M(node(2,:),node(2,:))=M(node(2,:),node(2,:))+mij(:,:);
        M(node(3,:),node(3,:))=M(node(3,:),node(3,:))+mij(:,:);
    end
end