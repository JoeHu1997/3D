function [K_L  Ma]=   KL(Ca,k,EA,order,ele_num,ele_info,ele_he,ele_cos,rho_water,Area)

    K_L     =   zeros(3*(order+ele_num),3*(order+ele_num));
    Ma      =   zeros(3*(order+ele_num),3*(order+ele_num)); 
    
    
    for no_ele=1:ele_num    %no of element : no_ele
        kij     =   EA*2/ele_he(no_ele)*k;
        ma      =   Ca*rho_water*Area*ele_he(no_ele)/2;
        
%% Assemble local coordinate into global coordinate
        T   =   zeros(order+1,(order+1)*3);
        R   =   [ele_cos(no_ele,1)    ele_cos(no_ele,2)   ele_cos(no_ele,3);];
        for     j   =   1:3
            loc         =   [1:order+1;j:3:(order+1)*3;];
            loc(2,:)    =   (loc(2,:)-1)*(order+1);
            loc         =   sum(loc);
            T(loc)      =   R(j);
            clear       loc
        end
        
        kij=T'*kij*T;
        ma = ma*T'*T;
        % node : the no of node assemble to global coordinate
        node=[];
        for i=1:order+1
           node=[node [1 2 3]+ 3*(ele_info(no_ele,i)-1)];
        end
        
        K_L(node,node)=K_L(node,node)+kij(:,:);
        Ma(node,node)=Ma(node,node)+ma(:,:);
    end
    Ma  =   -1*Ma;
end