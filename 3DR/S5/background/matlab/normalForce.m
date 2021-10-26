function [tN tF]= normalForce(k,EA,ele_X,u,ele_num,ele_he,ele_cos,ele_info,order,node_num,N0)
    tN   =   zeros(ele_num,1);
    tF  =   zeros(node_num*3,1);
    
    for i   =   1:ele_num
       kij  =   [];
       kij  =   EA*2/ele_he(i)*k;
       f    =   kij*(u(i,:)-ele_X(i,:))'+[-N0(i) N0(i)]';
       tN(i) =   abs(f(1));
       
       
       T   =   zeros(order+1,(order+1)*3);
       R   =   [ele_cos(i,1)    ele_cos(i,2)   ele_cos(i,3);];
       for j    =   1:3
            loc         =   [1:order+1;j:3:(order+1)*3;];
            loc(2,:)    =   (loc(2,:)-1)*(order+1);
            loc         =   sum(loc);
            T(loc)      =   R(j);
       end
        node=[];
        for j   =   1:order+1
           node=[node [1 2 3]+ 3*(ele_info(i,j)-1)];
        end
       tF(node,1)    =   tF(node,1)   +   T'*f;
    end
end