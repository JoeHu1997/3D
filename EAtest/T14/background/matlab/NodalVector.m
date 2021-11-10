%%  NodalVector
%   Calculate the vector of each segment
%   rZ:     vector along the cable    
function    [vT vN l rZ]  =   NodalVector(x,y,z,tU,node_num,tV,field_V)  
    
    vR      =   tV  -   field_V;    %relative velocity of each node
    rZ      =   zeros(node_num-1,3);    
    vr      =   zeros(node_num-1,3);
    l       =   zeros(node_num-1,1);
    X       =   x+tU(1:3:end);
    Y       =   y+tU(2:3:end);
    Z       =   z+tU(3:3:end);
    vT      =   zeros(node_num-1,3);
    vN      =   zeros(node_num-1,3);
    
    for     no  =   1:node_num-1
        tmp         =  [(X(no+1)-X(no))  (Y(no+1)-Y(no))  (Z(no+1)-Z(no))];
        l(no)       =  norm(tmp);
        rZ(no,:)    =  tmp/l(no);   %tangent vector
        vr(no,:)    =  (vR(1+(no)*3:3+(no)*3)+vR(1+(no-1)*3:3+(no-1)*3))/2;
        vT(no,:)    =  (rZ(no,:)*vr(no,:)')*rZ(no,:);
        vN(no,:)    =   vr(no,:)-vT(no,:);       
    end


end