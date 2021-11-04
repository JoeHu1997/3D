function [ele_cos u] =   updateInfo(x,y,z,tU,ele_num,ele_info,order)

    L       =   length(tU);
    X       =   x+tU(1:3:L);
    Y       =   y+tU(2:3:L);
    Z       =   z+tU(3:3:L);
    ele_cos =   zeros(ele_num,3);
    u       =   zeros(ele_num,order+1);
    
    
    for i   =   1:ele_num
        % obtain new radian of element in global cooridinate        
        L       =   ((X(ele_info(i,1))-X(ele_info(i,order+1)))^2+(Y(ele_info(i,1))-Y(ele_info(i,order+1)))^2+(Z(ele_info(i,1))-Z(ele_info(i,order+1)))^2)^0.5;
        deltaX  =   X(ele_info(i,order+1))-X(ele_info(i,1));
        deltaY  =   Y(ele_info(i,order+1))-Y(ele_info(i,1));
        deltaZ  =   Z(ele_info(i,order+1))-Z(ele_info(i,1));
        ele_cos(i,1)=deltaX/L;
        ele_cos(i,2)=deltaY/L;
        ele_cos(i,3)=deltaZ/L;
        
        % obtain relative u of each element ( assume X1=0 )        
        Pos             =   zeros(order+1,3);
        Pos(:,1)        =   X(ele_info(i,1:order+1));
        Pos(:,2)        =   Y(ele_info(i,1:order+1));
        Pos(:,3)        =   Z(ele_info(i,1:order+1));
        Pos             =   Pos-Pos(1,:);
        Pos             =   Pos'.^2;
        u(i,:)          =   sum(Pos).^0.5;
    end

end