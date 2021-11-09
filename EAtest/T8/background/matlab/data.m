%input  point position,order,number of element
function [x y z order node_num ele_info ele_he ele_num ele_X]=data(no)

%data form x y z 
%% poistion
    load(['ele'  num2str(no) '.mat'])
    node_num    =   length(x);
    
%% information of element   
    load(['node'  num2str(no) '.mat'])
    ele_num     =   length(ele_info(:,1)); % number of element
    order       =   length(ele_info(1,:))-2; % order
    %ele_he      =   zeros(ele_num,1);
    ele_X       =   zeros(ele_num,order+1);
    for     i   =   1:ele_num
        % obtain new radian of element in global cooridinate 
     %   ele_he(i)   =   ((x(ele_info(i,1))-x(ele_info(i,order+1)))^2+(y(ele_info(i,1))-y(ele_info(i,order+1)))^2+(z(ele_info(i,1))-z(ele_info(i,order+1)))^2)^0.5;
        % obtain relative u of each element ( assume X1=0 )        
        Pos             =   zeros(order+1,3);
        Pos(:,1)        =   x(ele_info(i,1:order+1));
        Pos(:,2)        =   y(ele_info(i,1:order+1));
        Pos(:,3)        =   z(ele_info(i,1:order+1));
        Pos             =   Pos-Pos(1,:);
        Pos             =   Pos'.^2;
        ele_X(i,:)      =   sum(Pos).^0.5;
    end
end