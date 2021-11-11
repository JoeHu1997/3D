function k  = local_k(shape_fuc,order)

%% Build local stiffness
    %kNL=2/ele_he(no_ele)*[shape_fuc(order).dphi'*shape_fuc(order).weight*shape_fuc(order).dphi];
    %kNL=2/ele_he(no_ele)*kg;
    %kL =EA*2/ele_he(no_ele)*[shape_fuc(order).dphi'*shape_fuc(order).weight*shape_fuc(order).dphi];
    %kL =EA*2/ele_he(no_ele)*kL
%    kij=    zeros(3*(order+1),3*(order+1));
    k  =   [];
    k  =   [shape_fuc(order).dphi'*shape_fuc(order).weight*shape_fuc(order).dphi];
    
%    for i   =   1:order+1
%       for j    =   1:order+1
%           kij(1+3*(i-1):3+3*(i-1),1+3*(j-1):3+3*(j-1))     =   kij(1+3*(i-1):3+3*(i-1),1+3*(j-1):3+3*(j-1))+k(i,j)*eye(3);
%           
%       end
%    end
end

% order大於2的shape_func須再做檢查