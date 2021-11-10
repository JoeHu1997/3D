function bF = seabedForce(bF,tU,ele_num,w,L,Fw,f)
    EA = -Fw(6)/f
    %N = 1;
    %B0 = fixed_point(3)- N *D;
    

for no_ele=2:ele_num
    if tU((no_ele*3)) < f
        bF((no_ele*3),1)= (bF((no_ele*3),1)+(tU(no_ele*3) * EA ))/2 ;
    elseif tU((no_ele*3))  >f && tU((no_ele*3))  < 0
        bF((no_ele*3),1)= (bF((no_ele*3),1)+tU(no_ele*3) * EA)/2  ;
    else
        bF((no_ele*3),1) = 0;
        
    end
end



