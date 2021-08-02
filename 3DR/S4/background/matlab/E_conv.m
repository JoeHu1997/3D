function convergence = E_conv(R0,deltaR,deltaU,deltaU1,RNORM,ETOL,RTOL,DTOL,fix_pt_node,float_pt_node)    %check the convergnece of energy  
    check = [0 0 0];
    if RNORM==0
        RNORM=1;
    end
    fix_pt_node    =   [fix_pt_node; float_pt_node;];
    deltaU1(fix_pt_node)    =   [];
    deltaU(fix_pt_node)     =   [];
    deltaR(fix_pt_node)     =   [];
    R0(fix_pt_node)     =   [];
    
    %   Convergence critria for energy
    var1    =  deltaU1'*(R0);
    var2    =  deltaU'*(deltaR);
    r       =   var2/var1;
    
    if  sum(deltaU)< 10^-8 || r   <   ETOL
        check(1)    =   1;
    end
    
    %   Convergence critria for force
    var     =   norm(deltaR);
    r       =   var/RNORM;
    
    if  r   <   RTOL
        check(2)    =   1;
    end
    
    %   Convergence critria for displacement
    var     =   norm(deltaU);
    
    if  var     <   DTOL
        check(3)    =   1;
    end

    
    %   Fulfill all convergence
    if  sum(check)==3
        convergence     =   true;
    else 
        convergence     =   false;
    end
    if sum(R0)==0
        convergence     =   true;
    end
    
end