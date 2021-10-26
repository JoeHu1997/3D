function [convergence var4]= E_conv(R0,deltaR,deltaU,deltaU1,DNORM,RNORM,ETOL,RTOL,DTOL,fix_pt_node,float_pt_node,timer1,tR,tU)    %check the convergnece of energy  
    check = [0 0 0];
    if RNORM==0
         RNORM=1;
     end
    
    fix_pt_node    =   [fix_pt_node;  float_pt_node;];
    deltaU1(fix_pt_node)    =   [];
    deltaU(fix_pt_node)     =   [];
    deltaR(fix_pt_node)     =   [];
    R0(fix_pt_node)     =   [];
    tR(fix_pt_node)     =   [];
    tU(fix_pt_node)     =   [];
    
    %Convergence critria for energy
    var1    =  deltaU1'*(R0);
    var2    =  deltaU'*(deltaR);
    r       =   var2/var1;
    
    if  sum(deltaU)< 10^-8 || r   <   ETOL
        check(1)    =   1;
    end
    
    %   Convergence critria for force
    RNORM = norm(tR);
    %RNORM = 0.01;
    var3     =   norm(deltaR);
    r       =   var3/RNORM;
    if  r   <   RTOL
        check(2)    =   1;
    end
    
    %   Convergence critria for displacement
    DNORM = norm(max(tU));
    var4     =   norm(deltaU);  
    r       =   var4/DNORM;
    if  r     <   DTOL
        check(3)    =   1;
    end
    %check
    
    %   Fulfill all convergence
    if timer1 >500
        check = [1 1 1];
    end 
    if  sum(check)==3
        convergence     =   true;
    else 
        convergence     =   false;
    end
    if sum(R0)==0
        convergence     =   true;
    end
    
end