
deltaA      =   zeros(node_num*3,1); 
deltaU      =   zeros(node_num*3,1);
deltaU1     =   zeros(node_num*3,1);
deltaR      =   zeros(node_num*3,1);
first       =   1;
convergence =   0;

[ele_cos u] =   updateInfo(x,y,z,tU,ele_num,ele_info,order);
[tN tF]     = normalForce(k,EA,ele_X,u,ele_num,ele_he,ele_cos,ele_info,order,node_num,N0);    

deltaR      =   tR  -   tF  -C*tV   -   Ma*tA;

dir         =   direction(deltaR,tR,node_num,fix_pt_node,float_pt_node,ele_cos);

%%
tA          =   initial_A(M0,node_num,fix_pt_node,tR,tF,float_pt_node);
                
RNORM       =   norm(tR);

while convergence  ==  0
      
    % Update information of internal force
    [ele_cos u] =   updateInfo(x,y,z,tU,ele_num,ele_info,order);
    [tN tF]     =   normalForce(k,EA,ele_X,u,ele_num,ele_he,ele_cos,ele_info,order,node_num,N0); 
    
    % Nonlinear calculation 
     deltaR = tR       +   (M0+Ma)*(a0*Uorig    +a2*Vorig     +a3*Aorig)  +   C*(a1*Uorig    +a4*Vorig  +a5*Aorig)-a0*(M0+Ma)*tU-a1*C*tU-tF; 
      
    [K_L Ma]    =   KL(Ca,k,EA,order,ele_num,ele_info,ele_he,ele_cos,rho_water,Area);
    K_NL        =   KNL(tN,k,ele_he,ele_num,ele_info,order);
    Keff        =   a0*(M0+Ma)     +   a1*C;
    
    K           =   Keff    +K_L     +    K_NL   ;
    
    if first    ==  1
        deltaU1 =   displacement(K,node_num,fix_pt_node,deltaR,float_pt_node,dir);
        deltaU  =   deltaU1;
        R0      =   deltaR;
        first   =   0;
    else
        deltaU  =   displacement(K,node_num,fix_pt_node,deltaR,float_pt_node,dir);
    end
    
    % check Convergence
    convergence =   E_conv(R0,deltaR,deltaU,deltaU1,RNORM,ETOL,RTOL,DTOL,fix_pt_node,float_pt_node)
    tU          =   tU  +   deltaU;
    deltaU  =   tU-Uorig;
    [a  v]  =   solNewmark(a0,a2,a3,a6,a7,Aorig,Vorig,deltaU);
    tA=a;
    tV=v;
    
end
clear a v 