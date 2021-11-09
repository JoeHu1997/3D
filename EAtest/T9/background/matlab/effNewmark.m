function    [ Keff    Reff ]   =    effNewmark(a0,a1,a2,a3,a4,a5,M,C,tR,tU,tV,tA)   
    
    Keff = a0*M     +   a1*C;
    Reff = tR       +   M*(a0*tU    +a2*tV     +a3*tA)  +   C*(a1*tU    +a4*tV  +a5*tA); 

end