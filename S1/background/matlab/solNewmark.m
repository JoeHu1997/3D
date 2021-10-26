function    [A  V]  =   solNewmark(a0,a2,a3,a6,a7,tA,tV,deltaU)
    A = a0*deltaU    -a2*tV  -a3*tA;
    V = tV           +a6*tA  +a7*A;
end