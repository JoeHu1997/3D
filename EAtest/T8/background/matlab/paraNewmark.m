function [a0    a1  a2  a3  a4  a5  a6  a7] =   paraNewmark(a,b,deltaT) 
    a0 =    1/(a*deltaT*deltaT);
    a1 =    b/(a*deltaT);
    a2 =    1/(a*deltaT);
    a3 =    1/(2*a)-1;
    a4 =    b/a-1;
    a5 =    deltaT/2*(b/a-2);
    a6 =    deltaT*(1-b);
    a7 =    b*deltaT;
end