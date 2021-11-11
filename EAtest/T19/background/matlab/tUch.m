function tUcheck = tUch(deltaU,node_num,f)    %are there any point laying
tUcheck  = 1;
   for i = 1:node_num-2
       if deltaU( 3*(i))< f
         
           tUcheck  = 0;
           
       end
       
   end
   
end