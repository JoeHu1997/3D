
call = call+1;	
cd ([root '/postProcessing/line' num2str(cable_no) '/0' ])

fid  =   fopen('U','r');
dataU =   textscan(fid,'%s %f %s',node_num,'headerline',call);


if First
	First	=	false;
	field_V	=	zeros(node_num*3,1);
else
	field_V(2:3:node_num*3,1)	=	dataU{2};
	for	i=1:node_num
		% fluid velocity
		tmp	=	dataU{1}{i};
		tmp(1)	=	[];
		tmp	=	str2num(tmp);
		if 	abs(tmp) > 500
			tmp=0;
            field_V(3*i-1,1)	=	tmp;
		end		
		field_V(3*i-2,1)	=	tmp;	
		tmp	=	dataU{3}{i};
		tmp(length(tmp))	=	[];
		tmp	=	str2num(tmp);
		if 	abs(tmp) > 500
			tmp=0;
		end	
		field_V(3*i,1)	=	tmp;
    end
end		
fclose(fid);
clear dataP dataU dataA fid tmp
cd (root) 

if FOAMcall    <   MATLABcall    		
	FOAMcall    =   FOAMcall+1;
    	call = call+1;
end

cd (root)
