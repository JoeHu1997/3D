cd ([root '/system'])
fid     =   ['pointLocation' num2str(cable_no) '.txt'];
fileID  =   fopen(fid,'w+');
fprintf(fileID,'( %.12f %.12f %.12f )\n',X_tmp');

fclose(fileID);


