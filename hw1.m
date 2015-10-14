% HW , Henning Schei
%clear all;
%close all; 
fileID = fopen('hw_host_ip_list.txt','r');
fSpec  = '%s';
data   = textscan(fileID, fSpec);
fclose(fileID);
[res] = pingstats(data{1}{1}, 2, ''); 
res = res';


