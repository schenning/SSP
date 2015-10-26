% HW , Henning Schei
%clear all;
%close all; 

%Collecting data from 3 different places


fileID = fopen('hw_host_ip_list.txt','r');
fSpec  = '%s';
data   = textscan(fileID, fSpec);
fclose(fileID);
res=zeros(100,3);
res(1:100,1) = pingstats('mercury.iet.ntnu.no' , 100, ''); 
res(1:100,2) = pingstats('venus.iet.ntnu.no', 100,'');
%res(1:100,3) = pingstats('folk.ntnu.no/henninsc',100,'');
res(1:100,3) = pingstats('archlinux.uib.no',100,'');

%%

% Plotting histograms
hist(res(:,1), 50);
xlabel 'msec'
ylabel '#ping'
title 'mercury.iet.ntnu.no'
figure;
hist(res(:,2),50);
xlabel 'msec'
ylabel '#ping'
title 'venus.iet.ntnu.no'
figure;
hist(res(:,3),50);
xlabel 'msec'
ylabel '#ping'
title 'archlinux.uib.no'

%% Problem i) Maximum Likelihood estimation 


dis = { 'Normal', 'Rayleigh','Exponential'}; % The Erlang dist was not a part of the possibilities

phat = zeros(100,3);
for i=1:3
    phat(:,i) = mle(res(:,i),'distribution','Rayleigh');
    
end

mu_n=0; ro_n=0; %mu_normaldist and ro_normaldist

% Need separate for-loop bc of two parameters in normal dist

for i=1:3
    [mu_n, ro_n] = mle(res(:,i), 'distribution', 'Normal');
end



















%res = res';

