
% HW , Henning Schei
%clear all;
%close all; 

%Collecting data from 3 different places

choise = 0;

if (choise)
    fileID = fopen('hw_host_ip_list.txt','r');
    fSpec  = '%s';
    data   = textscan(fileID, fSpec);
    fclose(fileID);
    res=zeros(1000,4);
    res(1:1000,1) = pingstats('mercury.iet.ntnu.no' , 1000, ''); % Norwegian Univ. of Science and Technology #1 in res
    res(1:1000,2) = pingstats('atalante.stanford.edu', 1000,''); % Standford University                      #2 in res
    res(1:1000,3) = pingstats('mx.vvsu.ru',1000,''); % Vladivostok State University ofEconimics and Service  #3 in res
    res(1:1000,4) = pingstats('197.255.176.1',1000,''); %Brazzaville                                         #didn't work

    %res(1:100,3) = pingstats('archlinux.uib.no',100,'');
    % 217.74.123.50 % ISP = Vladivostok SU
    % 197.214.128.4 Airtel Congo. % sykt treg
    % 197.255.176.1 Brazzaville   
    % 197.220.64.1 , Mogadishu    % relativ rask
else
    % unpack from .mat file pingmania.mat
    tmp  = load('pingmania.mat');
    res  = tmp.res(:,:); % fix 4.colunm 
end



%%

% Plotting histograms
hist(res(:,1), 500);
xlabel 'msec'
ylabel '#ping'
title 'mercury.iet.ntnu.no - Located in Trondheim, Norway';
figure;
hist(res(:,2),500);
xlabel 'msec'
ylabel '#ping'
title 'atlante.standford.edu - Located in California';
figure;
hist(res(:,3),500);
xlabel 'msec'
ylabel '#ping'
title 'mx.vvsu.ru - Located in Vladivostok, Russia.';



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

