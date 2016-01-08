
% HW, Statistical Signal Processing,
% Written by Henning Schei 

%Collecting data from 3 different places
choise = 0;
close all;
if (choise)
    res=zeros(1000,4);
    res(1:1000,1) = pingstats('mercury.iet.ntnu.no' , 1000, ''); % Norwegian Univ. of Science and Technology #1 in res
    res(1:1000,2) = pingstats('atalante.stanford.edu', 1000,''); % Standford University                      #2 in res
    res(1:1000,3) = pingstats('mx.vvsu.ru',1000,'');% Vladivostok State University of Econimics and Service  #3 in res
    res(1:1000,4) = pingstats('197.255.176.1',1000,''); %Brazzaville, Congo                                           
else
    
    % Check if data file exists
    if exist('pingmania.mat','file') == 2
        % unpack from .mat file pingmania.mat
        tmp  = load('pingmania.mat');
        res  = tmp.res(:,:); % fix 4.colunm
    else
        fprintf('You need to add til file pingmania.mat to the current directory\n');
        return
    end
    
end


data = res(1:1000,1); % Using the data from Norwegian Univ. of Science and Technology
%% Problem iii) Calculating the parameters for each of the distributions 

% Gaussian distribution

mu_G = sum(data)/length(data);
tmp=0;
for i = 1:length(data)
    tmp = tmp + power(data(i) - mu_G,2);
end
ro_G = tmp/length(data); 

% Rayleigh
ro_R = power(2*length(data),-1) * sum (data.^2);


% Erlang distribution
m = [1,2,3];
lambda_E = (m*length(data))/sum(data);

% Shifted exponential distribution
alpha_exp = min(data);
tmp=0;
for i =1:length(data)
    tmp = tmp + data(i) - alpha_exp; 
end

lambda_exp = length(data)/tmp;


% Shifted Rayleigh
 


%% Problem iv: Superimpose graphs of histograms and marginal denseties

% pdf's:

f_G = (1./sqrt(2.*pi.*ro_G)) .* exp(-power( (1:1000)-mu_G,2)./(2.*ro_G));
f_R = ((1:1000)./ro_R).*exp (-(power((1:1000),2))./(2.*ro_R ));
M = 1;
f_E1 = (power(lambda_E(1),M+1)./(factorial(M))).*power((1:1000),M).*exp(-lambda_E(1).*(1:1000));
M = 2;
f_E2 = (power(lambda_E(2),M+1)./(factorial(M))).*power((1:1000),M).*exp(-lambda_E(2).*(1:1000));
M = 3;
f_E3 = (power(lambda_E(3),M+1)./(factorial(M))).*power((1:1000),M).*exp(-lambda_E(3).*(1:1000));
f_exp = zeros(1,1000);
for i=1:1000
    if i< alpha_exp
        f_exp(i)=0;
    elseif i>=alpha_exp
        f_exp(i) = lambda_exp .* exp(-lambda_exp.*(i-alpha_exp));
    end
end
alpha_SR = min(data)+0.01;
f_SR = zeros(1,1000);
tmp1=0;
tmp2=0;

for j=1:1000
    tmp1 = tmp1 + data(j)-alpha_SR;
    tmp2 = tmp2 + 1/(data(j)-alpha_SR);
end

ro_SR=tmp1/tmp2;

for i =1:1000
    if i<alpha_SR
        f_SR(i)=0;
    elseif i>=alpha_SR
        f_SR(i)= ((i-alpha_SR)./ro_SR)*exp(-((i-alpha_SR).^2)./(2*ro_SR)); 
    end
end


figure; 
% Plotting histograms and margnial densisies

subplot(4,2,1)
hist(res(:,1), 500);
xlabel 'msec'
ylabel '#ping'
title 'mercury.iet.ntnu.no - Located in Trondheim, Norway';
subplot(4,2,2)
plot(f_G)
xlim([0 max(data)+100])
ylabel 'Pr'
title  'Gaussian distribution'
subplot(4,2,3)
plot(f_E1);
xlim([0 max(data)+100])
xlabel 'm'
title 'Erlang distribution m=1'
subplot(4,2,4)
plot(f_E2)
xlim([0 max(data)+100])
title 'Erlang distribution m=2'
subplot(4,2,5)
plot(f_E3)
xlim([0 max(data)+100])
title 'Erlang distribution m=3'
subplot(4,2,6)
plot(f_exp)
xlim([0 max(data)+100])
title 'Shifted exponential distribution'
subplot(4,2,7)
plot(f_R)
xlim([0 max(data)+100])
title 'Rayleigh distribution'
subplot(4,2,8)
plot(f_SR)
xlim([0 max(data)+100])
title 'Shifted Reyleigh distribution'





%% Problem v) 
%
% Finding the distribution that maximizes the likelihood: 
% Putting in the measured data in the pdf's



f_GI = (1./sqrt(2.*pi.*ro_G)) .* exp(-power( (data)-mu_G,2)./(2.*ro_G));
f_RI = ((data)./ro_R).*exp (-(power((data),2))./(2.*ro_R ));
f_E1I = (power(lambda_E(1),M+1)./(factorial(M))).*power((data),M).*exp(-lambda_E(1).*(data));
M = 2;
f_E2I = (power(lambda_E(2),M+1)./(factorial(M))).*power((data),M).*exp(-lambda_E(2).*(data));
M = 3;
f_E3I = (power(lambda_E(3),M+1)./(factorial(M))).*power((data),M).*exp(-lambda_E(3).*(data));


f_expI = lambda_exp .* exp(-lambda_exp.*(data-alpha_exp));    
alpha_SR = min(data)+0.01; % matlab does not like to divide by zero

tmp1=0;
tmp2=0;
% Finding ro_SR
for j=1:1000
    tmp1 = tmp1 + data(j)-alpha_SR;
    tmp2 = tmp2 + 1/(data(j)-alpha_SR);
end
ro_SR=tmp1/tmp2;
data_SR = zeros(1,length(data));
for i=1:length(data)
    data_SR(i) = data(i)-alpha_SR;
end

f_SRI = ((data_SR)./ro_SR).*exp(-((data_SR).^2)./(2.*ro_SR));      
f_SRI = f_SRI';
dists = [f_GI, f_RI, f_E1I, f_E2I, f_E3I, f_expI, f_SRI];
fists = {'Gaussian', 'Rayleigh', 'Erlang0','Erlang1','Erlang2', 'Expoinential', 'Shifted Rayleigh'};
loglikelihoods(1:7) = sum(log(dists(:,1:7)));

[lhat,idx] = max(loglikelihoods);


% Output shows that the shifted exponential distribution is the one that
% maximizes the likelihood function, and by obersvation this is also the pdf
% that visually fits best. 


%% Bigger plots of the distributions

figure;
hist(res(:,1), 500);
xlabel 'msec'
ylabel '#ping'
title 'mercury.iet.ntnu.no - Located in Trondheim, Norway';
figure;
plot(f_G)
xlim([0 max(data)+100])
ylabel 'Pr'
title  'Gaussian distribution'
figure;
plot(f_E1);
xlim([0 max(data)+100])
xlabel 'm'
title 'Erlang distribution m=1'
figure;
plot(f_E2)
xlim([0 max(data)+100])
title 'Erlang distribution m=2'
figure;
plot(f_E3)
xlim([0 max(data)+100])
title 'Erlang distribution m=3'
figure;
plot(f_exp)
xlim([0 max(data)+100])
title 'Shifted exponential distribution'
figure;
plot(f_R)
xlim([0 max(data)+100])
title 'Rayleigh distribution'
figure;
plot(f_SR)
xlim([0 max(data)+100])
title 'Shifted Reyleigh distribution'





