%% TP1, Statistical Signal Processing
% Henning Schei

%% Part I: Spectrum Estimation
% a) 
% To obtain a separation, be must have the bandwith $$ B_w $$
%
% $$ B_w  \leq M \\ \\ $$ to 
% the sinusodial frequency separation. 
% $$ B_w = \frac{f_s}{M} \\ $$, 
%
% when $f_s = 1$, 
% 
% $$ M \geq \frac{1}{abs(f_1 - f_2)} $$  
% 
% which gives M = 40
close all
Ts = 1; 
[y,sigpar] = sig(1024); 
[pxx,w] = periodogram(y);
plot(w,10*log10(pxx));

%b) 

N = 256;
y_b = sig(N);
N_mark = [64, 128, 256, 512, 1024];
w = window('boxcar',N);
for j=1:length(N_mark)
    subplot(3,2,j)
    periodo(y_b,N_mark(j)); 
end

% Zero padding gives higher resolution the periograms. As the plots show, 
% it appears in this case to be enough with N' = 512


%c) 

% Estimation of correlation sequence of the signal 
ryy = conv(y,y(end:-1:1));
ryy = ryy(length(y):end);


% d)  
[a,e,K]= levinson2(ryy);
%Peridogram(N=1024) and autoregressive spectrum

%periodo(y,1024);
%hold on
%plot(levinson2(ryy(1:20)),'r');
%hold off


init = [a;zeros(235,1)];
spc  = fft(init);
spc  = 1./spc;
spc  = spc.*conj(spc);
autorg(1:21)= e(1:21).* spc(1:21);
sigma=a;
f = 0:1/40:0.5;
g = autorg/sqrt(21);
g = 10*log10(g.^2);
x = 0:1:20;
subplot(1, 2, 1);
plot(x, sigma(1:21)');
xlabel('k');
ylabel('sigma');
title('Evolution of sigma');

subplot(1, 2, 2);
plot(K);
xlabel('k');
ylabel('PARCORS');
title('evolution of PARCORS');







%% Part II: Parameter Estimation
% Problem e)
%-------------------------------------------------------------



n=32; % number of samples 
f1 = 1/8; phi = 0; A1 = sqrt(2); sigma = 1;
yk = zeros(n,10);
vk = sigma*randn(n,10); % gen noise 
for i = 1:10
	for k = 0:n-1
		yk(k+1,i) = A1*cos(2*pi*f1*k + phi);
	end
end
yk = yk + vk; % add noise
m=n;
df = 1/m;

% Taking the DFT
Yl = zeros(m,10);
for j=1:10
	for l=1:m
		for k=0:m-1
			Yl(l,j) = Yl(l,j) + yk(k+1,j)*exp(1j*2*pi*(l/m)*k);
		end
	end 
end
cnt =1;
f1_chapau    = zeros(1,10); 
A1_chapau    = zeros(1,10);
phi_chapau   = zeros(1,10);
sigma_chapau = zeros(1,10);
  


for j=1:10
	%Maximum Likelihood estimates
	[l_chapau,l_idx ]     = findArgMax(Yl(:,j));
	f1_chapau(j)          = l_idx/m;
	A1_chapau(j)          = (2/n)* abs(Yl(l_idx,j));
	phi_chapau(j)         = angle(Yl(l_idx,j));

	tmp = 0;
	for k=0:n-1
		  tmp = tmp + power((yk(k+1,j) - A1_chapau(j) * cos(2*pi*f1_chapau(j)*k + phi_chapau(j))),2);
	end
	sigma_chapau(j) = (1/n)*tmp;




end



CRB_sigma = 2*power(sigma,4)/n;
CRB_A1    = 2*power(sigma,2)/n; 
CRB_phi   = 8*power(sigma,2)/(n*power(A1,2));
CRB_f1    = 24*power(sigma,2)/(power(n,2) * power(A1,2));
fprintf('Estimation of parameters: ML Estimates \n\n');
fprintf('f1_chapau:    mean:  %f,  variance:   %f  CRLB: %f\n', mean(f1_chapau), var(f1_chapau),CRB_f1);
fprintf('A1_chapau:    mean:  %f,  variance:   %f  CRLB: %f\n', mean(A1_chapau), var(A1_chapau), CRB_A1); 
fprintf('phi_chapau:   mean:  %f,  varinace:   %f  CRLB: %f\n', mean(phi_chapau),var(phi_chapau), CRB_phi);
fprintf('sigma_chapau: mean:  %f,  variance:   %f  CRLB: %f\n', mean(sigma_chapau), var(sigma_chapau), CRB_sigma);
fprintf('\n');   

% The results show that f1_chapau and phi_chapau go below the CRLB 
% and can be considered to be biased, A1_chapau and sigma_chapau is close
% to the CRLB, but results are varying .

% Problem f)
%--------------------------------------------------------------------
% Estimating r0,r1,r2
% 
n=32; % number of samples 
f1 = 1/8; phi = 0; A1 = sqrt(2); sigma = 1;
yk = zeros(n,10);
vk = sigma*randn(n,10); % gen noise 
for i = 1:10
	for k = 0:n-1
		yk(k+1,i) = A1*cos(2*pi*f1*k + phi);
	end
end
yk = yk + vk; % add noise

A1=zeros(1,10);
f1=zeros(1,10);
sigma = zeros(1,10);



for j=1:10	
	N=32;
	r0 = 0;
	r1 = 0;
	r2 = 0;
	for n=1:N-0
		r0 = r0 + yk(n,j)*yk(n,j);	
	end
	for n=1:N-1
		r1 = r1 + yk(n+1,j)*yk(n,j);
	end
	for n=1:N-2
		r2 = r2 + yk(n+2,j)*yk(n,j);
	end

	r0 = r0/length(1:N); r1 = r1/length(1:N-1); r2 = r2/length(1:N-2);

	X  = 0.5*(-r2 + sqrt(power(r2,2) + 8*power(r1,2)));
	A1(j) = sqrt(2*X);
	f1(j) = power(2*pi,-1) * acos(r1/X);
	sigma(j) = r0 -X; 
end


fprintf('Estimation of parameters: Covariance Matching \n\n');
fprintf('f1_chapau:    mean:  %f,  variance:   %f   CRLB: %f\n', mean(f1), var(f1), CRB_f1);
fprintf('A1_chapau:    mean:  %f,  variance:   %f   CRLB: %f\n', mean(A1), var(A1), CRB_A1); 
fprintf('sigma_chapau: mean:  %f,  variance:   %f   CRLB: %f\n', mean(sigma), var(sigma),CRB_sigma);

% The estimates to appear the have the same charasteristics, only with
% a slightly higher variance. Here does the f_









