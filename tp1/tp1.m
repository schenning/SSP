%TD1 SSP 
close all;
% Zero padding
% For any unibased estimator, you can't go below the CRLB


Ts = 1; 


%% Part I: Spectrum Estimation
% a) 


[y,sigpar] = sig(1024); 
[pxx,w] = periodogram(y);
plot(w,10*log10(pxx));

%b) 

N = 256;
y_b = sig(N);

N_mark = [64, 128, 256, 512, 1024];
w = window('boxcar',N);


%subplot(3,2,1)
for j=1:length(N_mark)
    subplot(3,2,j)
    periodo(y_b,N_mark(j));
    title '' 
end

% Zero padding gives higher resolution the periograms. 

%c) 

% Estimation of correlation sequence of the signal 
 

ryy = conv(y,y(end:-1:1));
ryy = ryy(length(y):end);
figure;

plot(ryy)

%% Part II: Parameter Estimation
%% Problem e)
%-------------------------------------------------------------
%clear all;



close all;

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

  
% let m be 2 x n;
%m  = 10*n;
m=n;
df = 1/m;
% Zero padding 
%yk = [yk; zeros(m,10)];

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
	[l_chapau,l_idx ]  = henning(Yl(:,j));
	f1_chapau(j)          = l_idx/m;
	A1_chapau(j)          = (2/n)* abs(Yl(l_idx,j));
	phi_chapau(j)         = angle(Yl(l_idx,j));

	tmp = 0;
	for k=0:n-1
		  tmp = tmp + power((yk(k+1,j) - A1_chapau(j) * cos(2*pi*f1_chapau(j)*k + phi_chapau(j))),2);
	end
	sigma_chapau(j) = (1/n)*tmp;


	%fprintf('Estimation of parameters: \n\n');
	%fprintf('f1_chapau:     %f\n', f1_chapau(j));
	%fprintf('A1_chapau:     %f\n', A1_chapau(j)); 
	%fprintf('phi_chapau:    %f\n', phi_chapau(j));
	%fprintf('sigma_chapau:  %f\n', sigma_chapau(j));

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


	fprintf('Estimation of parameters: \n\n');
	fprintf('f1_chapau:    mean:  %f,  variance:   %f   CRLB: %f\n', mean(f1), var(f1), CRB_f1);
	fprintf('A1_chapau:    mean:  %f,  variance:   %f   CRLB: %f\n', mean(A1), var(A1), CRB_A1); 
	fprintf('sigma_chapau: mean:  %f,  variance:   %f   CRLB: %f\n', mean(sigma), var(sigma),CRB_sigma);









