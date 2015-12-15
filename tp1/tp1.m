%TD1 SSP 
close all;
% Zero padding
% For any unibased estimator, you can't go below the CRLB


Ts = 1; 


%% Part 1: Spectrum Estimation
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
    
end

% Zero padding gives higher resolution the periograms. 

%c) 

% Estimation of correlation sequence of the signal 
 

ryy = conv(y,y(end:-1:1));
ryy = ryy(length(y):end);
figure;

plot(ryy)

 




