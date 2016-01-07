

    LEN = 20;

    if length(ryy) >= LEN
     ryy = ryy(1:LEN);
    end


    % initialization
    %A = zeros(1,LEN);
    %K = zeros(1,LEN);
    %delta = zeros(1,LEN);
    %sigma_f = zeros(1,LEN);
    AN = zeros(1,LEN);
    A(1) = 1;
    K(1) = -ryy(2) / ryy(1); 
    delta(1) = ryy(1);
    sigma_f(1) = ryy(1) + ryy(2)*A(2);%%?
    S=0;
    Ek= ryy(1);
    for k=2:LEN
        
        
        for j = 1:k
            S = S - A(j)*ryy(k+1-j); %possibly index minstake
        end
        S=S/Ek;
        
        for n = 1:(k+1)/2
            tmp = A(n) + S *A(k+1-n);
            A(k+1-n) = tmp;
        end
        Ek = Ek*(1-power(S,2));
    end
    
        
      
        
        
    
    
    
%     S=0;
%     
%     for n = 2:LEN-1
%         for j=1: n-1
%             S = S + ryy(j)*A(n-j);
%         end
%         S = S + ryy(n);
%         K = -S/sigma_f;
%         for j=1:n-1
%             AN(j) = A(j) * K*A(n-j);
%             disp(AN(j));
%         end
%         AN(n) = K;
%         sigma_f = sigma_f * (1-K^2);
%     end
    
    
        
            
    
%     for i =1:LEN-1
%         A(:,i+1) = [A(:,i) ; 0] + K(i+1) * [ 0 ; J_func(A(:,i)'.*A(:,i)')];
%         delta(i+1)  = ryy(i:-1:1)*A(i);
%         K(i+1) = - (delta(i+1))/(sigma_f(i));
%         sigma_f(i+1) = (1-power(K(i+1),2));
%     end

    
    
        
    
        
        
        
    
    
    

    