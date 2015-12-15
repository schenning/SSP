    function signa_f = levinson2 (ryy)

    LEN = 20;

    if length(ryy) >= LEN
     ryy = ryy(1:LEN);
    end


    %Initialization
    A = zeros(1,LEN);
    K = zeros(1,LEN);
    delta = zeros(1,LEN);
    sigma_f = zeros(1,LEN);
  
    A(1)=1;
    K(1) = -ryy(2) / ryy(1); 
    delta(1) = ryy(1);
    
    
    for i =1:LEN-1
        A(i+1) = [A(i) ; 0] + K(1) * [ 0 ; J_func(A(i))*A(i)];
        delta(i+1)  = ryy(i:-1:1)*A(i);
        K(n+1) = - (delta(n+1))/(sigma_f);
        sigma_f(i+1) = (1-power(K(i+1)));
        
        
        
    end
    
    
        
    
        
        
        
    
    
    

    