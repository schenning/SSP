function [a,k,Jo]=levdur(r)

	
    % Input: r(m)
	p=length(r)-1;
    a=zeros(p+1,1); % Prediction error filter
    k=zeros(p,1);
    Jo=zeros(p+1,1); 
    % Initialization
    J=r(1); Jo(1)=J;
    beta=r(2); 
	k(1)=-beta/J; 
	a(1)=k(1);
    J=J+beta*k(1); 
	Jo(2)=J;
    % Recursion
    for m=2:p
        beta=(r(2:m))'*flipud(a(1:m-1))+r(m+1);
        k(m)=-beta/J;
        a(1:m)=[(a(1:m-1))' 0]'+[(flipud(a(1:m-1)))' 1]'*k(m);
        J=J+beta*k(m); Jo(m+1)=J;
    end
    a(2:p+1)=a(1:p); a(1)=1;
