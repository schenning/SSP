function [a , e, k] = levinson2 ( r )

	p=length(r)-1;
    a=zeros(p+1,1); % Prediction error filter
    e=zeros(p,1);
    k=zeros(p+1,1); 
    % Initialization
    sigm=r(1); k(1)=sigm;
    delta=r(2); 
	e(1)=-delta/sigm; 
	a(1)=e(1);
    sigm=sigm+delta*e(1); 
	k(2)=sigm;
    % Recursion
    for m=2:p
        delta=(r(2:m))'*flipud(a(1:m-1))+r(m+1);
        e(m)=-delta/sigm;
        a(1:m)=[(a(1:m-1))' 0]'+[(flipud(a(1:m-1)))' 1]'*e(m);
        sigm=sigm+delta*e(m); k(m+1)=sigm;
    end
    a(2:p+1)=a(1:p); a(1)=1;
    
end