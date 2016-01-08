function [val,idx] = henning(seq)

% idea: div input seq in two
% search for the first peak. 

s1 = seq(1:end/2);
s2 = seq(end/2+1:end); 

[val,idx] = max(abs(s1)); 
end
	

