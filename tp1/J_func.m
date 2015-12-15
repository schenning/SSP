function out = J_func(matr)
   [r,c] = size(matr);
   out = fliplr(eye(r));
   