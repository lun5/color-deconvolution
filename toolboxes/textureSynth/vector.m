% [VEC] = vector(MTX)
% 
% Pack elements of MTX into a column vector.  Same as VEC = MTX(:)

function vec = vector(mtx)

vec = mtx(:);
