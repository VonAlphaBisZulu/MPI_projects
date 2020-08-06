% creates a column vector of *length, with ones at *indices
function ive = iv( len, indices,nan)
% len: length of vector
% indices: indices of positions that should be 1
% row or column vector
ive = zeros(len,1);
ive(indices) = 1;
if nargin > 2
    ive(ive==0) = nan();
end
end

