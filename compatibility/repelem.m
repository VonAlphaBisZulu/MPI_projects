function U = repelem(V,varargin)
%REPELEM Replicate elements of an array.
%   U = REPELEM(V,N), where V is a vector, returns a vector of repeated
%   elements of V.
%   - If N is a scalar, each element of V is repeated N times.
%   - If N is a vector, element V(i) is repeated N(i) times. N must be the
%     same length as V.
%
%   B = repelem(A, R1, ..., RN), returns an array with each element of A
%   repeated according to R1, ..., RN. Each R1, ..., RN must either be a
%   scalar or a vector with the same length as A in the corresponding
%   dimension.

% check if size is correct
try
    % if V = 1xn or nx1 and N = 1xn
    if size(V,1) == 1 && length(varargin) == 1
        varargin{2} = 1;
        varargin = flip(varargin);
    elseif size(V,2) == 1 && length(varargin) == 1
        varargin{2} = 1;
    end
    % if V = 1xn and N = 1x1
    for i = find(cellfun(@length,varargin) == 1)
        varargin{i} = varargin{i}*ones(1,getelements(size(V)',i));
    end
    % for all cases
    if all( size(V) == cellfun(@length,varargin) )
        [reps{1:numel(varargin)}] = ndgrid(varargin{:});
        U = cell(cellfun(@length,varargin));
        U(:) = arrayfun(@(x) V(x)*ones(cellfun(@(y) y(x),reps)) ,1:numel(V),'UniformOutput',false);
        U = cell2mat(U);
    else
        error('The length of the repition vectors doesn''t match with the size of the original matrix');
    end
catch
    error('Provide as many repition vectors as dimensions of the original matrix');
end
end