function [bs p q] = hmcc_hd(D,R,eta)
% Spatially decompose domain and range point clouds in preparation to form an
% H-matrix.
%
% [bs p q] = hmcc_hd(D,[R])
%   D is a 3xnD matrix of points and similarly for R. D corresponds to the
% influence points: in other words, the columns of the matrix. R corresponds
% to the influenced points: in other words, the rows. In many cases, D = R,
% in which case you do not have to specify R.
%   bs,p,q are inputs to hm('Compress'). These are the H-matrix blocks and
% the permutations associated with the range and domain spaces.

% AMB ambrad@cs.stanford.edu
  
  if (nargin < 3) eta = 3; end
  if (nargin < 2) R = []; end
  CheckX(D,'D');
  if (~isempty(R)) CheckX(R,'R'); end
  [bs p q] = hmcc_hd_mex(D,R,eta);
end

function CheckX(X,name)
  if (size(X,1) ~= 3)
    error(sprintf('%s is not 3xN',name));
  end
  [foo p] = sort(X(1,:));
  if (any(sum(abs(diff(X(:,p),1,2))) == 0))
    error(sprintf('%s contains duplicate points.',name));
  end
end
