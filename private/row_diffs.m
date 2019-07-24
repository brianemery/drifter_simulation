function dx = row_diffs(x)
% ROW DIFFS.M - differences of all possible row combinations
% dx = row_diffs(x)
% matrix of differences of all possible row combinations
%
% example:
% dx = row_diffs([64 8 2]'*ones(1,10))

% Copyright (C) 2009-2010 Brian Emery 
% 5 Nov 2009 from compute_variance_vs_separation.m

% inputs must have at least 2 rows
if size(x,1) < 2
    disp([mfilename ': not enough rows'])
    dx = NaN; return
end

% Row indexing of x for all row combinations
% this is combinatory logic!
rc = nchoosek(1:size(x,1),2);

% create matrix of diffs
dx = x(rc(:,1),:)-x(rc(:,2),:);

end