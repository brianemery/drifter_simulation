function sd = compute_sd(x,y)
% COMPUTE SEPERATION DISTANCE
% asd = compute_sd(x,y)
% given x and y positions (arbitrary distance units)
% of drifter clusters, each row a drifter and each column a new time,
% computes the TOTAL seperation distance vs time for each pair.
%
% Similar to relative dispersion calculation, but computes total distance
% rather than x and y distances.
%
% EXAMPLE:
% x = [1 2 3; 1 2 3; 1 2 3];
% y = [1 2 3; 1 4 8; 1 1.5 2];

% Brian Emery Jan '10

% matrix of differences of all possible row combinations, (ie compute 
% separation dist for all possible drifter pairs)
dx = row_diffs(x);
dy = row_diffs(y);

sd = sqrt(dx.^2 + dy.^2);

end
