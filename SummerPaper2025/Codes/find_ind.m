function ind_bin = find_ind(dist_vals, dval)
% FIND_IND Find index j such that dist_vals(j) <= dval(i) < dist_vals(j+1)
%
%   ind_bin = FIND_IND(dist_vals, dval)
%
%   Inputs:
%       dist_vals - Sorted support vector (length N)
%       dval      - Query values (any shape, will be reshaped to column)
%
%   Output:
%       ind_bin - Indices such that for each i:
%                 dist_vals(j) <= dval(i) < dist_vals(j+1)
%
%   Example:
%       dist_vals = [0; 1; 2; 3];
%       dval = [0.2; 1.0; 2.9; 4.0];
%       ind_bin = find_bin_indices(dist_vals, dval)

    % Ensure inputs are column vectors
    dist_vals = dist_vals(:);
    dval = dval(:);
    
    % Add +Inf as the rightmost edge
    edges = [-Inf;dist_vals; Inf];

    % Discretize according to: dist_vals(j) <= x < dist_vals(j+1)
    ind_bin = discretize(dval, edges);
    ind_bin = ind_bin + (ind_bin>=2).*(-1);
end

