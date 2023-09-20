function [fval, pos] = data_pareto_front(input)
% Source: www.mathworks.com/matlabcentral/fileexchange/45885-find_pareto_frontier
% by Sisi Ma, version 1.0.0.0
%--------------------------------------------------------------------------
% Identifies the pareto frontier of a set of points (this function
% considers smaller values more desirable, if larger values are
% more desirable, use negative of those values in input data)
%--------------------------------------------------------------------------
% Input: input, a matrix, each row corresponds to a point, each column
% corresponds to a dimension
%--------------------------------------------------------------------------
% Outputs:
% (1) fval: matrix, contains values of the point(s) on the pareto frontier.
% (2) pos: logical array with same number of rows as input matrix
%       1 indicates row of input is on pareto frontier
%       0 otherwise
%--------------------------------------------------------------------------
out=[];
input=unique(input,'rows','stable'); %Fixed bug: requires 'stable' option
for i = 1:size(input,1)
    
    c_data = repmat(input(i,:),size(input,1),1);
    t_data = input;
    t_data(i,:) = Inf(1,size(input,2));
    smaller_idx = c_data>=t_data;
    
    idx=sum(smaller_idx,2)==size(input,2);
    if ~nnz(idx)
        out(end+1,:)=input(i,:);
    end
end
pos = ismember(input,out,'rows');
fval = out;
end