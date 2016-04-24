function n = nnnan(x,dim)
% number of not nans (find() finds nans as well)
% optional, find down dimension dim
% N Killian 110907
if nargin<2, dim=[];end
if ~isempty(dim)
    if dim == 2, x = x';dim = 1;
    elseif dim>2, error('only 1 or 2 dim. works');end
        n = zeros(size(x,2),1);
    for k = 1:size(x,2)
        n(k) = length(find(~isnan(x(:,k))));
    end
else
    n = length(find(~isnan(x(:))));
end
n = makerow(n);