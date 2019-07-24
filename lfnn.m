function y=lfnn(x)
% LFNN.M
% lfnn(x) - length(find(~isnan(x)))
% Shortcut 
y=num2str(length(find(~isnan(x))));

if nargout==0
    disp(y)
end
end