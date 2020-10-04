function [ ] = Untitled4( A )
for k=1:size(A,1)
    fprintf('$%8i$ & $%8.2f$ & $%8.2f$ & $%8.2f$  & $%8.2f$ \\\\ \n', A(k,1), A(k,2), A(k,3), A(k,4), A(k,5))
end
end

