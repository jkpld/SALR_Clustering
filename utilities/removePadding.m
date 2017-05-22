function A = removePadding(A,pad_size)

    sz = size(A);
    N = length(sz);
    
    if N ~= numel(pad_size)
        pad_size = pad_size(1)*ones(1,N);
    end
    
    inds = cell(1,N);
    
    for i = 1:ndims(A)
        inds{i} = (pad_size(i)+1) : (sz(i)-pad_size(i));
    end

    A = A(inds{:});

end