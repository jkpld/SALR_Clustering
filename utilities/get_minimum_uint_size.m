function class_name = get_minimum_uint_size(A)

maxValue = max(A(:));

if maxValue <= intmax('uint8')
	class_name = 'uint8';
elseif maxValue <= intmax('uint16')
	class_name = 'uint16';
elseif maxValue <= intmax('uint32')
	class_name = 'uint32';	
elseif maxValue <= intmax('uint64')
	class_name = 'uint64';
elseif mod(maxValue,1) ~= 0
	class_name = 'NONE';
end


end