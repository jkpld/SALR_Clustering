function class_name = get_minimum_uint_class(A)
% GET_MINIMUM_UINT_CLASS Return the class name for the smallest unsigned
% integer class that will hold the input without overflow.
%
% class_name = get_minimum_uint_class(A)

% James Kapaldo
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


%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
