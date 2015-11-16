function output = periodic_extensions(input, extension_size)

output = input;

% extend columns
if extension_size(2) > 0
    % extend right column
    output = [output, input(:, 1:1:extension_size(2))];    
    % extend left column
    output = [input(:, end-(extension_size(2)-1):1:end), output];
end

if extension_size(1) > 0
    output = output.';
    input = output;
    output = [output, input(:, 1:1:extension_size(1))];    
    output = [input(:, end-(extension_size(1)-1):1:end), output];
    output = output.';
end

end