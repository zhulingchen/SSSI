function output = symmetric_extensions(input, extension_size)

output = input;

% extend columns
if extension_size(2) > 0
    % extend right column
    output = [output, input(:, end-1:-1:end-extension_size(2))];
    % extend left column
    output = [input(:, 1+extension_size(2):-1:2), output];
end

if extension_size(1) > 0
    output = output.';
    output = [output, output(:, end-1:-1:end-extension_size(1))];
    output = [output(:, 1+extension_size(1):-1:2), output];
    output = output.';
end

end