function out = fill_gaps(in_prev, in, in_next, fill_label)

in = in .* fill_label;

ext_len = 3;
filter_size = 1.7;

if (isempty(in_prev))
    % symmetric extension on top
    in = [in(1+ext_len:-1:2, :); in];
else
    in = [in_prev(end-ext_len+1:end, :); in];
end

if (isempty(in_next))
    % symmetric extension on bottom
    in = [in; in(end-1:-1:end-ext_len, :)];
else
    in = [in; in_next(1:ext_len, :)];
end

in = periodic_extensions(in, [0, ext_len]);

fill_label = periodic_extensions(fill_label, [0, ext_len]);
fill_label = symmetric_extensions(fill_label, [ext_len, 0]);

[r_ext, c_ext] = size(in);

% initialize the covariance parameters of the filtering kernel for each pixel location
C00 = ones(size(in));
C01 = zeros(size(in));
C11 = ones(size(in));

% define the size of the filtering kernel. A larger filtering kernel applies more filtering to the data.
h_mat_interpolation = filter_size * ones(size(in));

% evaluate the window size of the filtering kernel
ksize_mat = evaluate_ksize(h_mat_interpolation);

% apply filtering and inpainting using a uniform kernel. 
out_tmp = inpainting_order_zero(in, fill_label, h_mat_interpolation, ksize_mat, C00, C01, C11, ext_len);

% evaluate the image gradients - these are used to evaluate the filtering kernel direction and width
[gx, gy] = gradient(out_tmp);

% estimate the steering parameters using the input gradients
[C00, C01, C11, sigma_1, sigma_2, theta] = estimate_steering_parameters(gy, gx, ones(r_ext, c_ext), ones(r_ext, c_ext), 7, 0.1);

h_mat_interpolation = 2 * ones(size(in));
h_mat_interpolation(fill_label > 0) = filter_size;
ksize_mat = evaluate_ksize(h_mat_interpolation);

% h_mat_interpolation = zeros(r_ext, c_ext);
% h_mat_interpolation(fill_label == 0) = filter_size;
% [~, h_mat_interpolation] = high_pass_filter_data(h_mat_interpolation', 1);
% h_mat_interpolation = h_mat_interpolation';
% ksize_mat = evaluate_ksize(h_mat_interpolation);

out = inpainting_order_zero(in, fill_label, h_mat_interpolation, ksize_mat, C00, C01, C11, ext_len);

out = out(ext_len + 1:end - ext_len, ext_len + 1:end - ext_len);

end