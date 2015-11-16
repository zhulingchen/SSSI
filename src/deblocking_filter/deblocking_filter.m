function out = deblocking_filter(in, deblocking_label)

[r, c] = size(in);

h_label_ext_len = 1;
v_label_ext_len = 1;
filter_size = 1.25;

% pad deblocking label
deblocking_label = [zeros(h_label_ext_len, c); deblocking_label; zeros(h_label_ext_len, c)];
deblocking_label = [zeros(r + 2 * h_label_ext_len, v_label_ext_len), deblocking_label, zeros(r + 2 * h_label_ext_len, v_label_ext_len)];

[r_label, c_label] = find(deblocking_label);

for ii = 1:v_label_ext_len
    idx_left = (c_label-ii-1) * (r + 2 * h_label_ext_len) + r_label;
    deblocking_label(idx_left) = 1;
end

for ii = 1:v_label_ext_len
    idx_right = (c_label+ii-1) * (r + 2 * h_label_ext_len) + r_label;
    deblocking_label(idx_right) = 1;
end

for ii = 1:h_label_ext_len
    idx_up = (c_label-1) * (r + 2 * h_label_ext_len) + r_label - ii;
    deblocking_label(idx_up) = 1;
end

for ii = 1:h_label_ext_len
    idx_down = (c_label-1) * (r + 2 * h_label_ext_len) + r_label + ii;
    deblocking_label(idx_down) = 1;
end

deblocking_label = deblocking_label(1+h_label_ext_len:end-h_label_ext_len, 1+v_label_ext_len:end-v_label_ext_len);

image_ext_len = 1;
in = periodic_extensions(in, [0, image_ext_len]);
in = symmetric_extensions(in, [image_ext_len, 0]);
deblocking_label = periodic_extensions(deblocking_label, [0, image_ext_len]);
deblocking_label = symmetric_extensions(deblocking_label, [image_ext_len, 0]);

[r_ext, c_ext] = size(in);

% initialize the covariance parameters of the filtering kernel for each pixel location
C00 = ones(r_ext, c_ext);
C01 = zeros(r_ext, c_ext);
C11 = ones(r_ext, c_ext);

% define the size of the filtering kernel. A larger filtering kernel applies more filtering to the data.
h_mat_interpolation = filter_size * ones(r_ext, c_ext);

% evaluate the window size of the filtering kernel
ksize_mat = evaluate_ksize(h_mat_interpolation);

% apply filtering and inpainting using a uniform kernel.
out_tmp = inpainting_order_zero(in, ones(r_ext, c_ext), h_mat_interpolation, ksize_mat, C00, C01, C11, image_ext_len);

% evaluate the image gradients - these are used to evaluate the filtering kernel direction and width
[gx, gy] = gradient(out_tmp);

% estimate the steering parameters using the input gradients
[C00, C01, C11, sigma_1, sigma_2, theta] = estimate_steering_parameters(gy, gx, ones(r_ext, c_ext), ones(r_ext, c_ext), image_ext_len, 0.1);

h_mat_interpolation = zeros(r_ext, c_ext);
h_mat_interpolation(deblocking_label > 0) = filter_size;
[~, h_mat_interpolation] = high_pass_filter_data(h_mat_interpolation', 1);
h_mat_interpolation = h_mat_interpolation';
ksize_mat = evaluate_ksize(h_mat_interpolation);
out = inpainting_order_zero(in, ones(r_ext, c_ext), h_mat_interpolation, ksize_mat, C00, C01, C11, image_ext_len);

out = out(1+image_ext_len:end-image_ext_len, 1+image_ext_len:end-image_ext_len);

end
