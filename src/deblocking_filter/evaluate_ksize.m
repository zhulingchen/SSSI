function ksize_mat = evaluate_ksize(h_mat, distance_to_edge)

%max_distance_to_pixel = 101;
max_distance_to_pixel = 61;

ksize_mat = h_mat*8;
ksize_mat(ksize_mat > max_distance_to_pixel*2+1) = max_distance_to_pixel*2+1;

if exist('distance_to_edge', 'var') == 1
    ind = find(floor((ksize_mat-1)/2) < distance_to_edge);
    ksize_mat(ind) = distance_to_edge(ind)*2+5;
end

ksize_mat = ksize_mat - 3;
ksize_mat = round(ksize_mat/2)*2;
ksize_mat = ksize_mat + 3;