function [window_s] = hanning_window(perc_U,u1,u2,delta_U)

% ray-parameter Hanning window, figure 5b from Mallick and Frazer (1987)

wl = (u2 - u1) * perc_U;
p1 = u1;
p2 = u1 + wl;
p3 = u2 - wl;
p4 = u2;
nu = (u2 - u1) / delta_U + 1;
p = u1 - delta_U;
for i = 1 : nu
    p = p + delta_U;
    window_s(i) = 1;
    if p <= p2
       window_s(i) = 0.5 + 0.5 * cos(pi * (p2 - p)/(p2 - p1));
    end
    if p >= p3
       window_s(i) = 0.5 + 0.5 * cos(pi * (p - p3)/(p4 - p3));
    end
end
