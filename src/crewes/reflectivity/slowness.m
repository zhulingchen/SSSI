function [pslowness,sslowness,ss_freq] = slowness(n_Layers,Qp,Qs,wc_norm,phase,vp,vs)

% computetion of the squared complex frequency dependent velocities and slowness squared,

%for i = 1 : n_Layers
i = sqrt(-1);
 for n=1:n_Layers
    
    % p-wave
    v = vp(n);
    Qp1 = v / (pi * Qp(n));    
    Qp2 = v / (2 * Qp(n));
    r = v + Qp1 * log(wc_norm);  % Muller (1985) Equation 132% a trick with "w_complex=A_0*exp(j*angle)"
    %im = Qp2 - Qp1 * phase;
    im = Qp2 + Qp1 * phase; %checked with Eqn. 132 should be "+"
    aux1 = r * r;
    aux2 = im * im;
    aux = 1/((aux1 + aux2) ^ 2);
    aux3 = (aux1 - aux2) * aux;
    pslowness(n) = complex(aux3, -2 * r * im * aux);
    
    %for s-wave
    v = vs(n);
    Qs1 = v / (pi * Qs(n));
    Qs2 = v / (2 * Qs(n));
    r = v + Qs1 * log(wc_norm);  % Muller (1985) Equation 132
    im = Qs2 + Qs1 * phase;
    
    ss_freq(n) = complex(r * r - im * im,2 * r * im);
    sslowness(n) = 1 / ss_freq(n);
    
end