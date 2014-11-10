function[rp] = free_surface(a,b,SSlowness,uuC2,uuC,uC,wC,zs,alpha,beta)

%This routine computes the free-surface boundary conditions for 
%incident P and SV waves. 

% vertical slowness of the first layer
am = a(1);
bm = b(1);
a = alpha(1);
b= beta(1);
amI = complex(-imag(am),real(am));
bmI = complex(-imag(bm),real(bm));

% auxiliar quantities
aux1 = SSlowness(1) - uuC2;
aux12 = aux1 * aux1;
ambm = am * bm;
ambm4uu = 4 * ambm * uuC;

den = aux12 + ambm4uu;
den = 1 / den;   

% the coefficients
auxm1 = ambm4uu - aux12;
rp(1,1) = auxm1 * den;          % Rpp

auxm1 = bm * aux1;
auxm3 = 4 * auxm1 * uC;
auxm3 = auxm3*b/a;
rp(1,2) = den * auxm3;          % Rsp

auxm1 = am * aux1;
auxm3 = 4 * auxm1 * uC;
auxm3 = auxm3*a/b;
rp(2,1) = den * auxm3;          % Rps

rp(2,2) = -1.0*rp(1,1);             % Rss

% the phase shift matrix above the source (Muller (1985), Equation (23))
wThick = -  2 * wC * zs;
E(1,1) = exp(amI * wThick);
E(1,2) = exp((amI + bmI) * wThick * 0.5);
E(2,1) = E(1,2);
E(2,2) = exp(bmI * wThick);

% applying phase-shift (Muller (1985), Equation (22))
rp = rp .* E;