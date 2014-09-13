function C = fdct_wrapping_aux(C,nbscales)

% Goes from C++ format to Matlab format, C = C(end:-1:1);

for s=2:nbscales-1
    B = C{s};
    nq = length(B)/4;
    B1 = B(1:nq);    B2 = B(nq+1:2*nq);    B3 = B(2*nq+1:3*nq);    B4 = B(3*nq+1:4*nq);
    C{s} = [B3(end:-1:1), B2(end:-1:1), B1(end:-1:1), B4(end:-1:1)];
  end
  
