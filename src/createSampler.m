% createSampler
% create sampler matrix for each branch in parallel structured sampler
% Author: Lingchen Zhu
% Creation Date: 11/07/2013

function Phi = createSampler(N, M, method, repeat, seed)


if (nargin < 4)
    repeat = false;
end

if (nargin < 5)
    seed = 0;
end

if (repeat)
    rng(seed);
end

switch lower(method)
    case 'rand'
        Phi = sqrt(1/M) * randn(M, N);
    case 'uniform'
        Phi = zeros(M, N);
        R = floor(N / M);
        Phi(:, 1:R:N) = eye(M, M);
    case {'rdemod', 'aic'}
        R = N / M;
        vec = binornd(1, 0.5, [1, N]);
        vec(vec == 0) = -1;
        D = diag(vec);
        if (mod(N, M)) % R is not integer, maybe need more debugging
            [N0, M0] = rat(R);
            H0 = zeros(M0, N0);
            jjtmp = 1;
            for ii = 1:M0
                Rtmp = R;
                if (ii > 1)
                    H0(ii, jjtmp) = 1 - H0(ii-1, jjtmp);
                else
                    H0(ii, jjtmp) = 1;
                end
                Rtmp = Rtmp - H0(ii, jjtmp);
                for jj = jjtmp+1:N0
                    if (Rtmp > 1)
                        H0(ii, jj) = 1;
                    elseif (Rtmp > 0)
                        H0(ii, jj) = Rtmp;
                    else
                        jjtmp = jj-1;
                        break;
                    end
                    Rtmp = Rtmp - H0(ii, jj);
                end
            end
            H = kron(eye(M/M0, N/N0), H0);
        else % R is an integer
            H = kron(eye(M, M), ones(1, R));
        end
        Phi = H * D;
    otherwise
        error('Unknown compressive sensing based sampler');
end