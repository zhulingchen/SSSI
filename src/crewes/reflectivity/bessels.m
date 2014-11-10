function[J0,J1] = bessels(arg,rec,v1,v2,window_S)
        % THE BESSEL FUNCTION COMPUTATION
        
            if arg <= 2.75
                arg3 = arg / 3;
                x2 = arg3 * arg3;
                x4 = x2 * x2;
                x6 = x2 * x4;
                x8 = x4 * x4;
                x10 = x4 * x6;
                x12 = x6 * x6;
                     
                J0 = 1. - 2.2499997 * x2 + 1.2656208 * x4 - 0.3163866 * x6 + 0.0444479 * x8 - 0.0039444 * x10 + 0.00021 * x12;
                J1 = (0.5 - 0.56249985 * x2 + 0.21093573 * x4 - 0.03954289 * x6 + 0.00443319 * x8 - 0.00031761 * x10 + 0.00001109 * x12) * arg;
            else
                sqrx = sqrt(2. / (pi * arg));
                J0 = sqrx * cos(arg - pi * 0.25);
                J1 = sqrx * cos(arg - pi * 0.75);
            end
