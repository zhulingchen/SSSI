#define BLOCK_SIZE 16 // 32x32 grid (1024 threads total)

__global__ void rtm2d_kernel(double* fdm1, double *fdm2, double *fdm3,
        double *boundary, double *a, double *b,
        const int nz, const int nx, int bz, int it,
        double *data, const int nt) {
    
    int ixt = blockIdx.y * blockDim.y + threadIdx.y; //surface dim [x]
    int izt = blockIdx.x * blockDim.x + threadIdx.x; //depth dim [z]
    int ix = ixt;
    int iz = izt;
    int ez = bz;
    int xg = nx/BLOCK_SIZE; // number of windows in x direction
    int zg = nz/BLOCK_SIZE; // number of windows in z direction
    
    if (nx%BLOCK_SIZE > 0) xg = xg+1;
    if (nz%BLOCK_SIZE > 0) zg = zg+1;
    
    for (int i = 0; i < zg; i++) {
        iz = izt + (BLOCK_SIZE)*i;
        for (int j = 0; j < xg; j++) {
            ix = ixt + (BLOCK_SIZE)*j;
            
            // apply absorbing boundary condition on left
            if (ix >= 0 && ix < 20 && iz < nz) {
                fdm1[ix*nz+iz] = boundary[ix]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[ix]*fdm2[ix*nz+iz];
            }
            
            // apply absorbing boundary condition on right
            if (ix >= nx-20 && ix < nx && iz < nz) {
                fdm1[ix*nz+iz] = boundary[nx-1-ix]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[nx-1-ix]*fdm2[ix*nz+iz];
            }
            
            // apply absorbing boundary condition at depth
            if (bz>=nz-20 && iz >= nz-20 && iz < nz) {
                fdm1[ix*nz+iz] = boundary[nz-1-iz]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[nz-1-iz]*fdm2[ix*nz+iz];
            }
        } // j loop
    }// i loop
    
    //__syncthreads();
    
    for (int i = 0; i < zg; i++) {
        iz = izt + (BLOCK_SIZE)*i;
        for (int j = 0; j < xg; j++) {
            ix = ixt + (BLOCK_SIZE)*j;
            
            // computing grid depth (extent in z to solve)
            if (bz == nz)
                ez = nz-2;
            else
                ez = bz;
            
            // time extrapolation between iz and bz
            if (iz < bz && ix < nx)
                fdm3[ix*nz+iz] = fdm3[ix*nz+iz] - fdm1[ix*nz+iz];
            
            //time extrapolation over interior
            if(iz > 0 && iz < ez && ix > 0 && ix < nx-1)
                fdm2[ix*nz+iz] = b[ix*nz+iz]*fdm1[ix*nz+iz] + fdm2[ix*nz+iz]
                        + a[(ix+1)*nz+iz]*fdm1[(ix+1)*nz+iz]
                        + a[(ix-1)*nz+iz]*fdm1[(ix-1)*nz+iz]
                        + a[ix*nz+iz+1]*fdm1[ix*nz+iz+1]
                        + a[ix*nz+iz-1]*fdm1[ix*nz+iz-1];
            
            // time extrapolation at iz = 0
            if (iz == 0 && ix > 0 && ix < nx-1)
                fdm2[ix*nz] = b[ix*nz]*fdm1[ix*nz] + fdm2[ix*nz]
                        + a[(ix+1)*nz]*fdm1[(ix+1)*nz]
                        + a[(ix-1)*nz]*fdm1[(ix-1)*nz]
                        + a[ix*nz+1]*fdm1[ix*nz+1];
            
            if (iz > 0 && iz < ez && ix == 0)
                // time extrapolation at ix = 0
                fdm2[iz] = b[iz]*fdm1[iz] + fdm2[iz]
                        + a[nz+iz]*fdm1[nz+iz]
                        + a[iz+1]*fdm1[iz+1]
                        + a[iz-1]*fdm1[iz-1];
            
            if (iz > 0 && iz < ez && ix == nx-1)
                //time extrapolation at ix = nx-1
                fdm2[(nx-1)*nz+iz] = b[(nx-1)*nz+iz]*fdm1[(nx-1)*nz+iz] + fdm2[(nx-1)*nz+iz]
                        + a[(nx-2)*nz+iz]*fdm1[(nx-2)*nz+iz]
                        + a[(nx-1)*nz+iz+1]*fdm1[(nx-1)*nz+iz+1]
                        + a[(nx-1)*nz+iz-1]*fdm1[(nx-1)*nz+iz-1];
            
            if (bz == nz) {
                if (iz == nz-1 && ix > 0  && ix < nx)
                    // time extrapolation at iz = nz-1
                    fdm2[ix*nz+nz-1] = b[ix*nz+nz-1]*fdm1[ix*nz+nz-1] + fdm2[ix*nz+nz-1]
                            + a[(ix+1)*nz+nz-1]*fdm1[(ix+1)*nz+nz-1]
                            + a[(ix-1)*nz+nz-1]*fdm1[(ix-1)*nz+nz-1]
                            + a[ix*nz+nz-2]*fdm1[ix*nz+nz-2];
                
                if (iz == nz-1 && ix == 0)
                    // time extrapolation at corner (nz-1,0)
                    fdm2[nz-1] = b[nz-1]*fdm1[nz-1] + fdm2[nz-1]
                            + a[nz+nz-1]*fdm1[nz+nz-1] + a[nz-2]*fdm1[nz-2];
            }
            
            if (iz == 0 && ix == 0)
                // time extrapolation at corner (0,0)
                fdm2[0] = b[0]*fdm1[0] + fdm2[0]
                        + a[nz]*fdm1[nz] + a[1]*fdm1[1];
            
            if (iz == 0 && ix == nx-1)
                // time extrapolation at corner (0,nx-1)
                fdm2[(nx-1)*nz] = b[(nx-1)*nz]*fdm1[(nx-1)*nz] + fdm2[(nx-1)*nz]
                        + a[(nx-2)*nz]*fdm1[(nx-2)*nz]
                        + a[(nx-1)*nz+1]*fdm1[(nx-1)*nz+1];
        } // j loop
    } // i loop
    
//__syncthreads();
    
    for (int i = 0; i < zg; i++) {
        iz = izt + (BLOCK_SIZE)*i;
        for (int j = 0; j < xg; j++) {
            ix = ixt + (BLOCK_SIZE)*j;
            
            if (ix<nx && iz<nz) {
                // set up fdm for next iteration
                fdm1[ix*nz+iz] = fdm2[ix*nz+iz];
                fdm2[ix*nz+iz] = fdm3[ix*nz+iz];
                
                
                // insert surface boundary wavefield
                if (it > 2) {
                    if (iz>0)    fdm3[ix*nz+iz] = 0;
                    if (iz == 0) fdm3[ix*nz] = data[(it-3)*nx+ix];
                }
            }
            __syncthreads();
        } // j loop
    } // i loop
} //rtm2d_step