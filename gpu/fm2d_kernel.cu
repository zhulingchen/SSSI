#define BLOCK_SIZE 16 // 32x32 grid (1024 threads total)

__global__ void fm2d_kernel(double* fdm1, double *fdm2, double *fdm3,
        double *boundary, double *a, double *b,
        const int nz, const int nx) {
    int ixt = blockIdx.y * blockDim.y + threadIdx.y; //surface dim [x]
    int izt = blockIdx.x * blockDim.x + threadIdx.x; //depth dim [z]
    int ix = ixt;
    int iz = izt;
    int xg = nx/BLOCK_SIZE; // number of windows in x direction
    int zg = nz/BLOCK_SIZE; // number of windows in z direction
    
    if (nx%BLOCK_SIZE > 0) xg = xg+1;
    if (nz%BLOCK_SIZE > 0) zg = zg+1;
    
    // advance the region of interest
    for (int i = 0; i < zg; i++) {
        iz = izt + (BLOCK_SIZE)*i;
        for (int j = 0; j < xg; j++) {
            ix = ixt + (BLOCK_SIZE)*j;
            
            // finite differencing on interior
            //if (iz > 0 && ix > 0 && iz < nz-2 && ix < nx-2) {
            if (iz+1<nz && ix+1<nx && iz-1>=0 && ix-1>=0){
                fdm3[ix*nz+iz] = b[ix*nz+iz]*fdm2[ix*nz+iz]- fdm1[ix*nz+iz] +
                       a[ix*nz+iz]*(fdm2[(ix+1)*nz+iz] + fdm2[(ix-1)*nz+iz] +
                        fdm2[ix*nz+(iz+1)] + fdm2[ix*nz+(iz-1)]);
            }
            
            // finite differencing at ix = 0
            if (ix == 0 && iz < nz) {
                fdm3[iz] = b[iz]*fdm2[iz] - fdm1[iz] +
                        a[iz]*(fdm2[nz+iz] + fdm2[(iz+1)] + fdm2[(iz-1)]);
            }
            
            // finite differencing at ix = nx-1
            if (ix == nx-1 && iz < nz) {
                fdm3[(nx-1)*nz+iz] = b[(nx-1)*nz+iz]*fdm2[(nx-1)*nz+iz] - fdm1[(nx-1)*nz+iz] +
                        a[(nx-1)*nz+iz]*(fdm2[(nx-2)*nz+iz] + fdm2[(nx-1)*nz+iz+1]
                        + fdm2[(nx-1)*nz+iz-1]);
            }
           
            // finite differencing at iz = 0
            if (iz == 0 && ix < nx) {
                fdm3[ix*nz] = b[ix*nz]*fdm2[ix*nz] -  fdm1[ix*nz] +
                        a[ix*nz]*(fdm2[ix*nz+1] + fdm2[(ix+1)*nz] + fdm2[(ix-1)*nz]);
            }
            
            // finite differencing at iz = nz-1
            if (iz == nz-1 && ix < nx) {
                fdm3[ix*nz+(nz-1)]= b[ix*nz+nz-1]*fdm2[ix*nz+nz-1]- fdm1[ix*nz+nz-1] +
                        a[ix*nz+nz-1]*(fdm2[ix*nz+(nz-2)] + fdm2[(ix+1)*nz+nz-1] +
                        fdm2[(ix-1)*nz+nz-1]);
            }
        }
    }
    // finite differencing at four corners [0][0],[nz-1][0],[0][nx-1],[nz-1][nx-1]
    if (iz == 0 && ix == 0)
        fdm3[0] = b[0]*fdm2[0] -fdm1[0] + a[0]*(fdm2[1] + fdm2[nz]);
    if (iz == nz-1 && ix == 0)
        fdm3[nz-1] = b[nz-1]*fdm2[nz-1] -fdm1[nz-1] + 
                a[nz-1]*(fdm2[nz+nz-1] +fdm2[nz+(nz-2)]);
    if (iz == 0 && ix == nx-1)
        fdm3[(nx-1)*nz] = b[(nx-1)*nz]*fdm2[(nx-1)*nz] -fdm1[(nx-1)*nz] + 
                a[(nx-1)*nz]*(fdm2[(nx-1)*nz] + fdm2[(nx-1)*nz+1]);
    if (iz == nz-1 && ix == nx-1)
        fdm3[(nx-1)*nz+(nz-1)] = b[(nx-1)*nz+(nz-1)]*fdm2[(nx-1)*nz+(nz-1)] -fdm1[(nx-1)*nz+(nz-1)] +
            a[(nx-1)*nz+(nz-1)]*(fdm2[(nx-1)*nz+(nz-2)] +fdm2[(nx-1)*nz+(nz-1)]);
            
    __syncthreads();
    
    for (int i = 0; i < zg; i++) {
        iz = izt + (BLOCK_SIZE)*i;
        for (int j = 0; j < xg; j++) {
            ix = ixt + (BLOCK_SIZE)*j;
            
            // update fdm for next time iteration
            if (iz < nz && ix < nx) {
                fdm1[ix*nz+iz] = fdm2[ix*nz+iz];
                fdm2[ix*nz+iz] = fdm3[ix*nz+iz];
            }
            
            // apply absorbing boundary conditions to 3 sides [not surface]
            if (ix >= 0 && ix < 20){
                fdm1[ix*nz+iz] = boundary[ix]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[ix]*fdm2[ix*nz+iz];
            }
            
            if (ix >= nx-20 && ix < nx) {
                fdm1[ix*nz+iz] = boundary[nx-1-ix]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[nx-1-ix]*fdm2[ix*nz+iz];
            }
            
            if (iz >= nz-20 && iz < nz) {
                fdm1[ix*nz+iz] = boundary[nz-1-iz]*fdm1[ix*nz+iz];
                fdm2[ix*nz+iz] = boundary[nz-1-iz]*fdm2[ix*nz+iz];
            }     
       }
    }
}