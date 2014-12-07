// eikonal2d_cpp.cpp
// Solution of eikonal equation in 2D rectangular domain for a single source
// at (sz,sx), written in C++ for .mex
//
// V                 velocity model
// dx                grid spacing, assume dx=dy
// sx                x coordinates of the seismic events
// sz                z coordinates of the seismic events
// This matlab source file is free for use in academic research.
// iMax              maximum number of iterations of sweeping (related to the size of model and initial condition)
// All rights reserved.
//
// Written by Lingchen Zhu (zhulingchen@gmail.com)
// Center for Signal and Information Processing, Center for Energy & Geo Processing
// Georgia Institute of Technology
//
// Reference: H. Zhao, A fast sweeping method for Eikonal equations,
// Mathematics of computation, 74(250), pp. 603-627, 2004

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

// function T = eikonal2d(V, dx, sz, sx, iMax)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using namespace std;
    
    const int MAX_VELOCITY = 65535;
    
    // input
    double *p_V = mxGetPr(prhs[0]);
    int dx = (int)*mxGetPr(prhs[1]);
    int sz = (int)*mxGetPr(prhs[2]);
    int sx = (int)*mxGetPr(prhs[3]);
    int iMax = (int)*mxGetPr(prhs[4]);
    
    mwSize nz = mxGetM(prhs[0]);
    mwSize nx = mxGetN(prhs[0]);
    //mexPrintf("dx = %d, sz = %d, sx = %d, iMax = %d, nz = %d, nx = %d\n", dx, sz, sx, iMax, nz, nx);
    
    mxArray *m_Slow = mxCreateDoubleMatrix(nz, nx, mxREAL);
    mxArray *m_TOld = mxCreateDoubleMatrix(nz, nx, mxREAL);
    mxArray *m_TNew = mxCreateDoubleMatrix(nz, nx, mxREAL);
    
    
    double *p_Slow = mxGetPr(m_Slow);
    double *p_TOld = mxGetPr(m_TOld);
    double *p_TNew = mxGetPr(m_TNew);
    
    // ATTENTION: Matlab uses column-major ordering but C/C++ uses row-major ordering
    // so here C/MEX code is written with Matlab's column-major ordering convention
    for (int j = 0; j < nx; j++)
    {
        for (int i = 0; i < nz; i++)
        {
            p_Slow[i+j*nz] = 1 / p_V[i+j*nz];
            p_TOld[i+j*nz] = MAX_VELOCITY;
            p_TNew[i+j*nz] = MAX_VELOCITY;
        }
    }
    p_TOld[(sz-1) + (sx-1)*nz] = 0;
    
    
    double a, b, TT;
    
    for (int iter = 0; iter < iMax; iter++)
    {
        // first sweep
        for (int j = 0; j < nx; j++)
        {
            for (int i = 0; i < nz; i++)
            {
                if (i == 0)
                    // a=Told(2,j);
                    a = p_TOld[1 + j*nz];
                else if (i == nz-1)
                    // a=Told(nz-1,j);
                    a = p_TOld[(nz-2) + j*nz];
                else
                    // a=min(Told(i-1,j),Told(i+1,j));
                    a = min(p_TOld[(i-1) + j*nz], p_TOld[(i+1) + j*nz]);
                
                if (j == 0)
                    // b=Told(i,2);
                    b = p_TOld[i + nz];
                else if (j == nx-1)
                    // b=Told(i,nx-1);
                    b = p_TOld[i + (nx-2)*nz];
                else
                    // b=min(Told(i,j-1),Told(i,j+1));
                    b = min(p_TOld[i + (j-1)*nz], p_TOld[i + (j+1)*nz]);
                
                // if (fabs(a-b)<Slow(i,j)*dx)
                // ATTENTION: use fabs for floating point values
                if (fabs(a-b) < p_Slow[i+j*nz] * dx)
                    // TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
                    TT = 0.5 * (a + b + sqrt( (2 * p_Slow[i+j*nz] * p_Slow[i+j*nz] * dx * dx) - (a-b)*(a-b) ));
                else
                    // TT=min(a,b)+Slow(i,j)*dx;
                    TT = min(a, b) + (p_Slow[i+j*nz] * dx);
                
                // Tnew(i,j) = min(Told(i,j),TT);
                p_TNew[i+j*nz] = min(p_TOld[i+j*nz], TT);
            }
        }
        
        // Tnew(sz,sx)=0;
        p_TNew[(sz-1) + (sx-1)*nz] = 0;
        // Told=Tnew;
        for (int j = 0; j < nx; j++)
            for (int i = 0; i < nz; i++)
                p_TOld[i+j*nz] = p_TNew[i+j*nz];
        
        
        // second sweep
        for (int j = 0; j < nx; j++)
        {
            for (int i = nz-1; i >= 0; i--)
            {
                if (i == 0)
                    // a=Told(2,j);
                    a = p_TOld[1 + j*nz];
                else if (i == nz-1)
                    // a=Told(nz-1,j);
                    a = p_TOld[(nz-2) + j*nz];
                else
                    // a=min(Told(i-1,j),Told(i+1,j));
                    a = min(p_TOld[(i-1) + j*nz], p_TOld[(i+1) + j*nz]);
                
                if (j == 0)
                    // b=Told(i,2);
                    b = p_TOld[i + nz];
                else if (j == nx-1)
                    // b=Told(i,nx-1);
                    b = p_TOld[i + (nx-2)*nz];
                else
                    // b=min(Told(i,j-1),Told(i,j+1));
                    b = min(p_TOld[i + (j-1)*nz], p_TOld[i + (j+1)*nz]);
                
                // if (fabs(a-b)<Slow(i,j)*dx)
                // ATTENTION: use fabs for floating point values
                if (fabs(a-b) < p_Slow[i+j*nz] * dx)
                    // TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
                    TT = 0.5 * (a + b + sqrt( (2 * p_Slow[i+j*nz] * p_Slow[i+j*nz] * dx * dx) - (a-b)*(a-b) ));
                else
                    // TT=min(a,b)+Slow(i,j)*dx;
                    TT = min(a, b) + (p_Slow[i+j*nz] * dx);
                
                // Tnew(i,j) = min(Told(i,j),TT);
                p_TNew[i+j*nz] = min(p_TOld[i+j*nz], TT);
            }
        }
        
        // Tnew(sz,sx)=0;
        p_TNew[(sz-1) + (sx-1)*nz] = 0;
        // Told=Tnew;
        for (int j = 0; j < nx; j++)
            for (int i = 0; i < nz; i++)
                p_TOld[i+j*nz] = p_TNew[i+j*nz];
        
        
        // third sweep
        for (int j = nx-1; j >= 0; j--)
        {
            for (int i = nz-1; i >= 0; i--)
            {
                if (i == 0)
                    // a=Told(2,j);
                    a = p_TOld[1 + j*nz];
                else if (i == nz-1)
                    // a=Told(nz-1,j);
                    a = p_TOld[(nz-2) + j*nz];
                else
                    // a=min(Told(i-1,j),Told(i+1,j));
                    a = min(p_TOld[(i-1) + j*nz], p_TOld[(i+1) + j*nz]);
                
                if (j == 0)
                    // b=Told(i,2);
                    b = p_TOld[i + nz];
                else if (j == nx-1)
                    // b=Told(i,nx-1);
                    b = p_TOld[i + (nx-2)*nz];
                else
                    // b=min(Told(i,j-1),Told(i,j+1));
                    b = min(p_TOld[i + (j-1)*nz], p_TOld[i + (j+1)*nz]);
                
                // if (fabs(a-b)<Slow(i,j)*dx)
                // ATTENTION: use fabs for floating point values
                if (fabs(a-b) < p_Slow[i+j*nz] * dx)
                    // TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
                    TT = 0.5 * (a + b + sqrt( (2 * p_Slow[i+j*nz] * p_Slow[i+j*nz] * dx * dx) - (a-b)*(a-b) ));
                else
                    // TT=min(a,b)+Slow(i,j)*dx;
                    TT = min(a, b) + (p_Slow[i+j*nz] * dx);
                
                // Tnew(i,j) = min(Told(i,j),TT);
                p_TNew[i+j*nz] = min(p_TOld[i+j*nz], TT);
            }
        }
        
        // Tnew(sz,sx)=0;
        p_TNew[(sz-1) + (sx-1)*nz] = 0;
        // Told=Tnew;
        for (int j = 0; j < nx; j++)
            for (int i = 0; i < nz; i++)
                p_TOld[i+j*nz] = p_TNew[i+j*nz];
        
        
        // fourth sweep
        for (int j = nx-1; j >= 0; j--)
        {
            for (int i = 0; i < nz; i++)
            {
                if (i == 0)
                    // a=Told(2,j);
                    a = p_TOld[1 + j*nz];
                else if (i == nz-1)
                    // a=Told(nz-1,j);
                    a = p_TOld[(nz-2) + j*nz];
                else
                    // a=min(Told(i-1,j),Told(i+1,j));
                    a = min(p_TOld[(i-1) + j*nz], p_TOld[(i+1) + j*nz]);
                
                if (j == 0)
                    // b=Told(i,2);
                    b = p_TOld[i + nz];
                else if (j == nx-1)
                    // b=Told(i,nx-1);
                    b = p_TOld[i + (nx-2)*nz];
                else
                    // b=min(Told(i,j-1),Told(i,j+1));
                    b = min(p_TOld[i + (j-1)*nz], p_TOld[i + (j+1)*nz]);
                
                // if (fabs(a-b)<Slow(i,j)*dx)
                // ATTENTION: use fabs for floating point values
                if (fabs(a-b) < p_Slow[i+j*nz] * dx)
                    // TT=0.5*(a+b+sqrt(2*Slow(i,j)^2*dx^2-(a-b)^2));
                    TT = 0.5 * (a + b + sqrt( (2 * p_Slow[i+j*nz] * p_Slow[i+j*nz] * dx * dx) - (a-b)*(a-b) ));
                else
                    // TT=min(a,b)+Slow(i,j)*dx;
                    TT = min(a, b) + (p_Slow[i+j*nz] * dx);
                
                // Tnew(i,j) = min(Told(i,j),TT);
                p_TNew[i+j*nz] = min(p_TOld[i+j*nz], TT);
            }
        }
        
        // Tnew(sz,sx)=0;
        p_TNew[(sz-1) + (sx-1)*nz] = 0;
        // Told=Tnew;
        for (int j = 0; j < nx; j++)
            for (int i = 0; i < nz; i++)
                p_TOld[i+j*nz] = p_TNew[i+j*nz];
        
    }
    
    // output
    plhs[0] = m_TNew;
    
    // ATTENTION: Don't forget to free dynamic memory allocated by MXCREATE* functions (except for output arrays), otherwise memory leak will occur
    mxDestroyArray(m_Slow);
    mxDestroyArray(m_TOld);
}

