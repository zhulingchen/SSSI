Contourlet Toolbox (version 2.0, November 2003)
http://www.ifp.uiuc.edu/~minhdo/software/

This toolbox contains Matlab files that implement the contourlet
transform and its utility functions.  

The main functions are the following.  In addition, there are several demos 
(decdemo, nlademo, nlademo2, and denoisedemo) that provide examples on 
how to use these functions.  For an overview description of all functions 
provided in the toolbox, see file Contents.m or type 
	help <contourlet_toolbox_directory> 
at the Matlab command.

pdfbdec:    Contourlet or pyramidal directional filter bank (PDFB)
            decomposition.
            For example, suppose that x is a double matrix, then
                y = pdfbdec(x, '9-7', 'pkva', [0, 3, 3, 4]);
            will decompose x using a PDFB with '9-7' filter for the
            pyramidal decomposition stage and 'pkva' filter for the
            direction decomposition stage.  There will be 4 levels of
            pyramidal decomposition and the numbers directional
            decomposition at each pyramidal level (from coarse to fine)
            are: 0, 3, 3, and 4.  The output y has cell vector structure
            such that:
                y{1} is the lowpass subband image;
                y{2}{1}, y{2}{2}, y{2}{3} are 3 highpass subband images 
                    that together with y{1} is output of 2-D wavelet 
                    filter bank; 
                y{3}{1}, ..., y{3}{8} are 8 directional subband images at
                    the next finer pyramidal level;
                y{4}{1}, ..., y{4}{8} are 8 directional subband images at
                    the next finer pyramidal level;
                y{5}{1}, ..., y{5}{16} are 16 directional subband images at
                    the finest pyramidal level.
            For other available pyramidal and directional filter names,
            see functions dfilters and pfilters, respectively.
            Function pdfbdec can also be used for 2-D wavelet decomposition
            For example,
		y = pdfbdec(x, '9-7', '', zeros(1, 4))
            will decompose x into 4-level wavelet decomposition using 
            '9-7' filter.

pdfbrec:    Pyramidal directional filter bank (PDFB) reconstruction.
            This is the inverse of the last function, i.e.
                xrec = pdfbrec(y, '9-7', 'pkva');
            would reconstruct xrec = x.

showpdfb:   Show the PDFB output.

pdfb_tr:    Truncate the PDFB coefficients to the given subbands or 
            a number of most significant coefficients.

pdfb2vec:   Convert the PDFB output into a vector of coefficients.

vec2pdfb:   Undo the last function.

Note:       There is one mex file (resampc.c) in the Contourlet Toolbox
            that might need to be recompiled.
            This can be done by typing from the Matlab command window
            >> mex resampc.c  

            
References:

[1] M. N. Do and M. Vetterli, "The contourlet transform: an efficient 
    directional multiresolution image representation."
    submitted to IEEE Transactions Image  Processing, 2003.
    
[2] M. N. Do and M. Vetterli, "Contourlets." Beyond Wavelets, 
    G. V. Welland ed., Academic Press, 2003. 

[3] M. N. Do and M. Vetterli, "Contourlets: a directional multiresolution 
    image representation," IEEE International Conference on Image Processing
    Rochester, September 2002.
    
(Check http://www.ifp.uiuc.edu/~minhdo/publications/ for these and 
other related papers.)
