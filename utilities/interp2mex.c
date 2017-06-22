//BILININTERP2 2-D interpolation (table lookup)
// ZI = BILININTERP2(Z, XI, YI) interpolates to find ZI,
// the values for the underlying 2-D function Z
// at the points in matrices XI and YI.
// Out of range values are interpolated to nearest neighbor
// Only input matrices of class double is supported
// Other classes will give unpredictable results
//
// This function assumes that the spacing between all points is 1.
//
// https://www.mathworks.com/matlabcentral/newsreader/view_thread/68708

#include "mex.h"

void pixel(double *pimage, int mrows, int ncols, double x, double y,double *val)
{
    double nwpix, nepix, swpix, sepix; //Pixels surrounding point (x,y)
    double ewweight, nsweight; //Variables for linear interpolation
    double northval, southval;
    int ix, iy; //Truncated coordinates (x,y)

    //Subtract x and y coordinates since Matlab uses indexes 1->for images
    //while c uses indexes from 0->
    x--;
    y--;
    
    if ((x>=0) & (y>=0) & (x<=ncols-1) & (y<=mrows-1))
    {
        //Do bilinear interpolation for current pixel (x,y)
        
        //Get closest integer coordinates north west of current coordinates
        ix = (int)(x);
        iy = (int)(y);

        //Extract the north west, north east, south west and south east pixel values
        nwpix = *(pimage+(int)(ix * mrows + iy)); // linear indexing
        nepix = *(pimage+(int)((ix+1) * mrows + iy));
        swpix = *(pimage+(int)(ix * mrows + iy + 1));
        sepix = *(pimage+(int)((ix+1) * mrows + iy + 1));
        
        ewweight = x-ix;
        nsweight = y-iy;
    
        northval = nwpix + ewweight*(nepix-nwpix);
        southval = swpix + ewweight*(sepix-swpix);

        *val = northval + nsweight*(southval-northval);
    }
    else
    {
        // Do nearest neighbor extrapolation outside of the range
        if ((x<0) & ((y>=0) & (y<=mrows-1))) // left side of image
        {
            *val = *(pimage+(int)(y+0.5));
        }
        else if ((x>ncols-1) & ((y>=0) & (y<=mrows-1))) // right size of image
        {
            *val = *(pimage+(int)(y+0.5 + mrows*(ncols-1)));
        }
        else if ((y<0) & ((x>=0) & (x<=ncols-1))) // top of image
        {
            *val = *(pimage+(int)(mrows*( (int)(x+0.5) )));
        }
        else if ((y>mrows-1) & ((x>=0) & (x<=ncols-1))) // bottom of image
        {
            *val = *(pimage+(int)(mrows*( (int)(x+1.5) ) - 1));
        }
        else if ((x<0) & (y<0)) // upper left corner
        {
            *val = *(pimage);
        }
        else if ((x<0) & (y>mrows-1)) // bottom left corner
        {
            *val = *(pimage + (int)(mrows-1));
        }
        else if ((x>ncols-1) & (y<0)) // upper right corner
        {
            *val = *(pimage + (int)(mrows*(ncols-1)));
        }
        else if ((x>ncols-1) & (y>mrows-1)) // bottom right corner
        {
            *val = *(pimage + (int)(mrows*ncols - 1));
        }
        else
        {
            *val = 0; // Default to zero -- this should only happen if there is an error in the code above.
        }
    };
}
void interp2(double *pimage, int mrows, int ncols, double *px,double *py, int prows, int qcols, double *pval)
{
    int row, col, count=0;

    //Loop through each pair of x-y coordinates and lookup pixel value in input image
    for (col = 0; col<qcols; col++)
    {
        for (row = 0; row<prows; row++)
        {
            pixel(pimage, mrows, ncols, *(px+count), *(py+count),(pval+count));
            count++;
        }
    }

}

void mexFunction(
int nlhs, // Number of left hand side (output) arguments
mxArray *plhs[], // Array of left hand side arguments
int nrhs, // Number of right hand side (input) arguments
const mxArray *prhs[] // Array of right hand side arguments
)
{
    double *pimage, *px, *py, *pval; //Pointers to input data
    int mrows, ncols, prows, qcols; //Local vars to hold matrix sizes
    if(nrhs!=3)
        mexErrMsgTxt("Three inputs required");
    if(nlhs>1)
        mexErrMsgTxt("Only one output is supported");

    if (!mxIsDouble(prhs[0]) | !mxIsDouble(prhs[1]) | !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Only matrices of class double are supported");

    //Extract input variables image, and x/y coordinates
    pimage = mxGetPr(prhs[0]);
    px = mxGetPr(prhs[1]);
    py = mxGetPr(prhs[2]);

    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    prows = mxGetM(prhs[1]);
    qcols = mxGetN(prhs[1]);

    //Establish output matrix
    plhs[0] = mxCreateDoubleMatrix(prows,qcols, mxREAL);
    pval = mxGetPr(plhs[0]);

    //Run bilinear interpolation
    interp2(pimage, mrows, ncols, px, py, prows, qcols,pval);
}
