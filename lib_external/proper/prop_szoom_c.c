/* Routines to calculate magnification (using Lanczos interpolation) of arrays
* 2008 Sep 22  jek  Created prop_szoom_c and roundval routines
* 2014 Aug 21  gmg  Altered and added code for use within Matlab
*/

#include <stdio.h>
#include <math.h>
#include "mex.h"

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#define  K   13
#define  DK   6

/*--------------------------------------------------------------------*/
static double roundval( double x )
{
  if ( x > 0.0 ) 
    return( floor(x+0.5) );
  else
  return( -floor(-x+0.5) );
}

/*--------------------------------------------------------------------*/
void prop_szoom
(
  double        *image_in,      /* square 2D input array */
  int           size_in,        /* size of input array (pixels) */
  double        *image_out,     /* square 2D output array */
  int           size_out,       /* size of output array (pixels) */
  double        magn            /* magnification */
)
{
  double        val;
  double        x;
  double        y;
  double        x_in;
  double        x_phase;
  double        y_in;
  double        **sinc_table;
  int           i;
  int           ik;
  int           ikx;
  int           iky;
  int           ix;
  int           iy;
  int           x1;
  int           x2;
  int           x_pix;
  int           x_out;
  int           y1;
  int           y2;
  int           y_pix;
  int           y_out;

/* Precompute table of sinc kernel coefficients.  Because this routine *
 * only expands or contracts the square image symmetrically about the  *
 * center, just the kernel components for one axis are needed.	       */

  sinc_table = (double **)malloc( size_out * sizeof(double *) );
  for ( i = 0; i < size_out; ++i )
    sinc_table[i] = (double *)malloc( K * sizeof(double) );

  for ( x_out = 0; x_out < size_out; ++x_out )
  {
    x_in = (x_out - size_out / 2) / magn;
    x_phase = x_in - roundval(x_in);
    for ( ik = 0; ik < K; ++ik )
    {
      x = (ik - K / 2) - x_phase;
      if ( fabs(x) <= DK )
      {
        if ( x != 0.0 )
        {
          x = x * 3.141592653589793;
          sinc_table[x_out][ik] = sin(x) / x * sin(x / DK) / (x / DK);
        }
        else
          sinc_table[x_out][ik] = 1.0;
      }
      else
      {
        sinc_table[x_out][ik] = 0.0;
      }
    }
  }

  for ( y_out = 0; y_out < size_out; ++y_out )
  {
    y_in = (y_out - size_out / 2) / magn;
    y_pix = roundval(y_in) + size_in / 2;
    y1 = y_pix - K / 2;
    y2 = y_pix + K / 2;
    if ( (y1 < 0) || (y2 >= size_in) )
      continue;

    for ( x_out = 0; x_out < size_out; ++x_out )
    {
      x_in = (x_out - size_out / 2) / magn;
      x_pix = roundval(x_in) + size_in / 2;
      x1 = x_pix - K / 2;
      x2 = x_pix + K / 2;
      if ( (x1 < 0) || (x2 >= size_in) )
        continue;

      val = 0.0;
      iky = 0;
      for ( iy = y1; iy <= y2; ++iy )
      {
        ikx = 0;
        for ( ix = x1; ix <= x2; ++ix )
        {
          val = val + image_in[iy*(long)size_in+ix] * sinc_table[y_out][iky] * sinc_table[x_out][ikx];
          ++ikx;
        }
        ++iky;
      }
      image_out[y_out*(long)size_out+x_out] = val;
    }
  }

  for ( i = 0; i < size_out; ++i )
    free( sinc_table[i] );
  free( sinc_table );

  return;
} /* prop_szoom */

void mexFunction
(
  int           nlhs,           /* number of expected output mxArrays */
  mxArray       *plhs[],        /* array of pointers to the output mxArrays */
  int           nrhs,           /* number of input mxArrays */
  const mxArray *prhs[]         /* array of pointers to the input mxArrays */
)
{
  double        *image_in;      /* square 2D input array */
  int           size_in;        /* size of input array (pixels) */
  double        *image_out;     /* square 2D output array */
  int           size_out;       /* size of output array (pixels) */
  double        magn;           /* magnification */
  int           ia;             /* index of arguments */

/* Examine input (right-hand-side) arguments */
/*mexPrintf("There are %d right-hand-side argument(s).\n", nrhs);
  for (ia = 0; ia < nrhs; ia++)
  {
    mexPrintf("\tInput Arg %i is of type:\t%s\n", ia, mxGetClassName(prhs[ia]));
  }
*/

/* Examine output (left-hand-side) arguments */
/*mexPrintf("\nThere are %d left-hand-side argument(s).\n", nlhs);
*/

/* Get a pointer to the real data in the input array */
  image_in    = mxGetPr(prhs[0]);
/* Get the values of the scalar inputs */
  size_in     = mxGetScalar(prhs[1]);
  size_out    = mxGetScalar(prhs[2]);
  magn        = mxGetScalar(prhs[3]);
/*mexPrintf("    size_in:%4d   size_out:%4d          magnification:%8.3f\n",
    size_in, size_out, magn);
*/

/*  Create the output array */
  plhs[0] = mxCreateDoubleMatrix( (int)size_out, (int)size_out, mxREAL );
/*  Get a pointer to the real data in the output array */
  image_out   = mxGetPr(plhs[0]);

  prop_szoom( image_in, size_in, image_out, size_out, magn );
}
