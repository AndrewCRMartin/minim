/*************************************************************************
   File:     Emin.c
   Function: Simple demo of energy minimization by Steepest Descents
             in internal coordinates
   Date:     02.07.25
   Author:   Andrew C.R. Martin
   Licence:  MIT
*************************************************************************/

/* Includes                                                             */
#include <stdio.h>
#include <math.h>

/* Structure definitions                                                */
typedef struct _eval
{
   double E;
   double grad_rOH;
   double grad_thetaHOH;
} EVAL;

typedef struct _conf
{
   double rOH;
   double thetaHOH;
} CONF;

/* Prototypes                                                           */
EVAL *SDMinimization(CONF *conf, double c, int maxsteps);
EVAL *EnergyAndGradients(CONF *conf);

/************************************************************************/
/* Main function                                                        */
int main(int argc, char **argv)
{
   CONF conf;
   EVAL *eval;
   
   double stepSize = 0.005;
   int    maxsteps = 20;

   /* starting geometry*/
   conf.rOH      = 10.0;
   conf.thetaHOH = 180.0;

   eval = SDMinimization(&conf, stepSize, maxsteps);

   printf("%f %f %f\n", eval->E,conf.rOH,conf.thetaHOH);

   return(0);
}

/************************************************************************/
/* Steepest Descent Minimization
   -----------------------------
   I/O:     conf       Pointer to a CONF structure containing the
                       conformation
   Input:   stepSize   Minimization step size. Doesn't have to be constant
            maxsteps   Maximum number of steps
   Returns:            Pointer to a structure containing the energy and
                       the gradients

   02.07.25 Original   By: ACRM
*/
EVAL *SDMinimization(CONF *conf, double stepSize, int maxsteps)
{
   int i;
   EVAL *eval;
   
   for(i=0; i<maxsteps; i++)
   {
      eval = EnergyAndGradients(conf);
      if ((fabs(eval->grad_rOH) > (0.001/stepSize)) ||
          (fabs(eval->grad_thetaHOH) > (0.01/stepSize)))
      {
         conf->rOH      -= stepSize*eval->grad_rOH;
         conf->thetaHOH -= stepSize*eval->grad_thetaHOH;
      }
      else
      {
         break;
      }
   }

   if ((fabs(eval->grad_rOH) > (0.001/stepSize)) ||
       (fabs(eval->grad_thetaHOH) > (0.01/stepSize)))
      printf("not converged\n");
   else
      printf("converged in %d steps\n", i);

   return(eval);
}
       

/************************************************************************/
/* Energy and graident calculation
   -------------------------------
   Input:   conf       Pointer to a CONF structure containing the
                       conformation
   Returns:            Pointer to a structure containing the energy and
                       the gradients

   Note that the eval structure is declared as static so that it is a
   fixed location on the heap rather than an automatic variable on the
   stack (in which case the pointer would be invalid when the functio
   returns).

   Energy function:
   E = 2 * k_r((r - r_o)^2) + k_a((a - a_o)^2)
   where k_r is the radius (O-H distance) constant
         r   is the current radius
         r_o is the optimum radius
         k_a is the (H-O-H) angle constant
         a   is the current angle
         a_o is the optimum angle
   The multiplication by 2 is because there are two O-H bonds

   The gradients are given by differentiating the energy terms:
   - WRT the radius
     Gr = d(k_r((r - r_o)^2))/dr = 2*k_r(r - r_o)
   - WRT the angle
     Ga = d(k_a((a - a_o)^2))/dr = 2*k_a(a - a_o)

   In cartesian space, this would need partial differentials wrt x and y
   (and z if in 3 dimensions)
     
   02.07.25 Original   By: ACRM
*/
EVAL *EnergyAndGradients(CONF *conf)
{
   double k_rOH        = 50.0;
   double opt_rOH      = 0.95;
   double k_thetaHOH   = 50.0;
   double opt_thetaHOH = 104.5;

   static EVAL eval;
   eval.E = 2 *
            k_rOH *
            (conf->rOH - opt_rOH)*(conf->rOH - opt_rOH)
            +
            k_thetaHOH *
            (conf->thetaHOH - opt_thetaHOH)*(conf->thetaHOH - opt_thetaHOH);

   eval.grad_rOH      = 2*k_rOH*(conf->rOH - opt_rOH);

   eval.grad_thetaHOH = 2*k_thetaHOH*(conf->thetaHOH - opt_thetaHOH);

   return (&eval);
}
