#include <stdio.h>
#include <math.h>

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


EVAL *sd(CONF *conf, double c, int maxsteps);
EVAL *Eandg(CONF *conf);

int main(int argc, char **argv)
{
   CONF conf;
   EVAL *eval;
   
   double stepSize = 0.005;
   int    maxsteps = 20;

   /* starting geometry*/
   conf.rOH      = 10.0;
   conf.thetaHOH = 180.0;

   eval = sd(&conf, stepSize, maxsteps);

   printf("%f %f %f\n", eval->E,conf.rOH,conf.thetaHOH);

   return(0);
}

EVAL *sd(CONF *conf, double stepSize, int maxsteps)
{
   int i;
   EVAL *eval;
   
   for(i=0; i<maxsteps; i++)
   {
      eval = Eandg(conf);
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
       

EVAL *Eandg(CONF *conf)
{
   double k_rOH = 50.0;
   double opt_rOH = 0.95;
   double k_thetaHOH = 50.0;
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


