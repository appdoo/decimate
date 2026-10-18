// ============================================================================
//
// math3d.cpp : globale funktionen
//
// ============================================================================

#include "math3d.h"

char *_domain    = "argument domain error";
char *_sing      = "argument singularity";
char *_overflow  = "overflow range error";
char *_underflow = "underflow range error";
char *_tloss     = "total loss of precision";
char *_ploss     = "partial loss of precision";

bool _math_err;

bool ismatherr()
{
   return _math_err;
}

void setmatherrflag(bool flag)
{
   _math_err = flag;
}


int matherr(struct _exception *exc)
{
   char *type = "";

   // get string representation of error type
   switch (exc->type) {
      case _DOMAIN:  
         type = _domain;
         break;
      case _SING:
         type = _sing;
         break;
      case _OVERFLOW:
         type = _overflow;
         break;
      case _UNDERFLOW:
         type = _underflow;
         break;
      case _TLOSS:
         type = _tloss;
         break;
      case _PLOSS:
         type = _ploss;
         break;
   }

   _math_err = true;

   // Meldung ausgeben
   /*printf("Exception occurred in %s(%lf, %lf) : %lf\n %s\n", 
          exc->name, 
          exc->arg1, 
          exc->arg2, 
          exc->retval,
          type
         );*/

   return 0;
}