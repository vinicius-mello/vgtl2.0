#ifndef LBFGSB01_H
#define LBFGSB01_H

extern "C" {

	void setb01_(int * n, int * m, double * x, double * l, 
							double * u, double * f, double * g,
							double * factr, double * pgtol, double * wa,
							int * iwa, char * task, int * iprint, char * csave,
							int * lsave, int * isave, double * dsave);
				
}

/*
      subroutine setb01(n, m, x, l, u, f, g, factr, pgtol, wa, iwa,
     +                 task, iprint, csave, lsave, isave, dsave)
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, 
     +                 iwa(3*n), isave(44)
      double precision f, factr, pgtol, x(n), l, u, g(n),
     +                 wa(2*m*n+4*n+11*m*m+8*m), dsave(29)
*/		 

#endif
