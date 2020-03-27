/****************************************************************************
* Name: jlp_orbit1_gauss.cpp
* (for visual orbits)
* 
* JLP
* Version 08/06/2015
****************************************************************************/
#include "jlp_orbit1.h"

/*****************************************************************************
* From A. Tokovinin 
* LSQ1_GaussInversion with the pivot method
*****************************************************************************/
/*****************************************************************************
* Initialize private variables 
*****************************************************************************/
int JLP_Orbit1::LSQ1_GaussInversion(double alpha[NELEMENTS][NELEMENTS],
                                    const int n_free_oelements, 
                                    double beta[NELEMENTS])
{
double big, tmp, piv_inv;
int *ipiv, *indx_row, *indx_col;
int nmax = NELEMENTS, irow, icol; 
int i, j, k, l, ll, nn;

// Allocate memory space for row, column indices and pivots:
indx_row = new int[nmax];
indx_col = new int[nmax];
ipiv = new int[nmax];

nn = n_free_oelements;

// Initialize pivots:
for(i = 0; i < nn; i++) {
   ipiv[i] = 0;
  }

for(i = 0; i < nn; i++) {
  big = 0.;
  for(j = 0; j < nn; j++) {
   if(ipiv[j] != 1) {
      for(k = 0; k < nn; k++) {
           if(ipiv[k] == 0) {
              if(ABS(alpha[j][k]) >= big) {
                big = alpha[j][k];
                irow = j;
                icol = k;
                }
           } else if(ipiv[k] > 1) {
             fprintf(stderr, "LSQ1_GaussInversion/Fatal error: singular matrix\n");
             fprintf(stderr, " piv_%d > 1\n", k);
             return(-1);
           }
        } // EOF loop on k
     } // EOF ipiv[j] != 1
  } // EOF loop on j
   ipiv[icol] += 1;
   if(irow != icol) {
// Swap irow/icol:
        for(l = 0; l < nn; l++) {
          tmp = alpha[irow][l];
          alpha[irow][l] = alpha[icol][l];
          alpha[icol][l] = tmp;
          } // EOF loop on l
        tmp = beta[irow];
        beta[irow] = beta[icol];
        beta[icol] = tmp;
      } // EOF irow != icol
    indx_row[i] = irow;
    indx_col[i] = icol;
    if(alpha[icol][icol] == 0.) {
      fprintf(stderr, "LSQ1_GaussInversion/Fatal error: singular matrix\n");
      fprintf(stderr, " alpha(%d,%d) null\n", icol, icol);
      return(-1);
      }
     piv_inv = 1. / alpha[icol][icol];
     alpha[icol][icol] = 1.;
     for(l = 0; l < nn; l++) {
       alpha[icol][l] *= piv_inv; 
       }
     beta[icol] *= piv_inv;
     for(ll = 0; ll < nn; ll++) {
        if(ll != icol) {
          tmp = alpha[ll][icol];
          alpha[ll][icol] = 0.;
          for(l = 0; l < nn; l++) {
            alpha[ll][l] -= alpha[icol][l] * tmp; 
            }
          beta[ll] -= beta[icol] * tmp; 
          } // EOF ll != icol
       } // EOF loop on ll
} // EOF loop on i

for(l = nn - 1; l >= 0; l--) {
   if(indx_row[l] != indx_col[l]) {
     for(k = 0; k < nn; k++) {
        tmp = alpha[k][indx_row[l]];
        alpha[k][indx_row[l]] = alpha[k][indx_col[l]];
        alpha[k][indx_col[l]] = tmp;
     } // EOF loop on k
   } // EOF if(indx_row[l]... 
 } // EOF loop on l

delete[] indx_row;
delete[] indx_col;
delete[] ipiv;

return(0);
}
//        CALL GAUSSJ(ALPHA,NOEL2,10,BETA,1,1)
#ifdef TTT
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                WRITE(6,*)'GAUSSJ/Singular matrix (piv_',K,' lt 1)'
                STOP
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF(A(ICOL,ICOL).EQ.0.)THEN 
          WRITE(6,*)'GAUSSJ/Singular matrix: A_icol(',I,') null'
          STOP
        ENDIF
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
#endif
