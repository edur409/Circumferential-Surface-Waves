# "On the normal modes of free vibration of inhomogeneous and anisotropic elastic objects"

The first paper on what was later developed into the Resonant Ultrasound Spectroscopy (RUS) technique.  The notation in that paper carried on pretty much forward into the later works of Migliori, Zadler, and more.

You can run the example for a Sphere of copper here: <a target="_blank" href="https://colab.research.google.com/github/edur409/Circumferential-Surface-Waves/blob/main/RUS/RUS_Sphere.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

## Basic equations:

Start with the Lagrangian:

$$ L = \int_{V}(KE - PE)dV $$

> where $KE = \frac{1}{2}\rho \omega^2 u_{i}^2$ is the kinetic energy,

> and $PE = \frac{1}{2} C_{ijkl} u_{i,j} u_{k,l}$

Using basis functions of the form:

$$\Phi_{\lambda} = x^l y^m z^n$$ ,

the Lagrangian becomes:

$$ L = \frac{1}{2} \omega^2 a^T E a - \frac{1}{2} a^T \Gamma a $$

> where $E$ and $\Gamma$ are matrices of order $R$ given by the truncation condition:

$$ l + m + n \le N$$

> with:

$$ R = \frac{3(N + 1)(N + 2)(N + 3)}{6} $$

## E and $\Gamma$ matrices:

$$ E_{\lambda i \lambda' i'} = \delta_{ii'} \int_{V} \Phi_{\lambda, j} \rho \Phi_{\lambda', j'}$$

$$ \Gamma_{\lambda i \lambda' i'} = C_{ijkl} \int_{V} \Phi_{\lambda, j} \Phi_{\lambda', j'}$$

Taking the derivative of the Lagrangian with respect to each of the R components in $a_{i \lambda}$ and setting them to zero, we get the following eigenvalue equation: 

$$ \omega^2 E a = \Gamma a $$

## Python version of the Vissher et al. (1991) Fortran code

The following cell is a translation of the code.  Notice two things, or two typos: 1) RHO is declared but never used in creating the matrix $E$: just used at the end to scale the eigenvalues $W$. 2) When stating the counter $J$, $IG$ was used, instead of $IC$, as an array rather than a constant. 

The output cannot be directly compared with any result in the paper, however, the eigenvalue vector at the end has the first 6 entries made up of zeros, as stated by the authors.

The shape used in this calculation within function $F$ is that of a corner prism.  The stiffness matrix $C_{ijkl}$ is assumed isotropic with Lame parameters equal to 1.  Therefore, the entries $C_{11}$ and $C_{44}$ in Voigt's notation for the elastic matrix are $\lambda + 2\mu = 3.0$ and $\mu = 1.0$.

``` fortran
! The program is named for the basis functions it uses.
      PROGRAM XYZ
! The arrays are dimensioned 252 here, sufficienf to N = 6.
      DIMENSION GAMMA(252,252),E(252,252),W(252),FV(252),FW(252),
     &C(3,3,3,3),LB(252),MB(252),NB(252),IC(252)
! The data following is the elastic tensor Cijkl for our standard
! isotropic material. Any homogeneous anisotropy can be
! described by simply changing these data to include up to 81
! different elastic constants, for a general substance in the
! presence of magnetic fields.
      DATA C/3.,3*0.,1.,3*0.,1.,0.,1.,0.,1.,7*0.,1.,3*0.,1.,3*0.,1.,
     &0.,1.,5*0.,1.,3*0.,3.,3*0.,1.,5*0.,1.,0.,1.,3*0.,1.,3*0.,1.,
     &7*0.,1.,0.,1.,0.,1.,3*0.,1.,3*0.,3 /
      DATA RHO/1./
      TWOPI=2.*ACOS(-1.)
      PRINT*,"PLEASE INPUT NN"
      READ*,NN
! The next 16 lines assign an index IG to each basis funclion (6).
      IG=0
      DO I=1,3
      DO L=1,NN+1
      DO M=1,NN+1
      DO N=1,NN+1
      IF(L+M+N.GT.NN+3) GO TO 2
      IG=IG+1
      IC(IG)=I
      LB(IG)=L-1
      MB(IG)=M-1
      NB(IG)=N-1
    2 END DO
      END DO
      END DO
      END DO
      NR=IG
! In the next 13 statements the elements of the E and Gamma
! matrices are computed.
      DO IG=1,NR
      DO JG=IG,NR
      I=IC(IG)
      J=IC(JG) ! Typo in paper, they have IG(JG), but IG is not a vector!
      LS=LB(IG)+LB(JG)
      MS=MB(IG)+MB(JG)
      NS=NB(IG)+NB(JG)
      GAMMA(IG,JG)=
     & C(I,1,J,1)*FLOAT(LB(IG)*LB(JG))*F(LS-2,MS,NS)
     & +C(I,2,J,2)*FLOAT(MB(IG)*MB(JG))*F(LS,MS-2,NS)
     & +C(I,3,J,3)*FLOAT(NB(IG)*NB(JG))*F(LS,MS,NS-2)
     & +(C(I,1,J,2)*FLOAT(LB(IG)*MB(JG))+C(I,2,J,1)*
     & FLOAT(MB(IG)*LB(JG)))*F(LS-1,MS-1,NS)
     & +(C(I,1,J,3)*FLOAT(LB(IG)*NB(JG))+C(I,3,J,1)*
     & FLOAT(NB(IG)*LB(JG)))*F(LS-1,MS,NS-1)
     & +(C(I,2,J,3)*FLOAT(MB(IG)*NB(JG))+C(I,3,J,2)*
     & FLOAT(NB(IG)*MB(JG)))*F(LS,MS-1,NS-1)
      GAMMA(JG,IG)=GAMMA(IG,JG)
      IF(I.EQ.J) E(IG,JG)=F(LS,MS,NS)   
      E(JG,IG)=E(IG,JG)  
      END DO
      END DO            
! The next line solves the eigenvalue problem (12) using the
! EISPACK13,14 subroutine RSG.
      CALL RSG(252,NR,GAMMA,E,W,0,Z,FV,FW,IERR)
      PRINT*,IERR
! Now obtain the frequencies from the eigenvalues If Cijkl is in
! 10^12dynes/cm^2, \rho is in g/cm^3, and dimensions are in cm, then
! frequencies W are in MHz
      DO I=1,NR
      W(I)=SQRT(AMAX1(0.,W(I))/RHO)/TWOPI
      END DO
! The lowest 36 frequencies are printed out: the first 6 of
! these are always zero if GAMMA is positive.
      PRINT*,"FREQUENCIES FOR CORNER PRISM, NN= ",NN
      PRINT 101,(W(I),I=1,36)
  101 FORMAT(6G12.5)
      END
! Next is a function subprogram for f(p,q,r) for the corner prism
!(Eq. 15). It is straightforward here to substitute function
! subprograms for any of the other objects we have
! considered.
      FUNCTION F(IP,IQ,IR)
      DATA A,B,C/3*1./
      F=A**(IP+1)*B**(IQ+1)*C**(IR+1)*
     &FACT(IP)*FACT(IQ)*FACT(IR)/FACT(IP+IQ+IR+3)
      RETURN
      END
! Factorial subprogram follows.
      FUNCTION FACT(N)
      FACT=1.
      IF(N.LT.2) RETURN
      DO I=2,N
      FACT=FACT*FLOAT(I) ! 1 value for closing the DO LOOP added!
      END DO
      RETURN
      END     
```
