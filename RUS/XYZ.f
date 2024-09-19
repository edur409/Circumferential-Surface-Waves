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
