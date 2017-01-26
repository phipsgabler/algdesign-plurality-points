c---------------------------------------------------------------------------

        SUBROUTINE fdepth(U,V,N,X,Y,ALPHA,F,SDEP,HDEP)

c        Calcule la profondeur de Liu et de Tukey du point (u,v) 
c        dans l'echantillon (x,y). Pour plus d'informations, consultez
c        l'article disponible au http://win-www.uia.ac.be/u/statis/ :
c        Rousseuw, P.J., and Ruts, I. (1996), AS 307 : Bivariate location
c        depth, Applied Statistics (JRRSS-C), vol.45, 516-526

      REAL*8 U,V,X(n),Y(n),ALPHA(n)
      REAL*8 P,P2,EPS,D,XU,YU,ANGLE,ALPHK,BETAK,SDEP,HDEP
      INTEGER F(N),GI
      integer n,nums,numh,nt,i,nn,nu,ja,jb,nn2,nbad,nf,j,ki,k
      NUMS=0
      NUMH=0
      SDEP=0.0
      HDEP=0.0
      IF (N.LT.1) RETURN
      P=ACOS(-1.0)
      P2=P*2.0
      EPS=0.00000001
      NT=0
C
C  Construct the array ALPHA.
C
      DO 10 I=1,N
          D=SQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPS) THEN
              NT=NT+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (ABS(XU).GT.ABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      ALPHA(I-NT)=ASIN(YU)
                      IF(ALPHA(I-NT).LT.0.0) THEN
                          ALPHA(I-NT)=P2+ALPHA(I-NT)
                      ENDIF
                  ELSE
                      ALPHA(I-NT)=P-ASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      ALPHA(I-NT)=ACOS(XU)
                  ELSE
                      ALPHA(I-NT)=P2-ACOS(XU)
                  ENDIF
              ENDIF
              IF (ALPHA(I-NT).GE.(P2-EPS)) ALPHA(I-NT)=0.0
          ENDIF
  10  CONTINUE
      NN=N-NT
      IF (NN.LE.1) GOTO 60
C
C  Sort the array ALPHA.
C
      CALL SORT(ALPHA,NN)
C
C  Check whether theta=(U,V) lies outside the data cloud.
C
      ANGLE=ALPHA(1)-ALPHA(NN)+P2
      DO 20 I=2,NN
          ANGLE=MAX(ANGLE,(ALPHA(I)-ALPHA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPS)) GOTO 60
C
C  Make smallest alpha equal to zero,
C  and compute NU = number of alpha < pi.
C
      ANGLE=ALPHA(1)
      NU=0
      DO 30 I=1,NN
          ALPHA(I)=ALPHA(I)-ANGLE
          IF (ALPHA(I).LT.(P-EPS)) NU=NU+1
  30  CONTINUE
      IF (NU.GE.NN) GOTO 60
C
C  Mergesort the alpha with their antipodal angles beta,
C  and at the same time update I, F(I), and NBAD.
C
      JA=1
      JB=1
      ALPHK=ALPHA(1)
      BETAK=ALPHA(NU+1)-P
      NN2=NN*2
      NBAD=0
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPS).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=ALPHA(JA)
              ELSE
                  ALPHK=P2+1.0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
              NBAD=NBAD+K((NF-I),2)
              IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                      BETAK=ALPHA(JB+NU)-P
                  ELSE
                      BETAK=ALPHA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.0
              ENDIF
          ENDIF
  40  CONTINUE
      NUMS=K(NN,3)-NBAD
C
C  Computation of NUMH for halfspace depth.
C
      GI=0
      JA=1
      ANGLE=ALPHA(1)
      NUMH=MIN(F(1),(NN-F(1)))
      DO 50 I=2,NN
          IF(ALPHA(I).LE.(ANGLE+EPS)) THEN
              JA=JA+1
          ELSE
              GI=GI+JA
              JA=1
              ANGLE=ALPHA(I)
          ENDIF
          KI=F(I)-GI
          NUMH=MIN(NUMH,MIN(KI,(NN-KI)))
   50 CONTINUE
C
C  Adjust for the number NT of data points equal to theta:
C
  60  NUMS=NUMS+K(NT,1)*K(NN,2)+K(NT,2)*K(NN,1)+K(NT,3)
      IF (N.GE.3) SDEP=(NUMS+0.0)/(K(N,3)+0.0)
      NUMH=NUMH+NT
      HDEP=(NUMH+0.0)/(N+0.0)
      RETURN
      END
c------------------------------------------------------------------


      INTEGER FUNCTION K(M,J)
      integer m,j
      IF (M.LT.J) THEN
          K=0
      ELSE
          IF (J.EQ.1) K=M
          IF (J.EQ.2) K=(M*(M-1))/2
          IF (J.EQ.3) K=(M*(M-1)*(M-2))/6
      ENDIF
      RETURN
      END


c------------------------------------------------------------------

      SUBROUTINE SORT(B,N)
C  Sorts an array B (of length N<=1000) in O(NlogN) time.

      REAL*8 B(N),x(n)
      integer q(n),i,n
      
      call indexx(n,b,q)
      
      do 10 i=1,n
        x(i)=b(i)
10    continue
 
      do 20 i=1,n
        b(i)=x(q(i))
20    continue
      end

