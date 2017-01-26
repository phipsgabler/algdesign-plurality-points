c ====================================================================
c                        moy_tr.f
c ====================================================================

       subroutine getcol(p,n,x,m,y)

c        Extrait la colonne m d'une matrice n*p
c        p   : entier egal au nombre de dimensions
c        n   : entier egal au nombre de points
c        x   : Matrice n*p
c        m   : entier = indice de la colonne a extraire 
c        y   : vecteur de p elements = colonne desiree
c        Jean-Francois Plante ete 1999


       integer p,n,i,m
       real*8 x(n,p),y(n)

       do 10 i=1,n
           y(i)=x(i,m)
10     continue
       end

c----------------------------------------------------------------------


       subroutine liumed(x,y,n,maxdpth,med)

c        Calcule la mediane de Liu = moyenne des points de l'echantilon
c        dont la profondeur est maximale.
c        x   : Coordonnees en x des points
c        y   : Coordonnees en y des points
c        n   : entier egal au nombre de points
c        maxdpth: Profondeur de l'approximation de la mediane trouvee 
c        med : vecteur de 2 elements = approximation de la mediane
c        Jean-Francois Plante ete 1999

       integer n,i,k,cont,q(n)
       integer f(n),dpth(n),maxdpth
       real*8 x(n),y(n),med(2)
       real*8 sdep,hdep,sdp(n)
       real*8 alpha(n),u,v

       do 5 i=1,n
              u=x(i)
              v=y(i)
              call fdepth(u,v,n,x,y,alpha,f,sdep,hdep)
              sdp(i)=sdep
              dpth(i)=int(sdep*(K(N,3)+0.0)+.5)
5      continue

       call indexx(n,sdp,q)
       maxdpth=dpth(q(n))

       med(1)=0.
       med(2)=0.
       cont=0
       do 60 i=1,n
             if(dpth(i).eq.maxdpth) then
                   med(1)=med(1)+x(i)
                   med(2)=med(2)+y(i)
                   cont=cont+1
             endif
60     continue

       med(1)=med(1)/dble(cont)
       med(2)=med(2)/dble(cont)

       end



c----------------------------------------------------------------------


      subroutine iso3d(x,y,z,n,t,liu,xx,yy)

c        Calcule la profondeur de chaque point (xx(i),yy(j)) dans
c        l'echantillon (x(k),y(k)).
c        n   : entier egal au nombre de points-1 de l'echantillon
c        x   : Coordonnees en x des points de l'echantillon
c        y   : Coordonnees en y des points de la d'echantillon
c        t   : entier egal au nombre de points-1 de la discretisation
c        xx  : Coordonnees en x des points de la discretisation
c        yy  : Coordonnees en y des points de la discretisation
c        z   : z[i,j] = profondeur de (x(i),y(j))
c        liu : 0 => profondeur de Tukey 1=> prof. de Liu utilisee
c        Jean-Francois Plante ete 1999

      integer n,t,i,j,f(n),liu
      real*8 x(n),y(n),z(t+1,t+1),xx(t+1),yy(t+1)
      real*8 sdep,hdep,alpha(n)

      do 10 i=1,t+1
          do 5 j=1,t+1
             call fdepth(xx(i),yy(j),n,x,y,alpha,f,sdep,hdep)
             if (liu.ne.0) then
                 z(i,j)=sdep
             else
                 z(i,j)=hdep
             endif
5         continue
10    continue
      end


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
c------------------------------------------------------------------
      SUBROUTINE SORT2(B,I1,I2,n)
C
C  Sorts a integer array B of length N and permutes two real arrays 
C  I1 and I2 and one real array R in the same way. (descending order)
C
      INTEGER N,b(n),q(n),x(n),i
      REAL*8 I1(N),I2(n),II1(n),II2(n),bb(n)

      do 5 i=1,n
        bb(i)=dble(b(i))
5     continue

      call indexx(n,bb,q)
      
      do 10 i=1,n
        x(i)=b(i)
        II1(i)=I1(i)
        II2(i)=I2(i)
10    continue
 
      do 20 i=1,n
        b(i)=x(n+1-q(i))
        I1(i)=II1(n+1-q(i))
        I2(i)=II2(n+1-q(i))
20    continue

      END

c ====================================================================
c                        depth3d.f
c ====================================================================

c------------------------------------------------------------------

      SUBROUTINE STAND(N,X,Y,Z,U,V,W,XN,EPS,err)

c        voir depth3

      integer n,err
      real*8 X(N),Y(N),Z(N),U,V,W,XN(N),EPS


      CALL STAND1(N,X,U,XN,EPS,1,err)
      CALL STAND1(N,Y,V,XN,EPS,2,err)
      CALL STAND1(N,Z,W,XN,EPS,3,err)

      RETURN
      END

c------------------------------------------------------------------

      SUBROUTINE STAND1(N,X,U,XN,EPS,J,err)

c        voir depth3

      integer i,j,n,err,jn
      real*8 X(N),U,XN(N),EPS,findq
      real*8 QLOC,QSCA,AVE,VAR
      

      jn=0
      DO 10 I=1,N
         XN(I)=X(I)
 10   CONTINUE         

      IF ((2*INT(N/2)).EQ.N) THEN
         QLOC=FINDQ(XN,N,N/2)
         QLOC=(FINDQ(XN,N,(N/2)+1)+QLOC)/2.D0
      ELSE
         QLOC=FINDQ(XN,N,INT(N/2)+1)
      ENDIF
      DO 30 I=1,N
         XN(I)=DABS(X(I)-QLOC)
 30   CONTINUE
      IF ((2*INT(N/2)).EQ.N) THEN
         QSCA=FINDQ(XN,N,N/2)
         QSCA=(FINDQ(XN,N,(N/2)+1)+QSCA)/2.D0
      ELSE
         QSCA=FINDQ(XN,N,INT(N/2)+1)
      ENDIF
      IF (DABS(QSCA).LT.EPS) THEN
         AVE=0.D0
         DO 40 I=1,N
            AVE=AVE+X(I)
 40      CONTINUE
         AVE=AVE/(N+0.D0)
         VAR=0.D0
         DO 50 I=1,N
            VAR=VAR+(X(I)-AVE)*(X(I)-AVE)
 50      CONTINUE  
         IF (N.NE.1) VAR=VAR/(N-1.D0)
         IF (DABS(VAR).LT.EPS) THEN
            err=J+10
            QSCA=1.D0
         ELSE
            err=j
            QSCA=DSQRT(VAR)
         ENDIF
      ENDIF
      JN=JN+1
      DO 60 I=1,N
         X(I)=(X(I)-QLOC)/QSCA
 60   CONTINUE         
      U=(U-QLOC)/QSCA

      RETURN
      END

c------------------------------------------------------------------


      SUBROUTINE DEPTH3(N,U,V,W,X,Y,Z,ALPHA,F,XN,YN,EPS,NDIM,NDEP)

C
C  This program computes the halfspace depth of
C  a point (U,V,W) in a 3-dimensional data set (X,Y,Z) of size N.
C
C  The data set is read from a file in free format, measurements for
C  each new object should start on a new line.
C
C  N = the actual size of the data set.
C  X(MAXN) = measurements for the first variable.
C  Y(MAXN) = measurements for the second variable.
C  Z(MAXN) = measurements for the third variable.
C  U = first coordinate of the point of which the depth is computed.
C  V = second coordinate of the point of which the depth is computed.
C  W = third coordinate of the point of which the depth is computed.
C  NDEP = the halfspace depth.
C

c        Pour plus de renseignements, consultez l'article suivant,]
c        disponible au http://win-www.uia.ac.be/u/statis/ :
c        Rousseeuw, P.J. and Struyf, A. (1998), Computing location 
c        depth and regression depth in higher dimensions, Statistics
c        and Computing, vol.8, 193-203.  



      integer n,ndim,ndep,i,j,ntnul,ntpos,ntneg,nh
      real*8 U,V,W,X(N),Y(N),Z(N),ALPHA(N),XN(N),YN(N)
      real*8 A(2,3),B(3,2),EPS,DP
      INTEGER F(N)

C
C Make theta the center of the dataset.
C
      DO 10 I=1,N
         X(I)=X(I)-U
         Y(I)=Y(I)-V
         Z(I)=Z(I)-W
 10   CONTINUE
      NDIM=3
C
C Handle special cases where N is less or equal to 1.
C
      IF (N.LE.1) THEN
         IF ((N.EQ.1) .AND. (DABS(X(1)).LE.EPS) .AND. 
     +        (DABS(Y(1)).LE.EPS) .AND. 
     +        (DABS(Z(1)).LE.EPS)) THEN
            NDEP=1
         ELSE
            NDEP=0
         ENDIF
         RETURN
      ENDIF
C
C General case: initialize halfspace depth.
C
      NDEP=N
C
C Loop over all lines (theta,x(i)).
C
      DO 20 I=1,N
         IF ((DABS(X(I)).LE.EPS).AND.(DABS(Y(I)).LE.EPS).AND.
     +        (DABS(Z(I)).LE.EPS)) THEN
            GOTO 20
         ENDIF
C
C Calculate the matrix of the orthogonal projection on the plane through 
C theta, orthogonal to the line through theta and x(i).
C Let the third coordinate coincide with the line (theta,x(i)).
C
         IF (DABS(X(I)).GT.EPS) THEN
            B(2,1)=1.D0
            B(3,1)=1.D0
            B(1,1)=-(Y(I)+Z(I))/X(I)
         ELSEIF (DABS(Y(I)).GT.EPS) THEN
            B(1,1)=1.D0
            B(3,1)=1.D0
            B(2,1)=-(X(I)+Z(I))/Y(I)
         ELSE
            B(1,1)=1.D0
            B(2,1)=1.D0
            B(3,1)=-(X(I)+Y(I))/Z(I)
         ENDIF
         B(1,2)=B(2,1)*Z(I)-B(3,1)*Y(I)
         B(2,2)=B(3,1)*X(I)-B(1,1)*Z(I)
         B(3,2)=B(1,1)*Y(I)-X(I)*B(2,1)
      
         A(1,1)=(B(2,2)*Z(I)-Y(I)*B(3,2))
         A(1,2)=-(B(1,2)*Z(I)-X(I)*B(3,2))
         A(1,3)=(B(1,2)*Y(I)-B(2,2)*X(I))
         A(2,1)=-(B(2,1)*Z(I)-Y(I)*B(3,1))
         A(2,2)=(B(1,1)*Z(I)-X(I)*B(3,1))
         A(2,3)=-(B(1,1)*Y(I)-X(I)*B(2,1))
C
C Compute the new planar coordinates for all points.
C If a point collapses with theta, identify its position:
C     NTNUL = real ties, 
C     NTPOS = the original point lies on the positive side of 
C             the projection plane,
C     NTNEG = the original point lies on the negative side of 
C             the projection plane.
C
         NTNUL=0
         NTPOS=0
         NTNEG=0
         DO 30 J=1,N
            XN(J)=X(J)*A(1,1)+Y(J)*A(1,2)+Z(J)*A(1,3)
            YN(J)=X(J)*A(2,1)+Y(J)*A(2,2)+Z(J)*A(2,3)
            IF ((DABS(XN(J)).LE.EPS).AND.
     +           (DABS(YN(J)).LE.EPS)) THEN
               DP=X(J)*X(I)+Y(J)*Y(I)+Z(J)*Z(I)
               IF (DABS(DP).LE.EPS) THEN
                  NTNUL=NTNUL+1
               ELSEIF (DP.GT.EPS) THEN 
                  NTPOS=NTPOS+1
               ELSE 
                  NTNEG=NTNEG+1
               ENDIF
            ENDIF
 30      CONTINUE
         IF ((NTNUL+NTNEG+NTPOS).EQ.N) GOTO 50
C
C Compute the halfspace depth in two dimensions.
C
         CALL DEPTH2(0.D0,0.D0,N,XN,YN,ALPHA,F,NH,NTPOS,NTNEG,
     +        NTNUL,EPS,NDIM)
C
C Update the three-dimensional halfspace depth.
C
         NDEP=MIN0(NDEP,NH)
 20   CONTINUE
      RETURN
C
C All points and theta lie on one line.
C
 50   NDEP=MIN0(NTNUL+NTPOS,NTNUL+NTNEG)
      NDIM=1

      RETURN
      END

c------------------------------------------------------------------

      SUBROUTINE DEPTH2(U,V,N,X,Y,ALPHA,F,NH,NTPOS,NTNEG,NTNUL,EPS,NDIM)

      integer n,ntpos,ntneg,ntnul,ndim,numh,nt,nd,i,nn,nu,ja,jb,nn2,nf
      integer j,ki
      real*8 U,V,X(N),Y(N),ALPHA(N)
      real*8 P,P2,EPS,D,XU,YU,ANGLE,ALPHK,BETAK
      INTEGER F(N),GI,NH
      NUMH=0
      NH=0
      IF (N.LT.1) RETURN
      P=DACOS(-1.D0)
      P2=P*2.D0
      NT=0
      ND=0
C
C  Construct the array ALPHA.
C
      DO 10 I=1,N
          D=DSQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPS) THEN
              NT=NT+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (DABS(XU).GT.DABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      ALPHA(I-NT)=DASIN(YU)
                      IF(ALPHA(I-NT).LT.0.D0) THEN
                          ALPHA(I-NT)=P2+ALPHA(I-NT)
                      ENDIF
                  ELSE
                      ALPHA(I-NT)=P-DASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      ALPHA(I-NT)=DACOS(XU)
                  ELSE
                      ALPHA(I-NT)=P2-DACOS(XU)
                  ENDIF
              ENDIF
              IF (ALPHA(I-NT).GE.(P2-EPS)) ALPHA(I-NT)=0.D0
          ENDIF
  10  CONTINUE
      NN=N-NT
      IF (NN.LE.1) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C
C  Sort the array ALPHA.
C
      CALL SORT(ALPHA,NN)
C
C  Check whether theta=(U,V) lies outside the data cloud.
C
      ANGLE=ALPHA(1)-ALPHA(NN)+P2
      DO 20 I=2,NN
          ANGLE=DMAX1(ANGLE,(ALPHA(I)-ALPHA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPS)) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C
C  Make smallest alpha equal to zero,
C  and compute NU = number of alpha < pi.
C
      ANGLE=ALPHA(1)
      NU=0
      DO 30 I=1,NN
          ALPHA(I)=ALPHA(I)-ANGLE
          IF (ALPHA(I).LT.(P-EPS)) NU=NU+1
          IF ((DABS(ALPHA(I)).LE.EPS).OR.
     +         (DABS(ALPHA(I)-P).LE.EPS)) ND=ND+1
  30  CONTINUE
      IF (ND.EQ.NN) NDIM=2
      IF (NU.GE.NN) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C
C  Mergesort the alpha with their antipodal angles beta,
C  and at the same time update I, and F(I).
C
      JA=1
      JB=1
      ALPHK=ALPHA(1)
      BETAK=ALPHA(NU+1)-P
      NN2=NN*2
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPS).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=ALPHA(JA)
              ELSE
                  ALPHK=P2+1.D0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
              IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                      BETAK=ALPHA(JB+NU)-P
                  ELSE
                      BETAK=ALPHA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.D0
              ENDIF
          ENDIF
  40  CONTINUE
C
C  Compute the halfspace depth.
C  Correct for ties (which where no ties in three dimensions)
C  by considering some small deviations of the planes (without changing
C  their intersection line with the plane which contains the 2D-data).
C
      GI=0
      JA=1
      ANGLE=ALPHA(1)
      NUMH=MIN0(MIN0(F(1)+NTNEG,F(1)+NTPOS),
     +     MIN0(NN-F(1)+NTNEG,NN-F(1)+NTPOS))
      DO 50 I=2,NN
         IF(ALPHA(I).LE.(ANGLE+EPS)) THEN
            JA=JA+1
         ELSE
            GI=GI+JA
            JA=1
            ANGLE=ALPHA(I)
         ENDIF
         KI=F(I)-GI
         NUMH=MIN0(NUMH,MIN0(MIN0(KI+NTNEG,KI+NTPOS),
     +        MIN0(NN-KI+NTNEG,NN-KI+NTPOS)))
 50   CONTINUE
C
C  Adjust for the number NTNUL of data points equal to theta.
C
 60   NH=NUMH+NTNUL
      RETURN
      END

c ====================================================================
c                        moytrpd.f
c ====================================================================

      SUBROUTINE STANDPD(MAXN,MAXP,N,NP,X,T,EPS,err,NDEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer maxn,maxp,n,np,ndep,i,j,jn,err(np)
      real*8 X(MAXN,MAXP),T(NP),XN(N),EPS
      real*8 QLOC,QSCA,AVE,VAR

      JN=0
      DO 10 J=1,NP
         DO 20 I=1,N
            XN(I)=X(I,J)
 20      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QLOC=FINDQ(XN,N,N/2)
            QLOC=(FINDQ(XN,N,(N/2)+1)+QLOC)/2.D0
         ELSE
            QLOC=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         DO 30 I=1,N
            XN(I)=DABS(X(I,J)-QLOC)
 30      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QSCA=FINDQ(XN,N,N/2)
            QSCA=(FINDQ(XN,N,(N/2)+1)+QSCA)/2.D0
         ELSE
            QSCA=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         IF (DABS(QSCA).LT.EPS) THEN
            AVE=0.D0
            DO 40 I=1,N
               AVE=AVE+X(I,J)
 40         CONTINUE
            AVE=AVE/(N+0.D0)
            VAR=0.D0
            DO 50 I=1,N
               VAR=VAR+(X(I,J)-AVE)*(X(I,J)-AVE)
 50         CONTINUE  
            IF (N.NE.1) VAR=VAR/(N-1.D0)
            IF (DABS(VAR).LT.EPS) THEN
               IF (DABS(T(J)-X(1,J)).GT.EPS) NDEP=0
               NP=NP-1
               err(j)=-1
               GOTO 10
            ELSE
               err(j)=-2
               QSCA=DSQRT(VAR)
            ENDIF
         ENDIF
         JN=JN+1
         DO 60 I=1,N
            X(I,JN)=(X(I,J)-QLOC)/QSCA
 60      CONTINUE         
         T(JN)=(T(J)-QLOC)/QSCA
 10   CONTINUE

      RETURN
      END

c--------------------------------------------------------------------

      SUBROUTINE HDEPTH(N,P,NNP,NDIR,maxn,maxp,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,err,err2,NDEP,NSIN)
C
C  This program computes an approximation for the halfspace depth of
C  a point T in an NP-dimensional data set X of size N.
C  Missing values are not allowed.
C 
C  N = the size of the data set.
C  P = the number of variables in the data set.
C  NNP = the number of variables currently being used.
C  NDIR = the number of samples to draw.
C  X(N,P) = the full data set.
C  T(P) = point of which the depth is computed.
C  NDEP = final approximation for the halfspace depth.
C  JSAMP(P) = indices of points in the current sample.
C

c        Pour plus de renseignements, consultez l'article suivant,]
c        disponible au http://win-www.uia.ac.be/u/statis/ :
c        Rousseeuw, P.J. and Struyf, A. (1998), Computing location 
c        depth and regression depth in higher dimensions, Statistics
c        and Computing, vol.8, 193-203.  


      integer n,p,np,nnp,ndir,err,ndep,nsin,j,numh,nt,l
      integer err2(ndir),nnp1,maxn,maxp

      real*8 X(maxN,maxP),T(maxP),R(maxP),EPS,AVE(maxP)
      real*8 EVECS(maxP,maxP),EVALS(maxP),COV(maxP,maxP)
      INTEGER JSAMP(maxP),ierr
C
C  Initialize the number of singular samples.
C
      np=p
      err=0
      NSIN=0
C
C  Handle special case where N is equal to 1.
C
      IF (N.EQ.1)THEN
         IF (N.EQ.1) THEN
            DO 10 J=1,NP
               IF (DABS(X(1,J)-T(J)).GT.EPS) GOTO 15
 10         CONTINUE
            NDEP=1
            RETURN
         ENDIF
 15      NDEP=0
         RETURN
      ENDIF
C
C  Handle special case where NNP is equal to 1.
C
 25   IF (NNP.EQ.1) THEN
         NUMH=0
         NT=0
         DO 20 L=1,N
            IF (X(L,1).GT.(T(1)+EPS)) THEN
               NUMH=NUMH+1
            ELSEIF (X(L,1).GE.(T(1)-EPS)) THEN
               NT=NT+1
            ENDIF
 20      CONTINUE
         NDEP=MIN(NUMH+NT,N-NUMH)
         RETURN
      ENDIF
C
C  General case: call subroutine DEP.
C
      CALL DEP(N,NNP,NDIR,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN,err2)
C
C  If all points and theta are identified as lying on the same hyperplane,
C  reduce the dimension of the data set by projection on that hyperplane,
C  and compute the depth on the reduced data set.
C
      IF (NSIN.EQ.(-1)) THEN
         NSIN=0
         err=-1
         NNP1=NNP
         NNP=NNP-1
         CALL REDUCE(N,NNP,NNP1,MAXN,MAXP,X,T,R,EVECS,JSAMP,IERR)
         IF (IERR.LT.0) THEN
            err=-2
            GOTO 50
         ENDIF
         
         GOTO 25
      ENDIF

 50   RETURN

      END

c--------------------------------------------------------------------


      SUBROUTINE DEP(N,P,NDIR,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN,err)

      integer n,np,ndir,ndep,nsin,err(ndir),p,nrun,nran
      integer i,nsamp,j
      real*8 X(N,P),T(P),R(P)
      real*8 K,KT,EPS,RAN
      real*8 EVECS(P,P),EVALS(P),COV(P,P),AVE(P)
      INTEGER JSAMP(P),ierr,l,nt,numh
C
C  Initialize halfspace depth and random seed.
C
      np=p
      
      NDEP=N
      NRUN=0
      DO 100 NRAN=1,NDIR
C     
C  Draw a random sample of size np.
C     
         CALL RANDM(NRUN,RAN)
         I=N*RAN+1.
         IF(I.GT.N)I=N
         JSAMP(1)=I
         NSAMP=1
 20      CALL RANDM(NRUN,RAN)
         L=N*RAN+1.
         IF(L.GT.N)L=N
         DO 30 J=1,NSAMP
            IF(L.EQ.JSAMP(J)) GOTO 20
 30      CONTINUE
         NSAMP=NSAMP+1
         JSAMP(NSAMP)=L
         IF (NSAMP.LT.NP)GOTO 20
C     
C  Compute the covariance matrix of the sample.
C     
         DO 40 J=1,NP
            AVE(J)=0.D0
            DO 50 I=1,NP
                  AVE(J)=AVE(J)+X(JSAMP(I),J)
 50         CONTINUE
            AVE(J)=AVE(J)/NP
 40      CONTINUE
         DO 60 J=1,NP
            DO 70 L=1,J
               COV(J,L)=0.D0
               DO 80 I=1,NP
                  COV(J,L)=COV(J,L)+(X(JSAMP(I),J)-AVE(J))
     +                 *(X(JSAMP(I),L)-AVE(L))
 80            CONTINUE
               COV(J,L)=COV(J,L)/(NP-1)
               COV(L,J)=COV(J,L)
 70         CONTINUE
 60      CONTINUE
C     
C  Compute the eigenvalues and corresponding eigenvectors 
C  of the covariance matrix.
C     
         CALL EIGEN(NP,NP,COV,EVALS,EVECS,R,IERR)
         IF (IERR.NE.0) THEN
            err(nran)=ierr
            NSIN=NSIN+1
            GOTO 100
         ENDIF
         IF (EVALS(1).GT.EPS) THEN
            err(nran)=-1
            NSIN=NSIN+1
            GOTO 100
         ENDIF
C     
C  Test for singularity of the sample.
C     
         IF (EVALS(2).LE.EPS) THEN
            NSIN=NSIN+1
         ENDIF
C
C  Project all points on the line through theta with direction given by 
C  the eigenvector of the smallest eigenvalue, i.e. the direction 
C  orthogonal on the hyperplane given by the np-subset.
C  Compute the one-dimensional halfspace depth of theta on this line.
C         
         KT=0.D0
         NT=0
         DO 90 J=1,NP
            IF (DABS(EVECS(J,1)).LE.EPS) THEN
               NT=NT+1
            ELSE
               KT=KT+T(J)*EVECS(J,1)
            ENDIF
 90      CONTINUE
         IF (NT.EQ.NP) THEN
            err(nran)=-2
            NSIN=NSIN+1
            GOTO 100
         ENDIF

         NUMH=0
         NT=0
         DO 95 L=1,N
            K=0.D0
            DO 96 J=1,NP
               K=K+EVECS(J,1)*X(L,J)
 96         CONTINUE
            K=K-KT
            IF (K .GT. EPS) THEN
               NUMH=NUMH+1
            ELSEIF (K .GE. (0.D0-EPS)) THEN
               NT=NT+1
            ENDIF
 95      CONTINUE
C
C  If all projections collapse with theta, return to reduce the dimension.
C
         IF (NT.EQ.N) THEN
            NSIN=-1
            RETURN
         ENDIF
C
C  Update the halfspace depth.
C
         NDEP=MIN(NDEP,min(NUMH+NT,N-NUMH))
 100  CONTINUE
      RETURN
      END

c--------------------------------------------------------------------


      SUBROUTINE REDUCE(N,NNP,NNP1,MAXN,MAXP,X,T,R,EVECS,W,IERR)
   
      real*8 X(MAXN,MAXP),T(NNP1),R(NNP1)
      real*8 EVECS(NNP1,NNP1)
      INTEGER W(NNP),n,nnp,nnp1,maxn,maxp,ierr,i,j,io

      IERR=0
C
C  Invert matrix of base vectors EVECS.
C
      CALL VERT(EVECS,NNP+1,NNP+1,W,IERR)
      IF (IERR.LT.0) RETURN
C
C  Compute new NNP-dimensional coordinates for all points and theta.
C      
      DO 30 I=2,NNP+1
         R(I-1)=T(1)*EVECS(I,1)
         DO 31 J=2,NNP+1
            R(I-1)=R(I-1)+T(J)*EVECS(I,J)
 31      CONTINUE
 30   CONTINUE
      DO 32 I=1,NNP
         T(I)=R(I)
 32   CONTINUE

      DO 40 IO=1,N
         DO 41 I=2,NNP+1
            R(I-1)=X(IO,1)*EVECS(I,1)
            DO 42 J=2,NNP+1
               R(I-1)=R(I-1)+X(IO,J)*EVECS(I,J)
 42         CONTINUE
 41      CONTINUE
         DO 43 I=1,NNP
            X(IO,I)=R(I)
 43      CONTINUE
 40   CONTINUE
      RETURN
      END

c--------------------------------------------------------------------


      SUBROUTINE RANDM(NRUN,RAN)
CC   WE PROGRAMMED THIS GENERATOR OURSELVES BECAUSE WE WANTED IT
CC   TO BE MACHINE INDEPENDENT. IT SHOULD RUN ON MOST COMPUTERS 
CC   BECAUSE THE LARGEST INTEGER USED IS LESS THAN 2**30 . THE PERIOD 
CC   IS 2**16=65536, WHICH IS GOOD ENOUGH FOR OUR PURPOSES. 
      DOUBLE PRECISION RAN,RY
      INTEGER*4 NRUN,K
      NRUN=NRUN*5761+999
      K=NRUN/65536
      NRUN=NRUN-K*65536 
      RY=NRUN
      RAN=RY/65536.0
      RETURN
      END 

c--------------------------------------------------------------------


      function findq(aw,ncas,k)
cc  Finds the k-th order statistic of the array aw of length ncas.

      real*8 findq
      real*8 aw(ncas)
      real*8 ax,wa
      integer ncas,k,l,lr,jnc,j

      l=1
      lr=ncas
 20   if(l.ge.lr) goto 90
      ax=aw(k)
      jnc=l
      j=lr
 30   if(jnc.gt.j) goto 80
 40   if(aw(jnc).ge.ax) goto 50
      jnc=jnc+1
      goto 40
 50   if(aw(j).le.ax) goto 60
      j=j-1
      goto 50
 60   if(jnc.gt.j) goto 70
      wa=aw(jnc)
      aw(jnc)=aw(j)
      aw(j)=wa
      jnc=jnc+1
      j=j-1
 70   goto 30
 80   if(j.lt.k) l=jnc
      if(k.lt.jnc) lr=j
      goto 20
 90   findq=aw(k)
      return
      end

c--------------------------------------------------------------------


      SUBROUTINE VERT(V,LV,N,W,IERR)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module VERT from package NAPACK.
* Retrieved from NETLIB on Wed Feb 19 03:31:44 1997.
* ======================================================================
C
C      ________________________________________________________
C     |                                                        |
C     |                INVERT A GENERAL MATRIX                 |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
C     |                      ELEMENTS                          |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      DOUBLE PRECISION V(LV,1),S,T
      INTEGER W(1),I,J,K,L,M,N,P,lv,ierr
      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = DABS(V(L,L))
      DO 20 I = M,N
           T = DABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0. ) GOTO 120
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.
      S = 1./S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0. ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0. ) GOTO 120
      V(1,1) = 1./V(1,1)
      RETURN
120   IERR=-1
      RETURN
      END

c--------------------------------------------------------------------


      subroutine eigen(nm,n,a,w,z,fv1,ierr)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module RS from package EISPACK.
* Retrieved from NETLIB on Wed Nov 27 07:41:24 1996.
* ======================================================================
c
      integer n,nm,ierr
      real*8 a(nm,n),w(n),z(nm,n),fv1(n)

c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tql2.
c           the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
c --------------------------------------------------------------------
      real*8 function pythag(a,b)
      real*8 a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      real*8 p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c -------------------------------------------------------------------
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

c ---------------------------------------------------------------------
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),z(nm,n)
      real*8 f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h

  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end

c ====================================================================
c                        ojadepth.f
c ====================================================================

      subroutine ojaiso3d(x,z,n,t,xx,yy)

c        Calcule la profondeur de chaque point (xx(i),yy(j)) dans
c        l'echantillon (x(k),y(k)).
c        n   : entier egal au nombre de points-1 de l'echantillon
c        x   : Coordonnees en x des points de l'echantillon
c        y   : Coordonnees en y des points de la d'echantillon
c        t   : entier egal au nombre de points-1 de la discretisation
c        xx  : Coordonnees en x des points de la discretisation
c        yy  : Coordonnees en y des points de la discretisation
c        z   : z[i,j] = profondeur de (x(i),y(j))
c        Jean-Francois Plante ete 1999

      integer n,t,i,j
      real*8 x(n,2),z(t+1,t+1),xx(t+1),yy(t+1)
      real*8 odep,u(2)

      do 10 i=1,t+1
          u(1)=xx(i)
          do 5 j=1,t+1
             u(2)=yy(j)
             call ojadepth(x,u,2,n,odep)
             z(i,j)=odep
5         continue
10    continue
      end

c----------------------------------------------------------------------

      recursive subroutine reprow(p,n,k,l,sum,mat,x,u)
c used for the recursive calculation of Oja Depth
c p = #dim; n = sample size; k = column to replace; 
c l = where about in the data (what was the last line added)
c sum = sum of areas so far; mat = matrix with k-1 first columns filled
c x = data set; u = vector with respect to which depth is calculated.
      integer p,n,k,l
      real*8 x(n,p),u(p),one,sum,mat(p,p),col(p),v(p),det
      
      one=1.0
    
      if(k.eq.0) then
        call determinant(mat,p,det)
        sum=sum+dabs(det)
      else        

      do 10 i=l,(n-k+1)
        call getrow(p,n,x,i,col)
        call clv(p,one,u,-one,col,v)
        call putrow(p,p,mat,k,v)
        call reprow(p,n,k-1,i+1,sum,mat,x,u)
10    continue          
      endif
      
      end


c----------------------------------------------------------------------

      subroutine ojadepth(x,u,p,n,dpth)

      integer n,p,i,j
      real*8 x(n,p),u(p),dpth,sum,mat(p,p),aire
      real*8 coef,fac,det,one

c        Calcule la profondeur de u dans x selon la definition de Oja.
c        p   : entier egal au nombre de dimensions
c        n   : entier egal au nombre de points
c        x   : Matrice n*p contenant les n points en p dimensions
c        u   : vecteur de p composantes, coordonnees du point dont on
c              calcule la profondeur
c        dpth: Profondeur de u dans x
c        Jean-Francois Plante ete 1999

c Modified on June 27, 2008 to properly calulate depth for p>2
c modified the ouput 
      one=1.0
      sum=0.
      fac=dble(p)
      coef=dble(n)
      do 10 i=2,p
         fac=fac*dble(p-i+1)
         coef=coef*dble(n-i+1)
10    continue
      coef=coef/fac
      
      if (p.gt.2) then
        call reprow(p,n,p,1,sum,mat,x,u) 
        sum=sum/fac 
      endif
      
      
      if(p.eq.2) then
      do 50 i=1,n
         do 45 j = i,n
           det=(u(1)-x(i,1))*(u(2)-x(j,2))-(u(2)-x(i,2))*(u(1)-x(j,1))
           aire=dabs(det)/fac
           sum=sum+aire
45       continue
50    continue
      endif
      dpth=.5/(1.+sum/coef)
      end

c ------------------------------------------------------------------------

      subroutine determinant(x,p,det)

c        Calcule le determinant de x.
c        p   : entier egal au nombre de dimensions
c        x   : Matrice p*p contenant les n points en p dimensions
c        det : determinant de x
c        Jean-Francois Plante ete 1999


      integer p,ipvt(p),info
      real*8 x(p,p),det,d(2),work(p),s1

      if(p.eq.2) then
         det=x(1,1)*x(2,2)-x(1,2)*x(2,1)
      elseif (p.eq.3) then
      det = x(1,1)*x(2,2)*x(3,3)-x(1,1)*x(2,3)*x(3,2)-x(2,1)*x(1,2)*x(3,
     #3)+x(2,1)*x(1,3)*x(3,2)+x(3,1)*x(1,2)*x(2,3)-x(3,1)*x(1,3)*x(2,2)
      elseif (p.eq.4) then
      s1 = x(1,1)*x(2,2)*x(3,3)*x(4,4)-x(1,1)*x(2,2)*x(3,4)*x(4,3)-x(1,1
     #)*x(3,2)*x(2,3)*x(4,4)+x(1,1)*x(3,2)*x(2,4)*x(4,3)+x(1,1)*x(4,2)*x
     #(2,3)*x(3,4)-x(1,1)*x(4,2)*x(2,4)*x(3,3)-x(2,1)*x(1,2)*x(3,3)*x(4,
     #4)+x(2,1)*x(1,2)*x(3,4)*x(4,3)+x(2,1)*x(3,2)*x(1,3)*x(4,4)-x(2,1)*
     #x(3,2)*x(1,4)*x(4,3)-x(2,1)*x(4,2)*x(1,3)*x(3,4)+x(2,1)*x(4,2)*x(1
     #,4)*x(3,3)
      det= s1+x(3,1)*x(1,2)*x(2,3)*x(4,4)-x(3,1)*x(1,2)*x(2,4)*x(4,3)-x(
     #3,1)*x(2,2)*x(1,3)*x(4,4)+x(3,1)*x(2,2)*x(1,4)*x(4,3)+x(3,1)*x(4,2
     #)*x(1,3)*x(2,4)-x(3,1)*x(4,2)*x(1,4)*x(2,3)-x(4,1)*x(1,2)*x(2,3)*x
     #(3,4)+x(4,1)*x(1,2)*x(2,4)*x(3,3)+x(4,1)*x(2,2)*x(1,3)*x(3,4)-x(4,
     #1)*x(2,2)*x(1,4)*x(3,3)-x(4,1)*x(3,2)*x(1,3)*x(2,4)+x(4,1)*x(3,2)*
     #x(1,4)*x(2,3)
      else
      call dgefa(x,p,p,ipvt,info)
      call dgedi(x,p,p,ipvt,d,work,10)
      det = d(1) * 10.0**d(2)
      endif
      end

c ------------------------------------------------------------------------
        subroutine covmat(p,n,t,s)

c        Calcule un estimateur de la matrice de covariance de t
c        t   : Echantillon de n points en p dimensions (matrice n*p) 
c        p   : entier, nombre de dimension
c        n   : entier, nombre de points
c        s   : matrice p*p = estimateur = sum((Ti-Tbar)(Ti-Tbar)') 
c        Jean-Francois Plante, ete 1999

       integer p,n,i,j,k
       real*8 t(p,n),s(p,p),tbar(p)

c    Calcul de la moyenne

       call moycol(p,n,t,tbar)

c    Calcul Ti - Tbar

       do 20 i=1,p
           do 15 j=1,n
               t(i,j)=t(i,j)-tbar(i)
15         continue
20     continue

c   Calcul de S

      do 30 i=1,p
           do 25 j=i,p
               s(i,j)=0
               s(j,i)=0
               do 22 k=1,n
                   s(i,j)=s(i,j)+t(i,k)*t(j,k)/n
22             continue
               s(j,i)=s(i,j)
25         continue
30     continue

       end

c-------------------------------------------------------------------------------
        subroutine EQM(p,n,t,s)

c        Calcule un estimateur de l'erreur quadratique moyenne de t
c        t   : Echantillon de n points en p dimensions (matrice n*p)
c        p   : entier, nombre de dimension
c        n   : entier, nombre de points
c        s   : matrice p*p = estimateur = sum(Ti*Ti') 
c        Jean-Francois Plante, ete 1999

       integer p,n,i,j,k
       real*8 t(p,n),s(p,p)


c,tbar(p),tt(p,p)

c    Calcul de la moyenne

c       call moycol(p,n,t,tbar)
c       call ttprime(p,tbar,tt)

c    Calcul Ti - Tbar

c       do 20 i=1,p
c           do 15 j=1,n
c               t(i,j)=t(i,j)-tbar(i)
c15         continue
c20     continue

c   Calcul de S

      do 30 i=1,p
           do 25 j=i,p
               s(i,j)=0
               s(j,i)=0
               do 22 k=1,n
                   s(i,j)=s(i,j)+t(i,k)*t(j,k)/n
c+tt(i,j)
22             continue
               s(j,i)=s(i,j)
25         continue
30     continue

       end

c-------------------------------------------------------------------------------
       subroutine moycol(p,n,x,y)

c        Calcule la moyenne des lignes de  t
c        t   : Echantillon de n points en p dimensions (matrice n*p)
c        p   : entier, nombre de dimension
c        n   : entier, nombre de points
c        y   : vecteur de p elements : moyenne des lignes de t
c        Jean-Francois Plante, ete 1999

       integer p,n,i
       real*8 x(p,n),y(p)

       call sumcol(p,n,x,y)
       do 10 i=1,p
          y(i)=y(i)/n
10     continue
       end

c ------------------------------------------------------------------------------ 
       subroutine ttprime(p,t,x)

c        Calcule t*t'
c        t   : Vecteur de p elements
c        p   : entier, longueur de t
c        x   : matrice p*p, =t*t'
c        Jean-Francois Plante, ete 1999

       integer p,i,j
       real*8  t(p),x(p,p)

       do 20 i=1,p
           do 10 j=i,p
               x(i,j)=t(i)*t(j)
               x(j,i)=x(i,j)
10         continue
20     continue
       end
c ------------------------------------------------------------------------------
       subroutine sumcol(p,n,x,y)

c        Calcule la sommes des lignes de  t
c        t   : Echantillon de n points en p dimensions (matrice n*p)
c        p   : entier, nombre de dimension
c        n   : entier, nombre de points
c        y   : vecteur de p elements : somme des lignes de t
c        Jean-Francois Plante, ete 1999

       integer p,n,i,j
       real*8 x(p,n),y(p)

       do 20 i=1,p
           y(i)=0
           do 10 j=1,n
               y(i)=y(i)+x(i,j)
10         continue
20     continue
       end
c ------------------------------------------------------------------------------
       subroutine clv(p,a,u,b,v,x)

c        Calcule la combinaison lineaire x=a*u+b*v)
c        p = longueur des vecteurs 

       integer p,i
       real*8 a,b,u(p),v(p),x(p)

       do 100 i=1,p
              x(i)=a*u(i)+b*v(i)
100    continue

       end

c ------------------------------------------------------------------------------
       subroutine getrow(p,n,x,m,y)

c        Extrait la ligne m d'une matrice n*p
c        p   : entier egal au nombre de dimensions
c        n   : entier egal au nombre de points
c        x   : Matrice n*p
c        m   : entier = indice de la ligne a extraire 
c        y   : vecteur de p elements = colonne desiree
c        Jean-Francois Plante ete 1999

       integer p,n,i,m
       real*8 x(n,p),y(n)

       do 10 i=1,p
           y(i)=x(m,i)
10     continue
       end
c ------------------------------------------------------------------------------
       subroutine putrow(p,n,x,m,y)

c        Remplit la ligne m d'une matrice n*p
c        p   : entier egal au nombre de dimensions
c        n   : entier egal au nombre de points
c        x   : Matrice n*p, la colonne m est modifiee
c        m   : entier = indice de la ligne a remplir 
c        y   : vecteur de p elements a copier sur la ligne m
c        Jean-Francois Plante ete 1999

       integer p,n,i,m
       real*8 x(n,p),y(n)

       do 10 i=1,p
           x(m,i)=y(i)
10     continue
       end
c ------------------------------------------------------------------------------

      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),det(2),work(n)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

c-----------------------------------------------------------------------
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
c -----------------------------------------------------------------------
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(n),dy(n),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
c -----------------------------------------------------------------------
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(n)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
c ------------------------------------------------------------------------
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(n),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
c --------------------------------------------------------------------------
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(n),dy(n),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end

c ====================================================================
c                        ojamed.f
c ====================================================================

        SUBROUTINE OJAMED(X,Y,N,IWS,XMED,YMED,IFAULT,eps)
C
C  On exit (XMED,YMED) will contain the bivariate median of the
C  points (X(i),Y(i)), i=1...N.
C
C  The integer array IWS is used as workspace.
C
c        Pour plus de renseignements, consultez l'article suivant :
c        Niinimaa, A, Oja, H., Nyblom, J. (1992), {\it AS 277 : The 
c        Oja Bivariate Median}, Applied Statistics, vol.41, 611-617.

      REAL*8 X(N), Y(N), XMIN, XMAX, YMIN, YMAX, W
      real*8 YKL, XKL, SMALLD, SMALL, EPS, XMED, YMED, A

C
      INTEGER IWS(N), N,nerr,maxerr
      integer NZ, LL, NZERO, J, KQ, I, NT,L,K, NDIF, IFAULT

      REAL*8 TU, TL, T, WW, DP, D0
        nerr=0
        maxerr=n/200
        if(maxerr.lt.5) maxerr=5
        IFAULT=0
      call indexx(n, x, iws)
        XMIN=X(IWS(1))
        XMAX=X(IWS(N))
      call indexx(n, y, iws)
        YMIN=Y(IWS(1))
        YMAX=Y(IWS(N))
  
C
        IFAULT=0
        SMALL=EPS*(XMAX-XMIN)*(YMAX-YMIN)
        SMALLD=SMALL*N
C
        DO 2 NDIF=1,N
        WW=0.0
        DO 1 K=1,N-1

        DO 1 L=K+1,N
          NT=0
          XKL=X(K)-X(L)
          YKL=Y(K)-Y(L)
          IF(ABS(XKL)+ABS(YKL).LT.SMALL)GO TO 1
            DO 4 I=1,N
            W=(Y(I)-Y(L))*XKL-(X(I)-X(L))*YKL
            WW=WW+ABS(W)
C
C  ABS(W)/2 is the area of the triangle with vertices (X(i),Y(i)),
C  (X(j),Y(j)) and (X(k),Y(k)).
C
            IF(W.GT.0.0)NT=NT+1
4           CONTINUE
          IF((ABS(2*NT-N+2).LE.NDIF).AND.WW.GT.SMALLD) GOTO 3
1       CONTINUE
        IF(WW.LE.SMALLD) THEN
C
C  The data set is completely collinear.
C
            IF (MOD(N,2).EQ.0) THEN
              XMED=(X(IWS(N/2))+X(IWS(N/2+1)))/2
              YMED=(Y(IWS(N/2))+Y(IWS(N/2+1)))/2
              ELSE
              XMED=X(IWS((N+1)/2))
              YMED=Y(IWS((N+1)/2))
              END IF
                ifault=-1
                RETURN
            END IF

2         CONTINUE

3       CALL TUTL(X,Y,N,XMIN,YMIN,XMAX,YMAX,K,L,SMALL,TU,TL)

C Start the search along the line through the points  X(k),Y(k))
C and (X(k),Y9k)), dividing the N-2 other pints evenly.

100    KQ=1

C  Search for the minimum of the Oja objective function on the line
C  (X(k),Y(k)) to (X(l),Y(l)).

        DO 11 I=1,N-1
        DO 11 J=I+1,N
        WW=(X(K)-X(L))*(Y(J)-Y(I))-(Y(K)-Y(L))*(X(J)-X(I))
        IF(ABS(WW).LT.SMALL)GO TO 11
          T=((X(L)-X(J))*(Y(I)-Y(J))-(Y(L)-Y(J))*(X(I)-X(J)))/WW
        IF(T.LT.TU.AND.T.GT.TL)THEN
          CALL DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
            A=DP+D0
          IF(D0.GE.ABS(DP)-SMALLD) GO TO 400
          IF(DP+D0.LT.0.0)TL=T
          IF(DP+D0.GT.0.0)TU=T
        END IF
11      CONTINUE

C  The local minimum on the line (X(k),Y(k))-(X(l),Y(l)) was not
C  found.

        IFAULT=IFAULT+2

        IF(IFAULT.GT.10) then
           xmed=y(1)
           ymed=smalld
           RETURN
        endif

        L=N-1
        K=N

C  The next line has been found

400     KQ=0
        I=K
        J=L
        nerr=nerr+1
        if(nerr.gt.maxerr) goto 87
C
            LL=1
            NZ=NZERO
200         K=INT(IWS(LL)/REAL(N))
            L=MOD(IWS(LL),N)+1
            WW=(X(K)-X(L))*(Y(J)-Y(I))-(Y(K)-Y(L))*(X(J)-X(I))
            IF(ABS(WW).GE.SMALL) THEN
            T=((X(L)-X(J))*(Y(I)-Y(J))-(Y(L)-Y(J))*(X(I)-X(J)))/WW
            CALL DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
            IF (D0.LT.ABS(DP)-SMALL) THEN
            CALL TUTL(X,Y,N,XMIN,YMIN,XMAX,YMAX,K,L,SMALL,TU,TL)
              IF(DP+D0.LT.0.0)THEN
                TL=T
                ELSE
                TU=T
              END IF
                GO TO 100
            END IF
            LL=LL+1
            END IF
C
        IF(LL.EQ.NZ+1)THEN
            XMED=X(L)+T*(X(K)-X(L))
            YMED=Y(L)+T*(Y(K)-Y(L))
            RETURN
C
        END IF
C
        GO TO 200

87      call ojamed2(x,y,n,dble(1),dble(0),dble(0),xmed,ymed,eps)

        END
c----------------------------------------------------------------------
      subroutine ojadepth2(z,u,v,n,rep)

      real*8 z(n,2),u,v,w(2),rep
      integer n

      w(1)=u
      w(2)=v
      call ojadepth(z,w,2,n,rep)
      end

c----------------------------------------------------------------------
      subroutine ojamed2(wx,wy,n,dx,x0,y0,xx,yy,eps)

c        Calcule la mediane de Oja a eps pres
c        x  : composante en x des donnees
c        y  : composante en y des points
c        n  : nombre de points
c        dx : pas initial
c        x0 : Point de depart algo.
c        y0 : point de depart algo.
c        xx : abcisse mediane
c        yy : ordonee mediane
c        eps: precision

      integer n,i,gros
      real*8 wx(n),wy(n),z(n,2),xx,yy,eps,x0,y0,dx,x,y,stp,od,oldod,k

      do 10 i=1,n
            z(i,1)=wx(i)
            z(i,2)=wy(i)
10    continue

      k=dble(1)
      gros=2
      x=x0 
      y=y0
      stp=dble(1)

      do 20 
            
            call ojadepth2(z,x,y,n,od)
            call ojadepth2(z,x+k*dx,y,n,oldod)
            if(od.lt.oldod) stp=1
            call ojadepth2(z,x-k*dx,y,n,oldod)
            if(od.lt.oldod) stp=-1
            oldod=od
            do 15 

                  x=x+stp*k*dx
                  oldod=od
                  call ojadepth2(z,x+k*stp*dx,y,n,od)


                  if(od.le.oldod) goto 16
15          continue
16          continue

            call ojadepth2(z,x,y,n,od)
            call ojadepth2(z,x,y+k*dx,n,oldod)
            if(od.lt.oldod) stp=1
            call ojadepth2(z,x,y-k*dx,n,oldod)
            if(od.lt.oldod) stp=-1
            oldod=od
            do 25 
            

                  y=y+stp*k*dx
                  oldod=od
                  call ojadepth2(z,x,y+k*stp*dx,n,od)
                  if(od.le.oldod) goto 26
25          continue
26          continue

            if (2*dx*k.le.eps) goto 21
            k=k/dble(2)


20    continue
21    continue
      xx=x
      yy=y
      end

c -------------------------------------------------------------------
        SUBROUTINE DER(X,Y,N,K,L,T,IWS,SMALL,NZERO,KQ,DP,D0)
C
C  DP-D0 and DP+D0 are the directional derivatives of the Oja
C  objective function just before and after the point
C  (X(l)+t*(X(k)-X(l),Y(l)+t*(Y(k)-Y(l)).
C
        REAL*8 X(N), Y(N), DIF, YKL, XKL, SMALL, SMALLD
        REAL*8 TT, T, WW, DP, D0, SGN
        INTEGER IWS(N), NZERO, II, JJ, N, K, L, KQ
C
        smalld=1e-6

C
        DP=0.0
        D0=0.0
        NZERO=0
        XKL=X(K)-X(L)
        YKL=Y(K)-Y(L)
C
        DO 1 II=1,N-1
        DO 1 JJ=II+1,N
        WW=(Y(JJ)-Y(II))*XKL-(X(JJ)-X(II))*YKL
        IF(ABS(WW).LE.SMALL)GO TO 1
        TT=((X(L)-X(JJ))*(Y(II)-Y(JJ))-(Y(L)-Y(JJ))*(X(II)-X(JJ)))/WW
        DIF=T-TT
C
        IF(ABS(DIF).LE.SMALLD)THEN
            NZERO=NZERO+1
            NZERO=MIN(NZERO,N)
            IF(KQ.NE.0)IWS(NZERO)=N*II+JJ-1
            D0=D0+ABS(WW)
            ELSE
            DP=DP+ABS(WW)*SGN(DIF)
          END IF
  1       CONTINUE
          RETURN
          END
c -------------------------------------------------------------------
          FUNCTION SGN(X)
          REAL*8 SGN, X
          SGN=1.0
          IF (X.LT.0.0)SGN=-1.0
          RETURN
          END
c -------------------------------------------------------------------
        SUBROUTINE TUTL(X,Y,N,XMIN,YMIN,XMAX,YMAX,K,L,SMALL,TU,TL)
C
C  This subroutine calculates the upper (TU) and lower limit (TL) for
C  parameter T on the line (X(l)+t*(X(k)-X(l)), Y(l)+t*(Y(k)-Y(l))
C  inside the rectangle with the vectices (XMIN, YMIN), (XMIN, YMAX),
C  (XMAX, YMIN) and (XMAX, YMAX).
C
      REAL*8 X(N), Y(N), T1, T2, T3, T4, VBIG, SMALL
      real*8 XMIN, YMIN, XMAX, YMAX, TU, TL
      INTEGER K, L, N
C
        vbig=1e38
C
        T1=-VBIG
        T2=-T1
        IF(ABS(X(K)-X(L)).GT.SMALL)THEN
          T1=-(X(L)-XMIN)/(X(K)-X(L))
          T2=(XMAX-X(L))/(X(K)-X(L))
        END IF
        T3=-VBIG
        T4=-T3
        IF(ABS(Y(K)-Y(L)).GT.SMALL)THEN
          T3=-(Y(L)-YMIN)/(Y(K)-Y(L))
          T4=(YMAX-Y(L))/(Y(K)-Y(L))
        END IF
        TU=MIN(MAX(T1,T2),MAX(T3,T4))+SMALL
        TL=MAX(MIN(T1,T2),MIN(T3,T4))-SMALL
        RETURN
        END

c -------------------------------------------------------------------
      SUBROUTINE INDEXX(N,ARRIN,INDX)

      integer INDX(N)
      integer n,j,l,ir,indxt,i
      real*8 arrin(n),q

      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

c ====================================================================
c                        medctr.f
c ====================================================================

       subroutine medctr (m,n,x,norm,c,med,ifault,eps,maxit)
c
c        Algorithm AS 143  Appl. Statist. (1979) Vol. 28, No. 3
c
c        The Mediancentre
c        At least 10 times faster than AS 78.
c
c        Parameters:
c   
c        m       integer     input: number of dimensions
c        n       integer     input: number of sample points
c        x       real array  input: data, n points, m dimensions
c        norm    real        output: norm of gradient
c        c       real        output: multiplicity of m in S ; 0 if m not in S
c        med     real array  output: coordinates of the mediancentre
c        ifault  integer     output: -1 if m or n out of range
c                                     0 if matrix of second derivatives
c                                       is positive definite
c                                     1 if matrix of second derivatives
c                                       is not positive definite
c
      real*8 s,t,c,eps,zero,one,f,r,ang,rr,tmp
      real*8 x(n,m),med(m),md(m),g(m),q(m,m)
      integer l(n),m,n,ifault,i,nn,j,k,nt,j1,j2,itdone,maxit

c      dimension x(maxn,maxm), med(maxm), md(maxm), g(maxm),
c     *          q(maxm,maxm), l(maxn)
      real*8 norm
c
      data zero, one /0.0d0, 1.0d0/
c
      itdone=0
      ifault = -1
      if (m .lt. 1.or. n .lt. 1)
     *    return
      ifault = 0
      do 1 i = 1, m
    1 med(i) = zero
      nn = n
      do 2 i = 1, n
    2 l(i) = i
c
c        computation of (residual) mean
c
    3 f = one / dble(nn)
      do 5 i = 1, m
        s = zero
        g(i) = zero
        do 4 j = 1, nn
          k = l(j)
          s = s + x(k, i)
    4   continue
        md(i) = s * f
        med(i) = med(i) + md(i)
    5 continue
c
c        Computation of function, gradient, and norm of gradient
c
      c = zero
      f = zero
      s = zero
      do 9 i = 1, n
        t = zero
        do 6 j = 1, m
          x(i, j) = x(i, j) - md(j)
          t = t + x(i, j) * x(i, j)
    6   continue
        if (t .gt. zero) goto 7
        c = c + one
        goto 9
    7   t = sqrt(t)
        r = one / t
        s = s + r
        f = f + t
        do 8 j = 1, m
          g(j) = g(j) - r * x(i, j)
    8   continue
    9 continue
      norm = zero
      do 10 i = 1, m
        norm = norm + g(i) * g(i)
   10 continue
      norm = sqrt(norm)
c
c        Check for extremal points
c
      if (norm .le. c + eps) then

         return
      endif
      if (nn .eq. 1) goto 13
c
c        Simplex
c
      nt = nn
      nn = 0
      do 12 i = 1, nt
        k = l(i)
        ang = zero
        do 11 j = 1, m
          ang = ang + g(j) * x(k, j)
   11   continue
        if (ang .gt. zero) goto 12
        nn = nn + 1
        l(nn) = k
   12 continue
      if (nn .gt. 0) goto 3
c
c        Starting value
c
   13 r = (c / norm - one) / s
      do 14 i = 1, m
        md(i) = r * g(i)
   14 continue
c
c        Newton - Raphson procedure
c
   15 c = zero
      f = zero
      s = zero
      do 16 i = 1, m
        g(i) = zero
        med(i) = med(i) + md(i)
        do 16 j = i, m
          q(i, j) = zero
   16 continue
      do 19 i = 1, n
        t = zero
        do 17 j = 1, m
          x(i, j) = x(i, j) - md(j)
          t = t + x(i, j) * x(i, j)
   17   continue
        if (t .eq. zero) goto 19
        t = sqrt(t)
        r = one / t
        s = s + r
        rr = r * r * r
        f = f + t
        do 18 j = 1, m
          g(j) = g(j) - r * x(i, j)
          do 18 k = j, m
            q(j, k) = q(k, j) - rr * x(i, j) * x(i, k)
   18   continue
   19 continue
      norm = zero
      do 20 j = 1, m
        norm = norm + g(j) * g(j)
   20 continue
      norm = dsqrt(norm)
      if (norm .le. eps) then

         return
      endif
c
c        Cholesky and solution of equation
c
      ifault = 1
      if (q(1,1) + s .le. zero) then

         return
       endif
      q(1,1) = sqrt(q(1,1) + s) 
      do 21 i = 2, m
        q(i, 1) = q(1, i) / q(1, 1)
   21 continue
      do 25 j = 2, m
        tmp = s + q(j, j)
        j1 = j - 1
        do 22 i = 1, j1
          tmp = tmp - q(j, i) * q(j, i)
   22   continue
        if (tmp .le. zero) then

           return
        endif
        q(j, j) = sqrt(tmp)
        j2 = j + 1
        if (j2 .gt. m) goto 25
        do 24 k = j2, m
          tmp = q(j, k)
          do 23 i = 1, j1
            tmp = tmp - q(j, i) * q(k, i)
   23     continue
          q(k, j) = tmp / q(j, j)
   24   continue
   25 continue
      md(1) = -g(1) / q(1, 1)
      do 27 j = 2, m
        j1 = j - 1
        tmp = -g(j)
        do 26 i = 1, j1
          tmp = tmp - q(j, i) * md(i)
   26   continue
        md(j) = tmp / q(j, j)
   27 continue
      md(m) = md(m) / q(m, m)
      j1 = m - 1
      do 29 j2 = 1, j1
        i = m - j2
        tmp = md(i)
        k = i + 1
        do 28 j = k, m
          tmp = tmp - q(j, i) * md(j)
   28   continue
        md(i) = tmp / q(i, i)
   29 continue
      ifault = 0
      itdone=itdone+1
      if (itdone.gt.maxit) then 
         ifault=-2
         return
      endif
      goto 15
      end

c ------------------------------------------------------------------------

      SUBROUTINE MEDctr78(X, Y, N, IP, IT, IFAULT)
C
C     ALGORITHM AS 78  APPL. STATIST. (1974) VOL.23, NO.3
C
C     The mediancentre, generalising the median, is the point with
C     minimum total distance from a set of multivariate samples
C
      REAL*8  X(N, IP), Y(IP), corner,xx
      REAL*8  Z(ip), LAMBDA, LEPSI, LEPSR, LEPSD
      integer n,ip,it,ifault,icount,lcount,ll,ii,i,j,k,l,lc
      real*8 diam,s,epsi,u1,u2
C
      real*8  C(ip), COMP, D, DD, DELTA, EPSR, EPSD, SLAM,
     +  ZERO, ONE
C
      DATA LEPSD /0.0001/, LEPSR / 0.00001/, LEPSI /0.000001/
      DATA ICOUNT /100/, LCOUNT /50/
      DATA ZERO /0.D0/, ONE /1.D0/
C
C     Initial settings
C
      
      IFAULT = 0
      LL = 0
      II = 1

      IF (N .EQ. 1) GO TO 25
      IF (N .LE. 0 .OR. IP .LE. 0) GO TO 5
  

C
C     Calculate the diameter
C
      DIAM = 0.0
      DO 2 I = 2, N
        DO 2 J = 1, I-1
          S = 0.0
          DO 1 K = 1, IP
    1     S = S + (X(I,K) - X(J,K))**2
          DIAM = MAX(S, DIAM)
    2 CONTINUE
      DIAM = dSQRT(DIAM)
      EPSR = LEPSR * DIAM
      EPSI = LEPSI * DIAM
      EPSD = LEPSD * DIAM
C
C     Initial median centre = the centroid
C
      U1 = 1.0 / dble(N)
      DO 4 J = 1, IP
        S = 0.0
        DO 3 I = 1, N
    3   S = S + X(I,J)
        Y(J) = S * U1
    4 CONTINUE
      IT = ICOUNT
      IF (IP .LE. 50) GO TO 6
    5 IFAULT = 1
      IT = 0
      RETURN
C
C     Main iterative loop
C
    6 DO 23 L = 1, ICOUNT
C
C     Direction cosines and resultant
C
        CORNER = 0.0
        DO 7 J = 1, IP
    7   C(J) = ZERO
        LL = L
        DO 11 I = 1, N
          D = ZERO
          DO 8 J = 1, IP
    8     D = D + (X(I,J) - Y(J))**2
          DD = SQRT(D)
          IF (DD .GT. EPSD) GO TO 9
          CORNER = CORNER + 1.0
          II = I
          GO TO 11
    9     D = ONE / DD
          DO 10 J = 1, IP
   10     C(J) = C(J) + (X(I,J) - Y(J)) * D
   11   CONTINUE
        D = ZERO
        DO 12 J = 1, IP
   12   D = D + C(J)**2
        D = SQRT(D)
        DD = D
C
C     Tests for zero resultant or degenerate solution
C
        IF (CORNER .EQ. 0.0) GO TO 13
        IF (D .LE. CORNER) GO TO 25
        D = D - CORNER
   13   IF (D .LE. EPSR) GO TO 24
        DD = ONE / DD
        DO 14 J = 1, IP
   14   C(J) = C(J) * DD
C
C     Step by bisection to give zero component at lambda
C
        U1 = 0.0
        U2 = DIAM
        DO 20 LC = 1, LCOUNT
          COMP = ZERO
          LAMBDA = 0.5 * (U1 + U2)
          SLAM = LAMBDA * LAMBDA
          DO 15 J = 1, IP
   15     Z(J) = Y(J) + LAMBDA * C(J)
          DO 17 I = 1, N
            DELTA = ZERO
            D = SLAM
            DO 16 J = 1, IP
              XX = X(I,J)
              D = D - (XX - Y(J))**2
              DELTA = DELTA + (XX - Z(J))**2
   16       CONTINUE
            DD = SQRT(DELTA)
            IF (DD .LT. EPSD) GO TO 21
            COMP = COMP - (D + DELTA) / DD
   17     CONTINUE
          IF (COMP .GT. ZERO) GO TO 18
          U2 = LAMBDA
          GO TO 19
   18     U1 = LAMBDA
   19     IF ((U2 - U1) .LE. EPSI) GO TO 21
   20   CONTINUE
C
   21   DO 22 J = 1, IP
   22   Y(J) = Y(J) + C(J) * LAMBDA
   23 CONTINUE
C
   24 IT = LL
      RETURN
   25 IT = -LL
      DO 26 J = 1, IP
   26 Y(J) = X(II,J)
      RETURN
      END

c ====================================================================
c                        Tukmed.f
c ====================================================================

      subroutine HALFMED(x,y,n,dpth,tukmed,xcont,ycont,ncont,
     + dpths,ndpth,tm,maxnum,err,eps,dithfactor,maxdith,mustdith,
     + missing,factor)
C
C  Version: 18 December 1998
C
c        modifie ete 1999, Jean-Francois Plante
c
c        Calcule la mediane de Tukey ou les contours de profondeurs de 
c        l'echantillon (x,y).
c        x,y : Coordonnees de l'echantillon
c        n   : Nbre de points
c        dpth: Profondeur de la mediane de Tukey
c        tukmed: Coordonnees de la mediane de Tukey
c        xcont,ycont: Coordonnees des contours de profondeur
c        ncont: vecteur d'entiers, nbre de points pour chacun des contour
c        dpths: vecteur d'entiers = profondeur des contours a calculer
c        ndpth: longueur de dpths
c        tm  : 0 => Calcule mediane, 1 => Calcule contours
c        maxnum: floor(4 * n * sqrt(n) + 1)
c        err,missing: voir Tukmed et Isod de depthtools.S
c        eps : tolerance numerique
c        dithfactor,maxdith, mustdith : voir depthtools.ps
c        factor : voir depthtools.ps
c
c        Pour plus de renseignements, consultez les articles suivants,
c        disponible au http://win-www.uia.ac.be/u/statis/ : 
c        Rousseeuw, P.J. and Ruts, I. (1998), {\it Constructing the bivariate
c        Tukey median}, Statistica Sinica, vol.8, 828-839.
c        Ruts, I. and Rousseeuw, P.J. (1996), {\it Computing depth contours
c        of bivariate point clouds}, Computational Statistics and Data
c        Analysis, vol.23, 153-168. 

      integer n,dpth,maxnum,iv,nceil,dith,err,kk,ndpth
      INTEGER NCIRQ(N),MCIRQ(N),NRANK(N),F(N),count,sdep
      integer JLV(N),JRV(N),maxdith,xind,mustdith,maxd
      integer dpths(ndpth),ncont(ndpth),t,missing(ndpth)
      INTEGER IND1(n*(n-1)/2),IND2(n*(n-1)/2),nrun
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM),KORNR(MAXNUM,4)
      INTEGER LEFT,KOUNT,NUM,hdep1,khulp,moredith
      INTEGER I,J,L,M,empty,ib,ie,le,tm,notingp
      REAL*8 X(N),Y(N),WX(N),WY(N),dpf(n),dithfactor
      REAL*8 ANGLE(n*(n-1)/2),D(n*(n-1)/2),alpha(maxnum),beta(n)
      REAL*8 PI,PI2,EPS,xcord1,ycord1,xsum,ysum,xcordp,ycordp
      REAL*8 XCORD,YCORD,E1,E2,F1,F2,G1,G2,ANG1
      real*8 sum,tukmed(2),wx1(n),wy1(n),rand(2)
      real*8 xcont(n*(n-1)/2),ycont(n*(n-1)/2),factor
 

      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      nrun=maxnum

      dith=0
      xind=0
      t=1
      maxd=int(real(n)/2)  

      DO 5 I=1,N
         WX(I)=X(I)
         WY(I)=Y(I)
         wx1(i)=x(i)
         wy1(i)=y(i)
         NCIRQ(I)=I
         MCIRQ(I)=I
 5    CONTINUE

C
C  To test whether the data points are in general position,
C  we first check whether any two points coincide.
C  When we sort the points with respect to the x-coordinate, we
C  permute NCIRQ and MCIRQ in the same way, to initialize them.
C
      notingp=0
1      CALL SORT3(WX,NCIRQ,MCIRQ,WY,N)
      I=0
10    I=I+1
      IF (I+1.GT.N) GOTO 15
      J=I+1
      IF (WX(I).NE.WX(J)) GOTO 10
      IF (WY(I).EQ.WY(J)) THEN
         GOTO 209
      ELSE
         if (wy(i).lt.wy(j)) then
            iv=mcirq(j)
            mcirq(j)=mcirq(i)
            mcirq(i)=iv
         endif
         IF (J+1.LE.N) THEN
            L=J+1
            IF (WX(I).EQ.WX(L)) THEN
               GOTO 209
            ELSE
               GOTO 10
            ENDIF
         ENDIF
      ENDIF
C    
C  Compute all the angles formed by pairs of data points.
C
15    M=((N*(N-1))/2)
      L=1
      DO 20 I=1,N
         DO 25 J=I+1,N
            IF (X(I).EQ.X(J)) THEN
               ANGLE(L)=PI2
            ELSE
               ANGLE(L)=DATAN((Y(I)-Y(J))/(X(I)-X(J)))
               IF (ANGLE(L).LE.0.0) ANGLE(L)=ANGLE(L)+PI
            ENDIF
            IND1(L)=I
            IND2(L)=J
            L=L+1
25    CONTINUE
20    CONTINUE
C 
C  Sort all the angles and permute IND1 and IND2 in the same way.
C  To avoid using several SORT-routines, we will always permute
C  two integer arrays and one real array.
C
      CALL SORT3(ANGLE,IND1,IND2,D,M)
C 
C  Test whether any three points are collinear
C
      LEFT=1
30    ANG1=ANGLE(LEFT)
      DO 35 J=LEFT+1,M
         IF (ANGLE(J).GT.ANG1) THEN
            LEFT=J
            GOTO 30
         ELSE
            DO 36 I=LEFT,J-1
               IF ((IND1(I).EQ.IND1(J)).or. 
     +                  (IND1(I).EQ.IND2(J))) THEN 
               GOTO 209
               ENDIF
               IF ((IND2(I).EQ.IND1(J)).or.
     +            (IND2(I).EQ.IND2(J))) THEN
               GOTO 209
               ENDIF
36          CONTINUE
         ENDIF
35    CONTINUE
      goto 37
C
C     If the data are not in general position, use dithering
C
209   if (mustdith.eq.0) then
         err = -1
         goto 610
      endif
210   notingp=1

      if (dith.ge.maxdith.and.tm.eq.0) then
         err=1
         goto 610
      endif
     
      do 211 i=1,n
         ncirq(i)=i
         mcirq(i)=i
         call randm(nrun,rand(1))
         call randm(nrun,rand(2))
         wx(i)=wx(i)+(rand(1)-.5)*eps*dithfactor
         x(i)=wx(i)
         wy(i)=wy(i)+(rand(2)-.5)*eps*dithfactor
         y(i)=wy(i)
      
211   continue
      err=-1
      dith=dith+1
      moredith=0
      goto 1
37    continue

      if (tm.eq.1) then
         goto 200
      endif

C
C     Calculation of the Tukey median
C     
 
      xsum=0
      ysum=0
      tukmed(1)=0
      tukmed(2)=0
      if (n.le.3) then
         do 41 i=1,n
         xsum=xsum+x(i)
         ysum=ysum+y(i)
41       continue
         tukmed(1)=xsum/n
         tukmed(2)=ysum/n
         dpth=1
         goto 610   
      endif

      ib=nceil(n,3)
      ie=int(real(n)/2)
171     le=ie-ib
        if (le.eq.0) goto 185
        if (n.gt.150) then
           khulp=nceil(n,3)-20
        else
           khulp=0
        endif
      CALL ISODEPTH(N,M,X,Y,N,n*(n-1)/2,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,ib+nceil(le,2),khulp,EMPTY,moredith,eps
     +  ,tm,factor)
        if (moredith.eq.1) then
            if(mustdith.ne.0) then
                goto 210
            else
                err=ib+nceil(le,2)
            endif
        endif
        if (empty.eq.1) ie=ib+nceil(le,2)
        if (empty.eq.0) ib=ib+nceil(le,2)
        if (le.eq.1) goto 185
        goto 171

179     empty=0
        ie=int(real(n)/2)
180     le=ie-ib
        if (le.eq.0) goto 185
      CALL ISODEPTH(N,M,X,Y,n,n*(n+1)/2,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,ib+nceil(le,2),khulp,EMPTY,moredith,eps
     +  ,tm,factor)
        if (moredith.eq.1)  then
            if(mustdith.ne.0) then
                goto 210
            else
                err=ib+nceil(le,2)
            endif
        endif

        if (empty.eq.1) ie=ib+nceil(le,2)
        if (empty.eq.0) ib=ib+nceil(le,2)
        if (le.eq.1) goto 185
        goto 180
185   CALL ISODEPTH(N,M,X,Y,n,n*(n+1)/2,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,ib,khulp,EMPTY,moredith,eps
     +  ,tm,factor)
        if (moredith.eq.1)  then
            if(mustdith.ne.0) then
                goto 210
            else
                err=ib
            endif
        endif
      if(moredith.eq.1) goto 610
      dpth=ib
C
C  Scan KORNR and write the coordinates of the vertices to the file DK.DAT.
C
186   KOUNT=0
      I=1
      E1=Y(KORNR(I,2))-Y(KORNR(I,1))
      F1=X(KORNR(I,1))-X(KORNR(I,2))
      G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +   -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
      E2=Y(KORNR(I,4))-Y(KORNR(I,3))
      F2=X(KORNR(I,3))-X(KORNR(I,4))
      G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +   -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
      XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
      YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)

      if (tm.ne.0) then 
         xcont(xind+kount+1)=xcord
         ycont(xind+kount+1)=ycord
      endif

      if (tm.eq.0) then
         wx(kount+1)=xcord
         wy(kount+1)=ycord
         xsum=xcord
         ysum=ycord
      endif
      xcord1=xcord
      ycord1=ycord
      xcordp=xcord
      ycordp=ycord
      KOUNT=KOUNT+1
      I=I+1
      if (num.eq.1) goto 195 
190   IF ((KORNR(I,1).EQ.KORNR(I-1,1)).AND.(KORNR(I,2).EQ.KORNR(I-1,2))
     +.AND.(KORNR(I,3).EQ.KORNR(I-1,3).AND.KORNR(I,4).EQ.KORNR(I-1,4)))
     +THEN
        I=I+1
      ELSE
        IF ((KORNR(I,1).EQ.KORNR(1,1)).AND.(KORNR(I,2).EQ.KORNR(1,2))
     +    .AND.(KORNR(I,3).EQ.KORNR(1,3).AND.KORNR(I,4).EQ.KORNR(1,4))) 
     +        THEN
          GOTO 195
        ELSE
          E1=Y(KORNR(I,2))-Y(KORNR(I,1))
          F1=X(KORNR(I,1))-X(KORNR(I,2))
          G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +       -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
          E2=Y(KORNR(I,4))-Y(KORNR(I,3))
          F2=X(KORNR(I,3))-X(KORNR(I,4))
          G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +       -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
          XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
          YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)
          if (((dabs(xcord-xcordp).lt.eps).and.
     +        (dabs(ycord-ycordp).lt.eps)).or. 
     +       ((dabs(xcord-xcord1).lt.eps).and.
     +        (dabs(ycord-ycord1).lt.eps))) then
             i=i+1
          else
             xcordp=xcord
             ycordp=ycord

      if (tm.ne.0) then 
         xcont(xind+kount+1)=xcord
         ycont(xind+kount+1)=ycord
      endif


      if (tm.eq.0) then
             wx(kount+1)=xcord
             wy(kount+1)=ycord
             xsum=xsum+xcord
             ysum=ysum+ycord
      endif 
           KOUNT=KOUNT+1
          I=I+1
          endif
        ENDIF
      ENDIF
      IF (I.NE.(NUM+1)) GOTO 190


      count=0
      if (moredith.eq.1) then
         j=kount-1

         do 194 i=0,j
            if (tm.ne.0) then 
               xcord=xcont(xind+i+1)
               ycord=ycont(xind+i+1)
            endif

           if (tm.eq.0) then
              xcord=wx(i+1)
              ycord=wy(i+1)
           endif

           call tukdepth(xcord,ycord,n,x,y,beta,f,DPF,
     +  JLV,JRV,HDEP1,sdep,eps)
c           hdep1=int(hd*dble(n)+.5)

           if (kk.le.hdep1) then


              if (tm.ne.0) then 
                xcont(xind+count+1)=xcont(xind+i+1)
                ycont(xind+count+1)=ycont(xind+i+1)
              endif

              if (tm.eq.0) then
                 wx(count+1)=wx(i+1)
                 wy(count+1)=wy(i+1)
              endif 

              count=count+1

           else
              kount=kount-1
           endif

194      continue



c            IF (X(I).EQ.X(J)) THEN
c               ANGLE(L)=PI2
c            ELSE
c               ANGLE(L)=DATAN((Y(I)-Y(J))/(X(I)-X(J)))
c               IF (ANGLE(L).LE.0.0) ANGLE(L)=ANGLE(L)+PI
c            ENDIF



      endif


c    
c     Calculation of the center of gravity
c     
195   if(tm.ne.0) then
         ncont(t)=kount
         xind=xind+kount
         goto 300
      endif

      if (tm.eq.0) then
         if (kount.gt.1) then
         do 205 i=1,kount
         wx(i)=wx(i)-(xsum/kount)
         wy(i)=wy(i)-(ysum/kount)
205      continue
         sum=0
         tukmed(1)=0
         tukmed(2)=0
         do 206 i=1,kount-1
         sum=sum+dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i))
         tukmed(1)=tukmed(1)+
     +  ((wx(i)+wx(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(i)+wy(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
206      continue
         sum=sum+dabs(wx(kount)*wy(1)-wx(1)*wy(kount))
         tukmed(1)=tukmed(1)+
     +  ((wx(kount)+wx(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(kount)+wy(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(1)=(tukmed(1)/(3*sum))+(xsum/kount)
         tukmed(2)=(tukmed(2)/(3*sum))+(ysum/kount)
         else
         tukmed(1)=xsum
         tukmed(2)=ysum
         endif
         call tukdepth(tukmed(1),tukmed(2),
     +              n,x,y,beta,f,DPF,JLV,JRV,HDEP1,sdep,eps)
c           hdep1=int(hd*dble(n)+.5)
          
          if ((hdep1.lt.ib).and.(n.gt.150)) then
             khulp=khulp-50
             goto 179
          endif
          if ((hdep1.lt.ib).and.(n.le.150)) then
            khulp=khulp-10
            goto 179
         endif
      endif
      goto 610

200   kk=dpths(t)
      if((n.eq.1).and.(kk.eq.1)) goto 201
      if(kk.gt.maxd.and.dble(kk).gt.factor*dble(n/2)) then
         kount=0
         missing(t)=-1
         goto 195
      endif
201   if (n.le.3) then
         do 202 i=1,n
            xcont(xind+i)=x(i)
            ycont(xind+i)=y(i)
            
202      continue

         ncont(t)=n
         xind=xind+n
         goto 300
      endif

      if(n.gt.150) khulp=nceil(n,3)-20

      CALL ISODEPTH(N,M,X,Y,N,n*(n-1)/2,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,kk,khulp,EMPTY,moredith,eps
     +  ,tm,factor)

      if (empty.eq.1) then
          missing(t)=-1
          if(kk.lt.maxd.and.dble(kk).gt.factor*dble(n/2)) maxd=kk
          goto 300
      endif

      if (moredith.eq.1)  then
          if(mustdith.ne.0) then
             if (dith.ge.maxdith) then
                 missing(t)=-3
                 goto 300 
             else
                 goto 210
             endif
          else
             missing(t)=-2
          endif
      else
          missing(t)=0
      endif

      if (empty.eq.1) then
          ncont(t)=0

      endif
      if (empty.eq.0) goto 186
300   t=t+1
      dith=0

      if(t.le.ndpth) goto 200

610   END
     
      INTEGER FUNCTION NCEIL(M,J)
      integer m,j
      IF (MOD(M,J).EQ.0) THEN
         NCEIL=INT(REAL(M)/J)
      ELSE
         NCEIL=NINT(REAL(M)/J+0.5)
      ENDIF
      RETURN
      END

c ---------------------------------------------------------------------

      SUBROUTINE ISODEPTH(N,M,X,Y,maxn,maxm,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,K,khulp,EMPTY,moredith,eps,tm,factor)
C
C  Computes the depth contour of depth k. This subroutine was described
C  in: Ruts, I., and Rousseeuw, P.J. (1996). Computing depth contours
C  of bivariate point clouds. Computational Statistics and Data Analysis
C  23, 153-168.
C
      integer maxn,maxm,maxnum,jlv,jrv,firstiw2,tm
      INTEGER NCIRQ(N),MCIRQ(N),NRANK(N),F(N)
      DIMENSION JLV(N),JRV(N)
      INTEGER IND1(M),IND2(M),hdep2,hdep3,hdep4,hdep5
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM),KORNR(MAXNUM,4)
      INTEGER KON,KONTROL,NDATA,NDK,HALT,halt2,jj,JFULL,EMPTY
      INTEGER IV,IW1,IW2,NEXT,JFLAG,KOUNT,NUM
      INTEGER HDEP1,I,J,K,L,M,N,moredith,sdep
      INTEGER zoek,iw2nu,zold,khulp,NLOOP,iw1old,iw2old
      integer overschr,firstiw1,iw1oldloop
      REAL*8 X(N),Y(N),BETA(N)
      REAL*8 ANGLE(M),D(M),ALPHA(MAXNUM),DPF(N)
      REAL*8 PI,PI2,EPS,factor
      REAL*8 XCORD,YCORD,ANG1,m1,m2,xcord1,ycord1
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0

      empty=0
C
C   (Re)initialize NCIRQ and NRANK
C
      DO 45 I=1,N
         NCIRQ(I)=MCIRQ(I)
45    CONTINUE
      DO 50 I=1,N
         IV=NCIRQ(I)
         NRANK(IV)=I
 50   CONTINUE
C
C  Let the line rotate from zero to ANGLE(1)
C
      KOUNT=1
      HALT=0
      if (angle(1).gt.pi2) then
         l=1
         CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                          K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
         halt=1
      endif
      L=2
 60   KONTROL=0
      IF ((PI.LE.(ANGLE(L)+PI2)).AND.((ANGLE(L)-PI2).LT.ANGLE(1))) THEN
         CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                          K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
         KONTROL=1
      ENDIF
      L=L+1
      IF (KONTROL.EQ.1) HALT=1
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.1)) THEN
         JFLAG=1
         GOTO 79
      ENDIF
      IF (((HALT.EQ.1).AND.(KONTROL.EQ.0)).OR.(L.EQ.M+1)) THEN
         GOTO 70
      ELSE
         GOTO 60
      ENDIF
 70   if (l.gt.1) then
         JFLAG=L-1
      else
         jflag=m
      endif
      J=0
C
C  In case the first switch didn't occur between zero and ANGLE(1),
C  look for it between the following angles.
C
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.0)) THEN
         HALT=0
         halt2=0
 73      J=J+1
         if (j.eq.m+1) j=1
         L=J+1
         if (l.eq.m+1) l=1
 75      KONTROL=0
         IF ((ANGLE(L)+PI2).LT.PI) THEN
            ANG1=ANGLE(L)+PI2
         ELSE
            ANG1=ANGLE(L)-PI2
         ENDIF
         if (j.eq.m) then
            jj=1
            if (halt2.eq.0) angle(1)=angle(1)+pi
         else
            jj=j+1
         endif
         IF ((ANGLE(J).LE.ANG1).AND.(ANG1.LT.ANGLE(jj))) THEN
         if (angle(1).gt.pi) angle(1)=angle(1)-pi
            CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                       K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
            KONTROL=1
         ENDIF
         if (angle(1).gt.pi) angle(1)=angle(1)-pi
         IF (L.NE.M) THEN
            L=L+1
         ELSE
            L=1
         ENDIF
         IF (KONTROL.EQ.1) HALT=1
         IF ((HALT.EQ.1).AND.(KONTROL.EQ.0)) THEN
            if (halt2.eq.1) goto 101
            if (l.gt.1) then
               jflag=l-1
            else
               jflag=m
            endif
            GOTO 79
         ELSE
            IF (L.EQ.jj) THEN
               if (jj.eq.1) halt2=1
               GOTO 73
            ELSE
               GOTO 75
            ENDIF
         ENDIF
      ENDIF
C
C  The first switch has occurred. Now start looking for the next ones,
C  between the following angles.
C
79    DO 80 I=J+1,M-1
         L=JFLAG
 90      KONTROL=0
         IF ((ANGLE(L)+PI2).LT.PI) THEN
            ANG1=ANGLE(L)+PI2
         ELSE
            ANG1=ANGLE(L)-PI2
         ENDIF
         IF ((ANGLE(I).LE.ANG1).AND.(ANG1.LT.ANGLE(I+1))) THEN
            CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +                  ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
            KONTROL=1
         ENDIF
         IF (KONTROL.EQ.0) THEN
            JFLAG=L
         ELSE
            IF (L.NE.M) THEN
               L=L+1
            ELSE
               L=1
            ENDIF
            GOTO 90
         ENDIF
 80   CONTINUE
      L=JFLAG
C
C  Finally, look for necessary switches between the last angle and zero.
C
100   KONTROL=0
      IF ((ANGLE(L)+PI2).LT.PI) THEN
         ANG1=ANGLE(L)+PI2
      ELSE
         ANG1=ANGLE(L)-PI2
      ENDIF
      IF ((ANGLE(M).LE.ANG1).AND.(ANG1.LT.PI)) THEN
         CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +               ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
         KONTROL=1
      ENDIF
      IF (KONTROL.EQ.1) THEN
         IF (L.NE.M) THEN
             L=L+1
         ELSE
             L=1
         ENDIF
         GOTO 100
      ENDIF 
101      NUM=KOUNT-1
C  
C  Sort the NUM special k-dividers. 
C  Permute KAND1, KAND2 and D in the same way.
C

      if (dble(k).gt.dble(n/2)*factor) then
c ---------------
      CALL SORT3(ALPHA,KAND1,KAND2,D,NUM)
      IW1=1
      IW2=2
      JFULL=0
      NDK=0
      zoek=0
      zold=0
      firstiw1=0
      firstiw2=0
      nloop=0
      iw1old=0
      iw1oldloop=0
      iw2old=0
      moredith=0
120   NDATA=0
C
C  Compute the intersection point.
C
      IF (DABS(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +         +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))).LT.EPS) THEN
         IW2=IW2+1
         IF (IW2.EQ.NUM+1) IW2=1
         if ((zoek.ne.0).and.(iw2.eq.iw2nu)) then
            zoek=0
            zold=iw1
         endif
         GOTO 120
      ENDIF
      XCORD=(DCOS(ALPHA(IW2))*D(IW1)-DCOS(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +                   +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2)))
      YCORD=(-DSIN(ALPHA(IW2))*D(IW1)+DSIN(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))
     +                   +DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1)))
C 
C  Test whether the intersection point is a data point. 
C  If so, adjust IW1 and IW2.
C
      IF ((KAND1(IW1).EQ.KAND1(IW2)).OR.(KAND1(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND1(IW1)
      IF ((KAND2(IW1).EQ.KAND1(IW2)).OR.(KAND2(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND2(IW1)
      IF (NDATA.NE.0) THEN
        iv=0
125       NEXT=IW2+1
            iv=iv+1
          IF (NEXT.EQ.(NUM+1)) NEXT=1
          if (next.ne.iw1) then
          IF ((NDATA.EQ.KAND1(NEXT)).OR.(NDATA.EQ.KAND2(NEXT))) THEN
            IW2=IW2+1
            IF (IW2.EQ.(NUM+1)) IW2=1
         if ((zoek.ne.0).and.(iw2.eq.iw2nu)) then
            zoek=0
            zold=iw1
         endif
            GOTO 125
          ENDIF
          endif
         if (iv.eq.(num-1)) then
            num=1
                  KORNR(1,1)=KAND1(IW1)
                  KORNR(1,2)=KAND2(IW1)
                  KORNR(1,3)=KAND1(IW2)
                  KORNR(1,4)=KAND2(IW2)
            return
         endif
      ENDIF
      IF (IW2.EQ.NUM) THEN
         KON=1
      ELSE
         KON=IW2+1
      ENDIF
      if (kon.eq.iw1) kon=kon+1
      if (kon.eq.num+1) kon=1
C
C  Test whether the intersection point lies to the left of the special 
C  k-divider which corresponds to ALPHA(KON). If so, compute its depth.
C
      IF ((DSIN(ALPHA(KON))*XCORD-DCOS(ALPHA(KON))*YCORD
     +     -D(KON)).le.eps) THEN

      CALL TUKDEPTH(XCORD,YCORD,N,X,Y,BETA,F,DPF,JLV,JRV,HDEP1,sdep,eps)

c      hdep1=int(hd*dble(n)+.5)

      if ((ndk.eq.1).and.(iw1.gt.iw2).and.(iw2.lt.firstiw1).and.
     +    (hdep1.eq.k)) then
          iw2=firstiw1
          goto 120
      endif
      IF (HDEP1.eq.K) NDK=1
         if ((ndk.eq.1).and.(firstiw1.eq.0)) firstiw1=iw1

      IF (HDEP1.ne.K) THEN
      if ((hdep1.lt.k).and.(ndk.eq.1)) then
C
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
         IW2=IW2+1
         IF (IW2.EQ.(NUM+1)) IW2=1
         GOTO 120
      ENDIF

      if (ndk.eq.1) zoek=0
      if (iw1.ne.zold) then
         if ((zoek.eq.0).and.(ndk.eq.0).and.(hdep1.ge.khulp)) then
            zoek=iw1
            iw2nu=iw2
         endif
         if (zoek.ne.0) then
            iw2=iw2+1
            if (iw2.eq.num+1) iw2=1
            if (iw2.eq.iw1) iw2=iw2+1
            if (iw2.eq.num+1) iw2=1
            if (iw2.eq.iw2nu) then
               zoek=0
               zold=iw1
               goto 1111
            endif
            goto 120
         endif
      endif

      ENDIF
C
C  Store IW1 and IW2 in KORNR. If KORNR has already been filled, check whether 
C  we have encountered this intersection point before.
C
1111   continue
       if (ndk.eq.1) then
       if (iw1.ne.overschr) overschr=0
       if ((iw1.eq.iw2old).and.(iw2.eq.iw1old)) then
          moredith=1
          goto 170
       endif
       if ((iw1.gt.iw2).and.(iw2.gt.firstiw1)) then
          moredith=1
          goto 170
       endif
       endif
       IF ((IW2.GT.IW1).AND.(JFULL.EQ.0)) THEN
            DO 130 I=IW1,IW2-1
               KORNR(I,1)=KAND1(IW1)
               KORNR(I,2)=KAND2(IW1)
               KORNR(I,3)=KAND1(IW2)
               KORNR(I,4)=KAND2(IW2)
130         CONTINUE
         ELSE
            IF (IW2.GT.IW1) THEN
               DO 140 I=IW1,IW2-1
          IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +   .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +    THEN
                  GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
             ENDIF
140            CONTINUE
            ELSE
               JFULL=1
               DO 150 I=IW1,NUM
                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
150            CONTINUE
               DO 160 I=1,IW2-1
        IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +  .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +       THEN
                  GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
             ENDIF
160            CONTINUE
            ENDIF
         ENDIF
         iw1old=iw1
         iw2old=iw2
      ELSE
C
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
         IW2=IW2+1
         IF (IW2.EQ.(NUM+1)) IW2=1
         if ((zoek.ne.0).and.(iw2.eq.iw2nu)) then
             zoek=0
             zold=iw1
         endif
         if (iw1oldloop.eq.iw1) nloop=nloop+1
         if (iw1.ne.iw1oldloop) nloop=0
         iw1oldloop=iw1
         if (nloop.gt.num) then
            moredith=1
            goto 170
         endif
         GOTO 120
      ENDIF
C
C  Look for the next vertex of the convex figure.
C
      IW1=IW2
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 120
 

c ------------
      else
c ------------

      CALL SORT3(ALPHA,KAND1,KAND2,D,NUM)
      
      IW1=1
      IW2=2
      JFULL=0
      NDK=0
1200   NDATA=0
C
C  Compute the intersection point.
C
      IF (DABS(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +         +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))).LT.EPS) THEN
         IW2=IW2+1
         IF (IW2.EQ.NUM+1) IW2=1
         GOTO 1200
      ENDIF
      XCORD=(DCOS(ALPHA(IW2))*D(IW1)-DCOS(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +                   +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2)))
      YCORD=(-DSIN(ALPHA(IW2))*D(IW1)+DSIN(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))
     +                   +DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1)))
C 
C  Test whether the intersection point is a data point. 
C  If so, adjust IW1 and IW2.
C

      IF ((KAND1(IW1).EQ.KAND1(IW2)).OR.(KAND1(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND1(IW1)
      IF ((KAND2(IW1).EQ.KAND1(IW2)).OR.(KAND2(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND2(IW1)
      IF (NDATA.NE.0) THEN
        iv=0
1250       NEXT=IW2+1
            iv=iv+1
          IF (NEXT.EQ.(NUM+1)) NEXT=1
          if (next.ne.iw1) then
          IF ((NDATA.EQ.KAND1(NEXT)).OR.(NDATA.EQ.KAND2(NEXT))) THEN
            IW2=IW2+1
            IF (IW2.EQ.(NUM+1)) IW2=1
            GOTO 1250
          ENDIF
          endif
          if(iv.eq.(num-1)) then
            num=1
                  KORNR(1,1)=KAND1(IW1)
                  KORNR(1,2)=KAND2(IW1)
                  KORNR(1,3)=KAND1(IW2)
                  KORNR(1,4)=KAND2(IW2)
            return
          endif
      ENDIF
      IF (IW2.EQ.NUM) THEN
         KON=1
      ELSE
         KON=IW2+1
      ENDIF
       if (kon.eq.iw1) kon=kon+1
       if (kon.eq.num+1) kon=1

C
C  Test whether the intersection point lies to the left of the special 
C  k-divider which corresponds to ALPHA(KON). If so, compute its depth.
C
      IF ((DSIN(ALPHA(KON))*XCORD-DCOS(ALPHA(KON))*YCORD
     +     -D(KON)).LE.eps) THEN
      CALL TUKDEPTH(XCORD,YCORD,N,X,Y,BETA,F,DPF,JLV,JRV,HDEP1,sdep,eps)
c      hdep1=int(hd*dble(n)+.5)


      IF (HDEP1.EQ.K) NDK=1
      IF (HDEP1.NE.K) THEN
      CALL TUKDEPTH(XCORD-EPS*2,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +  JLV,JRV,HDEP2,sdep,eps)
c      hdep2=int(hd*dble(n)+.5)
      CALL TUKDEPTH(XCORD+EPS*2,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +  JLV,JRV,HDEP3,sdep,eps)
c      hdep3=int(hd*dble(n)+.5)
      CALL TUKDEPTH(XCORD-EPS*2,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +  JLV,JRV,HDEP4,sdep,eps)
c      hdep4=int(hd*dble(n)+.5)
      CALL TUKDEPTH(XCORD+EPS*2,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +  JLV,JRV,HDEP5,sdep,eps)
c      hdep5=int(hd*dble(n)+.5)

      IF ((NDK.EQ.0).AND.
     +    ((HDEP1.ge.K).OR.(HDEP2.ge.K).OR.(HDEP3.ge.K)
     +      .OR.(HDEP4.ge.K).OR.(HDEP5.ge.K))) THEN 
      NDK=1
      ENDIF
      IF ((HDEP1.LT.K).AND.(HDEP2.LT.K)
     +   .AND.(HDEP3.LT.K).AND.(HDEP4.LT.K)
     +   .AND.(HDEP5.LT.K).AND.(NDK.EQ.1)) THEN
C
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
         IW2=IW2+1
         IF (IW2.EQ.(NUM+1)) IW2=1
         GOTO 1200
      ENDIF
      ENDIF
C
C  Store IW1 and IW2 in KORNR. If KORNR has already been filled, check whether 
C  we have encountered this intersection point before.
C
         IF ((IW2.GT.IW1).AND.(JFULL.EQ.0)) THEN
            DO 1300 I=IW1,IW2-1
               KORNR(I,1)=KAND1(IW1)
               KORNR(I,2)=KAND2(IW1)
               KORNR(I,3)=KAND1(IW2)
               KORNR(I,4)=KAND2(IW2)
1300         CONTINUE
         ELSE
            IF (IW2.GT.IW1) THEN
               DO 1400 I=IW1,IW2-1
          IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +   .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +    THEN
                  GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
             ENDIF
1400            CONTINUE
            ELSE
               JFULL=1
               DO 1500 I=IW1,NUM
                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
1500            CONTINUE
               DO 1600 I=1,IW2-1
       IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +  .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +       THEN
                  GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

                  KORNR(I,1)=KAND1(IW1)
                  KORNR(I,2)=KAND2(IW1)
                  KORNR(I,3)=KAND1(IW2)
                  KORNR(I,4)=KAND2(IW2)
             ENDIF
1600            CONTINUE
            ENDIF
         ENDIF
      ELSE
C
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
         IW2=IW2+1
         IF (IW2.EQ.(NUM+1)) IW2=1
         GOTO 1200
      ENDIF
C
C  Look for the next vertex of the convex figure.
C 
      IW1=IW2
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 1200

c ------------
      endif
170   if (ndk.eq.0) empty=1
      if(dble(k).le.dble(n/2)*factor.and.empty.eq.1) moredith=1
      RETURN
      END

c ---------------------------------------------------------------------

      SUBROUTINE SORT3(B,I1,I2,R,N)
C
C  Sorts a real array B of length N and permutes two integer arrays 
C  I1 and I2 and one real array R in the same way.
C
      INTEGER i,N,q(n),I1(n),I2(n),II1(n),II2(n)
      REAL*8 R(n),rr(n),b(n),x(n)

      call indexx(n,b,q)
      
      do 10 i=1,n
        x(i)=b(i)
        II1(i)=I1(i)
        II2(i)=I2(i)
        rr(i)=r(i)
10    continue
 
      do 20 i=1,n
        b(i)=x(q(i))
        I1(i)=II1(q(i))
        I2(i)=II2(q(i))
        r(i)=rr(q(i))
20    continue

      END

c ---------------------------------------------------------------------

      SUBROUTINE ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,
     +           ALPHA,ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
C
C  Updates NCIRQ and NRANK, detects the special k-dividers and stores 
C  their angles and the constant terms of their equations.
C
      integer m,maxnum
      INTEGER NCIRQ(N),NRANK(N),IND1(M),IND2(M)
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM)
      INTEGER KOUNT,K,L,N,IV,IV1,IV2,D1,D2
      REAL*8 X(N),Y(N),ANGLE(M),D(M)
      REAL*8 ALPHA(MAXNUM),DUM,PI,PI2
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      D1=IND1(L)
      IV1=NRANK(D1)
      D2=IND2(L)
      IV2=NRANK(D2)
      IV=NCIRQ(IV1)
      NCIRQ(IV1)=NCIRQ(IV2)
      NCIRQ(IV2)=IV
      IV=IV1
      NRANK(D1)=IV2
      NRANK(D2)=IV
         IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +      .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))
     +      .OR.((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1))) 
     +      .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
            IF (ANGLE(L).LT.PI2) THEN
               DUM=ANGLE(L)+PI2
            ELSE
               DUM=ANGLE(L)-PI2
            ENDIF
            IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +         .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))) THEN
               IF (DUM.LE.PI2) THEN
                  ALPHA(KOUNT)=ANGLE(L)+PI
               ELSE
                  ALPHA(KOUNT)=ANGLE(L)
               ENDIF
            ENDIF
            IF (((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1)))
     +        .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
               IF (DUM.LE.PI2) THEN
                  ALPHA(KOUNT)=ANGLE(L)
               ELSE
                  ALPHA(KOUNT)=ANGLE(L)+PI
               ENDIF
            ENDIF
            KAND1(KOUNT)=IND1(L)
            KAND2(KOUNT)=IND2(L)
            D(KOUNT)=DSIN(ALPHA(KOUNT))*X(IND1(L))
     +                -DCOS(ALPHA(KOUNT))*Y(IND1(L))
            KOUNT=KOUNT+1
         ENDIF
      RETURN
      END
      
c ---------------------------------------------------------------------

      SUBROUTINE tukDEPTH(U,V,N,X,Y,BETA,F,DPF,JLV,JRV,HDEP,sdep,epsi)
C
C  Computes the halfspace depth of a point. This subroutine was described
C  in: Rousseeuw, P.J. and Ruts, I. (1996). AS 307:  Bivariate location 
C  depth, Applied Statistics (JRSS-C) 45, 516-526.
C
      integer n,jlv,jrv,numh,nz,i,nn,nu,ja,jb,nn2,nbad,nf,j,k,ki
      REAL*8 U,V,BETA(N),X(N),Y(N),DPF(N)
      REAL*8 P,P2,EPSI,D,XU,YU,ANGLE,ALPHK,BETAK
      INTEGER F(N),GI,HDEP,sdep
      DIMENSION JLV(N),JRV(N)
      NUMH=0
      HDEP=0
      IF (N.LT.1) RETURN
      P=DACOS(DBLE(-1.0))
      P2=P*2.0
      
      NZ=0
C
C  Construct the array BETA.
C

      DO 10 I=1,N
          D=DSQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPSI) THEN
              NZ=NZ+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (DABS(XU).GT.DABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      BETA(I-NZ)=DASIN(YU)
                      IF(BETA(I-NZ).LT.0.0) THEN
                          BETA(I-NZ)=P2+BETA(I-NZ)
                      ENDIF
                  ELSE
                      BETA(I-NZ)=P-DASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      BETA(I-NZ)=DACOS(XU)
                  ELSE
                      BETA(I-NZ)=P2-DACOS(XU)
                  ENDIF
              ENDIF
              IF (BETA(I-NZ).GE.(P2-EPSI)) BETA(I-NZ)=0.0
          ENDIF
  10  CONTINUE
      NN=N-NZ
      IF (NN.LE.1) GOTO 60
C
C  Sort the array BETA.
C
      DO 15 I=1,NN
      DPF(I)=DBLE(F(I))
15    CONTINUE
      CALL SORT3(BETA,F,F,DPF,NN)
C
C  Check whether Z=(U,V) lies outside the data cloud.
C
      ANGLE=BETA(1)-BETA(NN)+P2
      DO 20 I=2,NN
          ANGLE=DMAX1(ANGLE,(BETA(I)-BETA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPSI)) GOTO 60
C
C  Make smallest BETA equal to zero,
C  and compute NU = number of BETA < PI.
C
      ANGLE=BETA(1)
      NU=0
      DO 30 I=1,NN
          BETA(I)=BETA(I)-ANGLE
          IF (BETA(I).LT.(P-EPSI)) NU=NU+1
  30  CONTINUE
      IF (NU.GE.NN) GOTO 60
C
C  Mergesort the BETA with their antipodal angles,
C  and at the same time update I, F(I), and NBAD.
C
      JA=1
      JB=1
      ALPHK=BETA(1)
      BETAK=BETA(NU+1)-P
      NN2=NN*2
      NBAD=0
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPSI).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=BETA(JA)
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
                      BETAK=BETA(JB+NU)-P
                  ELSE
                      BETAK=BETA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.0
              ENDIF
          ENDIF
  40  CONTINUE
C
C  Computation of NUMH for halfspace depth.
C
      GI=0
      JA=1
      ANGLE=BETA(1)
      NUMH=MIN0(F(1),(NN-F(1)))
      DO 50 I=2,NN
          IF(BETA(I).LE.(ANGLE+EPSI)) THEN
              JA=JA+1
          ELSE
              GI=GI+JA
              JA=1
              ANGLE=BETA(I)
          ENDIF
          KI=F(I)-GI
          NUMH=MIN0(NUMH,MIN0(KI,(NN-KI)))
   50 CONTINUE
C
C  Adjust for the number NZ of data points equal to Z=(U,V):
C
   60 NUMH=NUMH+NZ
      HDEP=NUMH
      
      IF (N.GE.3) SDEP=(NUMS+0.0)/(K(N,3)+0.0)
      RETURN
      END


c ====================================================================
c                        deeploc.f
c ====================================================================

      subroutine deeplocstand(maxn,maxp,n,np,x,xn,eps,locsca,err)
      implicit double precision (a-h,o-z)
      integer maxn,maxp,n,np,err(np),jn,j,i
      double precision x(maxn,maxp),xn(n),eps
      double precision qloc,qsca,ave,var,locsca(maxp,2)

      jn=0
      do 10 j=1,np
         do 20 i=1,n
            xn(i)=x(i,j)
 20      continue
         if ((2*int(n/2)).eq.n) then
            qloc=findq(xn,n,n/2)
            qloc=(findq(xn,n,(n/2)+1)+qloc)/2.d0
         else
            qloc=findq(xn,n,int(n/2)+1)
         endif
         do 30 i=1,n
            xn(i)=dabs(x(i,j)-qloc)
 30      continue
         if ((2*int(n/2)).eq.n) then
            qsca=findq(xn,n,n/2)
            qsca=(findq(xn,n,(n/2)+1)+qsca)/2.d0
         else
            qsca=findq(xn,n,int(n/2)+1)
         endif
         if (dabs(qsca).lt.eps) then
            ave=0.d0
            do 40 i=1,n
               ave=ave+x(i,j)
 40         continue
            ave=ave/(n+0.d0)
            var=0.d0
            do 50 i=1,n
               var=var+(x(i,j)-ave)*(x(i,j)-ave)
 50         continue  
            if (n.ne.1) var=var/(n-1.d0)
            if (dabs(var).lt.eps) then
               if (np.ne.1) then
                  np=np-1
                  err(j)=-1
                  goto 10
               endif
            else
               err(j)=-2
               qsca=dsqrt(var)
            endif
         endif
         jn=jn+1
         locsca(jn,1)=qloc
         locsca(jn,2)=qsca
         do 60 i=1,n
            x(i,jn)=(x(i,j)-qloc)/qsca
 60      continue         
 10   continue

      return
      end

c ---------------------------------------------------------------------

      subroutine deepest(maxn,maxp,n,np,ndir,x,eps,nddpst,dpstM,
     +     nstp,ntry,nalt,err,errc,n4)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  Author: Anja Struyf
CC          Department of Mathematics and Computer Science
CC          University of Antwerp (UIA)
CC          Universiteitsplein 1
CC          B-2610 Wilrijk-Antwerpen
CC          Belgium
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  This program computes approximates the deepest location
CC  relative to given a data set, using the halfspace location depth.
CC  The accompanying paper "High-dimensional computation of the deepest 
CC  location", by Anja Struyf and Peter J. Rousseeuw can be obtained 
CC  from our website
CC            http://win-www.uia.ac.be/u/statis/index.html
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      implicit double precision (a-h,o-z)
      integer maxn,maxp,n,np,ndir,n4,err(n4),errc
      integer nddpst,nstp,ntry,nalt,nrun,j,l,numh,nt,mindep,nceil2,i
      integer nsa,ngen,ind,jdir,nmin,nmax,nsin,ji,nsamp,ierr
      integer nid,nidalt,ndstep,ndold,nstep,inalt,nrankl,nrankg
      integer lj,ntave,nback,mold,indM,indMold

      double precision x(maxn,maxp),eps,stepsM(2*(np+2),np),xn(n)
      double precision cov(np,np),ave(np),evecs(np,np),evals(np),d1(np)
      double precision utx(n,ndir),utxsort(n,ndir),dpstM(np),u(ndir,np)
      integer jsamp(np)

cc  initialize the random seed.
      nrun=0

cc  handle special case where n is equal to 1.
      if (n.eq.1)then
         do 1 j=1,np
            dpstM(j)=x(1,j)
 1       continue
         nddpst=1
         return
      endif

cc  handle special case where np is equal to 1.
      if (np.eq.1) then
         do 2 l=1,n
            xn(l)=x(l,1)
 2       continue
         dpstM(1)=dpmedian(xn,n)
         numh=0
         nt=0
         do 3 l=1,n
            if (x(l,1).gt.(dpstM(1)+eps)) then
               numh=numh+1
            elseif (x(l,1).ge.(dpstM(1)-eps)) then
               nt=nt+1
            endif
 3       continue
         nddpst=min0(numh+nt,n-numh)
         return
      endif

cc  general case.

cc  initialize the minimal depth of the deepest location.
      mindep=nceil2(dble(n/(np+1)),eps)
      
cc compute the coordinate-wise median, used as first approximation 
cc for the deepest location.
      do 4 j=1,np
         do 5 i=1,n
            xn(i)=x(i,j)
 5       continue
         dpstM(j)=dpmedian(xn,n)
 4    continue

cc  construct ndir unit vectors: these are the directions used to compute
cc  depths and moving directions.
cc  1) add coordinate axes
      do 6 l=1,np
         do 7 j=1,np
            if (j.eq.l) then
               u(l,j)=1.d0
            else
               u(l,j)=0.d0
            endif
 7       continue
 6    continue

cc  2) add vectors connecting data points with the coordinate-wise median
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      if (n.le.int(ndir/4)) then
         nsa=n
         ngen=0
      else
         nsa=int(ndir/4)
         ngen=1
      endif
      ind=np

      do 8 jdir=1,nsa
         if (ngen.eq.1) then
            call randm(nrun,ran)
            i=int(n*ran+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
         else
            jsamp(1)=jdir
         endif
         utj=0.d0
         do 9 j=1,np
            u(ind+1,j)=x(jsamp(1),j)-dpstM(j)
            utj=utj+u(ind+1,j)*u(ind+1,j)
 9       continue
         utj=dsqrt(utj)
         if (utj.gt.eps) then
            ind=ind+1
            do 10 j=1,np
               u(ind,j)=u(ind,j)/utj
 10         continue
         endif
 8    continue

cc  add vectors connecting two data points
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      if ((n*(n-1)/2).le.(int(ndir/2)-ind)) then
         nsa=(n*(n-1))/2
         ngen=0
      else
         nsa=int(ndir/2)-ind
         ngen=1
      endif
      nmin=1
      nmax=2

      do 11 jdir=1,nsa
         if (ngen.eq.1) then
            call randm(nrun,ran)
            i=int(n*ran+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
 12         call randm(nrun,ran)
            l=int(n*ran+1.d0)
            if(l.gt.n)l=n
            if(l.eq.jsamp(1)) goto 12
            jsamp(2)=l
         else
            jsamp(1)=nmin
            jsamp(2)=nmax
            if (nmax.eq.n) then
               nmin=nmin+1
               nmax=nmin+1
            else
               nmax=nmax+1
            endif
         endif
         utj=0.d0
         do 13 j=1,np
            u(ind+1,j)=x(jsamp(1),j)-x(jsamp(2),j)
            utj=utj+u(ind+1,j)*u(ind+1,j)
 13      continue
         utj=dsqrt(utj)
         if (utj.gt.eps) then
            ind=ind+1
            do 14 j=1,np
               u(ind,j)=u(ind,j)/utj
 14        continue
         endif
 11   continue

cc  add vectors perpendicular to hyperplanes through np data points
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      nsin=0
      dmax=1.d0
      if (np.gt.int(n/2)) then 
         p=n-np
      else
         p=np
      endif
      j=p
 15   dmax=dmax*dble(n-p+j)/dble(p-j+1)
      if (dmax.gt.dble(ndir-ind)) then
         nsa=ndir-ind
         ngen=1
         goto 21
      endif
      j=j-1
      if (j.ge.1) goto 15
      nsa=int(dmax)
      ngen=0

 21   do 100 jdir=1,nsa
         if (ngen.eq.0) then
            if (jdir.eq.1) then
               do 16 j=1,np
                  jsamp(j)=j
 16            continue
            else
               j=np
 17            if (jsamp(j).lt.(n-np+j)) then
                  jsamp(j)=jsamp(j)+1
                  do 18 ji=(j+1),np
                     jsamp(ji)=jsamp(ji-1)+1
 18               continue
               else
                  j=j-1
                  goto 17
               endif
            endif
         else
            call randm(nrun,ran)
            i=int(n*ran+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
            nsamp=1
 19         call randm(nrun,ran)
            l=int(n*ran+1.d0)
            if(l.gt.n)l=n
            do 20 j=1,nsamp
               if(l.eq.jsamp(j)) goto 19
 20         continue
            nsamp=nsamp+1
            jsamp(nsamp)=l
            if (nsamp.lt.np)goto 19
         endif
cc  compute the covariance matrix of the sample.
         do 30 j=1,np
            ave(j)=0.d0
            do 40 i=1,np
               ave(j)=ave(j)+x(jsamp(i),j)
 40         continue
            ave(j)=ave(j)/np
 30      continue
         do 50 j=1,np
            do 60 l=1,j
               cov(j,l)=0.d0
               do 70 i=1,np
                  cov(j,l)=cov(j,l)+(x(jsamp(i),j)-ave(j))
     +                 *(x(jsamp(i),l)-ave(l))
 70            continue
               cov(j,l)=cov(j,l)/(np-1)
               cov(l,j)=cov(j,l)
 60         continue
 50      continue
cc  compute the eigenvalues and corresponding eigenvectors 
cc  of the covariance matrix.
         call eigen(np,np,cov,evals,evecs,d1,ierr)
         if (ierr.ne.0) then
            err(jdir)=ierr
            nsin=nsin+1
            goto 100
         endif
         if (evals(1).gt.eps) then
            err(jdir)=-1
            nsin=nsin+1
            goto 100
         endif
cc  test for singularity of the sample.
         if (evals(2).le.eps) then
            nsin=nsin+1
         endif
cc  determine the direction orthogonal to the sample.
         nt=0
         do 80 j=1,np
            if (dabs(evecs(j,1)).le.eps) nt=nt+1
 80      continue
         if (nt.eq.np) then
            err(jdir)=-2
            nsin=nsin+1
            goto 100
         endif
         ind=ind+1
         do 90 j=1,np
            u(ind,j)=evecs(j,1)
 90      continue
 100  continue
      
cc search the deepest location.
cc initialize.
      ndir=ind
      do 101 jdir=1,ndir
         do 102 i=1,n
            utx(i,jdir)=0.d0
            do 103 j=1,np
               utx(i,jdir)=utx(i,jdir)+u(jdir,j)*x(i,j)
 103        continue
            xn(i)=utx(i,jdir)
 102        continue
         call sort(xn,n)
         do 104 i=1,n
            utxsort(i,jdir)=xn(i)
 104     continue
 101  continue
      do 105 j=1,np
         stepsM(1,j)=dpstM(j)
 105  continue
      indM=1
      nid=0
      nidalt=0
      ndstep=-1
      ndold=-1
      nt=0
      nstep=0
      inalt=2
     
cc Start iterations.
 110  ndold=max0(ndstep,ndold)
      nstep=nstep+1
      ndstep=n+1
      do 120 l=1,ndir
cc Compute rank of new deepest in direction l.
         utj=0.d0
         do 130 j=1,np
            utj=utj+u(l,j)*stepsM(indM,j)
 130     continue
         do 140 i=1,n
            xn(i)=utxsort(i,l)
 140     continue
         call irank(utj,xn,n,eps,nrankl,nrankg)
         if (nrankg.lt.ndstep) then
            ndstep=nrankg
            do 150 j=1,np
               ave(j)=u(l,j)
 150        continue
            lj=l
            nt=1
         elseif (nrankg.eq.ndstep) then
            do 160 j=1,np
               ave(j)=ave(j)+u(l,j)
 160        continue            
            nt=nt+1
         endif
         if (nrankl.lt.ndstep) then
            ndstep=nrankl
            do 170 j=1,np
               ave(j)=-u(l,j)
 170        continue
            lj=-l
            nt=1
         elseif (nrankl.eq.ndstep) then
            do 180 j=1,np
               ave(j)=ave(j)-u(l,j)
 180        continue            
            nt=nt+1
         endif
 120  continue
      ntave=0
      do 190 j=1,np
         ave(j)=ave(j)/(nt+0.d0)
         if (dabs(ave(j)).lt.eps) ntave=ntave+1
 190  continue            
      if (ntave .eq. np) then
         do 195 j=1,np
            if (lj.gt.0) ave(j)=u(lj,j)
            if (lj.lt.0) ave(j)=-u(lj,j)
 195     continue
      endif

      if (ndstep.ge.ndold) then
cc New deepest point found.
         nddpst=ndstep
         do 200 j=1,np
            dpstM(j)=stepsM(indM,j)
 200     continue
         if (inalt.ne.2) inalt=2

         if (ndold.eq.ndstep) then
            nid=nid+1
            nidalt=nidalt+1
         else 
            nid=0
            nidalt=0
         endif

      elseif (inalt.ne.2) then 
cc Alternative directions loop is active, but a bad direction was tried, 
cc hence try another direction.
         nid=nid+1
         goto 220
      else
         nidalt=nidalt+1
         nid=nid+1
      endif

      if (nidalt.ge.nalt) then
cc Save present configuration, and enter alternative directions loop.
         nback=ndstep
         do 210 j=1,np
            evals(j)=ave(j)
            evecs(j,1)=stepsM(indM,j)
 210     continue
         goto 220
      endif
      
      goto 300

cc Try some alternative directions, connecting the latest approximation
cc of the deepest location with one of the earlier approximations.
 220  nidalt=0
      if (inalt.eq.2*(np+2)) then
cc Restore old configuration.
         ndstep=nback
         if (indM.eq.2*(np+2)) then
            indM=1
         else
            indM=indM+1
         endif
         do 230 j=1,np
            ave(j)=evals(j)
            stepsM(indM,j)=evecs(j,1)
 230     continue            
         if (ndold.ge.ndstep) nid=nid+1
         inalt=2
         goto 300
      endif
      
      if (indM.eq.1) then
         indM=2*(np+2)
      else
         indM=indM-1
      endif
      if ((indM+inalt).gt.2*(np+2)) then
         mold=indM-2*(np+2)+inalt
      else
         mold=indM+inalt
      endif
      do 245 j=1,np
         ave(j)=stepsM(mold,j)-stepsM(indM,j)
 245  continue
      inalt=inalt+1
      goto 400
      
 300  if (nid.ge.ntry) then
         errc=nstep
         return
      endif

cc Take a step in the computed direction.
 400  if (nstep.ge.nstp) then
         errc=-nstep
         return
      endif
      utdpst=0.d0
      do 240 j=1,np
         utdpst=utdpst+ave(j)*ave(j)
 240  continue
      utdpst=dsqrt(utdpst)
      do 250 j=1,np
         ave(j)=ave(j)/utdpst
 250  continue
      utdpst=0.d0
      do 410 j=1,np
         utdpst=utdpst+ave(j)*stepsM(indM,j)
 410  continue
      do 420 i=1,n
         xn(i)=0.d0
         do 430 j=1,np
            xn(i)=xn(i)+ave(j)*x(i,j)
 430     continue
 420  continue
      call sort(xn,n)
      call irank(utdpst,xn,n,eps,nrankl,nrankg)
      
      if (ndstep.ge.int(n/2)) then
         return
      elseif (nrankg.lt.mindep) then
         utdpst=-(utdpst-xn(n+1-mindep))
      else
         utdpst=-(utdpst-xn(n-nrankg))
      endif
      indMold=indM
      if (indM.eq.2*(np+2)) then
         indM=1
      else
         indM=indM+1
      endif   
      do 440 j=1,np
         stepsM(indM,j)=stepsM(indMold,j)+utdpst*ave(j)
 440  continue
      goto 110

      end

c -------------------------------------------------------------------

      function dpmedian(aw,ncas)
cc  Finds the median of the array aw of length ncas.
      implicit double precision (a-h,o-z)
      integer ncas
      double precision dpmedian
      double precision aw(ncas),qloc

      if ((2*int(ncas/2)).eq.ncas) then
         qloc=findq(aw,ncas,ncas/2)
         qloc=(findq(aw,ncas,(ncas/2)+1)+qloc)/2.d0
      else
         qloc=findq(aw,ncas,int(ncas/2)+1)
      endif
      dpmedian=qloc
      end

c --------------------------------------------------------------------

      subroutine irank(u,aw,n,eps,indle,indge)
      implicit double precision (a-h,o-z)
      integer n,indle,indge,indl,indg,imin,imax,j
      double precision aw(n),u,eps

      if (u.lt.(aw(1)-eps)) then
         indge=n
         indle=0
         return
      elseif (u.le.(aw(1)+eps)) then
         indge=n
         indle=1
         indl=1
         goto 200
      endif
      if (u.gt.(aw(n)+eps)) then
         indge=0
         indle=n
         return
      elseif (u.ge.(aw(n)-eps)) then
         indge=1
         indle=n
         indg=n
         goto 100
      endif
      imin=1
      imax=n
 50   if ((imax-imin).eq.1) then
         indge=n-imin
         indle=imin
         return
      endif
      j=int((imax+imin)/2)
      if (u.lt.(aw(j)-eps)) then
         imax=j
      elseif (u.gt.(aw(j)+eps)) then
         imin=j
      else
         indge=n-j+1
         indle=j
         indl=j
         indg=j
         goto 100         
      endif
      goto 50
 100  if (dabs(aw(indg-1)-u).le.eps) then
         indge=indge+1
         indg=indg-1
         goto 100
      endif
      if (indle.eq.n) return
 200  if (dabs(aw(indl+1)-u).le.eps) then
         indle=indle+1
         indl=indl+1
         goto 200
      endif
      end

c ---------------------------------------------------------------------

      function nceil2(u,eps)
      implicit double precision (a-h,o-z)
      integer nceil2

      nceil2=int(u)
      if (dabs(dble(nceil2-u)).gt.eps) nceil2=nceil2+1
      end

c ---------------------------------------------------------------------


