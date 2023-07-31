C********************************TOP***************************
C     UMAT FOR ABAQUS/STANDARD
C     A ISOTROPIC ELASTO-VISCOPLASTIC MODEL FOR 316 STAINLESS STEEL
C     THIS CODE WAS DEVELOPED BY XIAOWEI YUE
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      
      PARAMETER (M=3,N=3,ID=3,ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 SIX=6.D0, NINE=9.D0, TOLER=1.D-5)
      
      DIMENSION XIDEN(M,N),XNV(6),DPSTRAN(6), STRESSOLD(6),
     1DESTRAN(6), DSTRESS(6), XNDIR(M,N),ESTR(M,N),
     2STR(M,N),DSTR(M,N),EDSTR(M,N),DPSTRN(M,N),
     3BSTRESS1(6),BSTRESS2(6),BSTRESS(6),BSTRESSOLD(6),ESTRESS(6),
     4BSTR1(M,N),BSTR2(M,N),BSTR(M,N),BSTR1T(M,N),BSTR2T(M,N)
      

      DDSDDE = DDSDDE * 0.0D0
      
      E = PROPS(1)
      XNUE = PROPS(2)
      YIELD = PROPS(3)
      C1 = PROPS(4)
      C2 = PROPS(5)
      a1 = PROPS(6)
      a2 = PROPS(7)
      b = PROPS(8)
      Q = PROPS(9)
      KK = PROPS(10)
      vn = PROPS(11)
       
      DO I=1,M
          DO J=1,N
              STR(I,J)=0.0
              BSTR(I,J)=0.0
              BSTR1(I,J)=0.0
              BSTR2(I,J)=0.0
              ESTR(I,J)=0.0
          END DO
      END DO

      
      DO K=1,6
          BSTRESS1(K)=STATEV(K)
      END DO
      

      DO K=1,6
          BSTRESS2(K)=STATEV(K+6)
      END DO
      

      
      p = STATEV(13)
      R = STATEV(14)
      

      
      EBULK3 = E/(ONE-TWO*XNUE)
      XK = EBULK3/THREE
      EG2 = E/(ONE+XNUE)
      EG = EG2/TWO
      ELAM = (EBULK3-EG2)/THREE
      
      DO K1 = 1,3
          DO K2 = 1,3
              DDSDDE(K2,K1) = ELAM
          END DO
          DDSDDE(K1,K1) = EG2+ELAM
      END DO
      
      DDSDDE(4,4)=EG
      DDSDDE(5,5)=EG
      DDSDDE(6,6)=EG
      
      
      DO 50 I=1,M
          DO 50 J=1,N
              IF(I .EQ. J)THEN
                  XIDEN(I,J)=1.0D0
              ELSE
                  XIDEN(I,J)=0.0D0
              END IF
50    CONTINUE

      
      DO K=1,6
          STRESSOLD(K)=STRESS(K)
          BSTRESS(K)=BSTRESS1(K)+BSTRESS2(K)
          BSTRESSOLD(K)=BSTRESS(K)
      END DO
      
      
      CALL KMLT1(DDSDDE,DSTRAN,DSTRESS,NTENS)
      DO K=1,NTENS
          STRESS(K)=STRESS(K)+DSTRESS(K)
      END DO

      
      DO K=1,NTENS
          ESTRESS(K)=STRESS(K)-BSTRESS(K)
      END DO
      

      
      DO K = 1,3
          ESTR(K,K)=ESTRESS(K)
      END DO
      ESTR(1,2)=ESTRESS(4)
      ESTR(2,1)=ESTRESS(4)
      ESTR(2,3)=ESTRESS(5)
      ESTR(3,2)=ESTRESS(5)
      ESTR(3,1)=ESTRESS(6)
      ESTR(1,3)=ESTRESS(6)
      
      
      DO K = 1,3
          BSTR1(K,K)=BSTRESS1(K)
      END DO
      BSTR1(1,2)=BSTRESS1(4)
      BSTR1(2,1)=BSTRESS1(4)
      BSTR1(2,3)=BSTRESS1(5)
      BSTR1(3,2)=BSTRESS1(5)
      BSTR1(3,1)=BSTRESS1(6)
      BSTR1(1,3)=BSTRESS1(6)
      
      
      DO K = 1,3
          BSTR2(K,K)=BSTRESS2(K)
      END DO
      BSTR2(1,2)=BSTRESS2(4)
      BSTR2(2,1)=BSTRESS2(4)
      BSTR2(2,3)=BSTRESS2(5)
      BSTR2(3,2)=BSTRESS2(5)
      BSTR2(3,1)=BSTRESS2(6)
      BSTR2(1,3)=BSTRESS2(6)
      
      CALL KDEVIA(ESTR,XIDEN,EDSTR)  
      
      CALL KEFFP(EDSTR,EPJ)
      
      DO I=1,3
          DO J=1,3
              XNDIR(I,J)=EDSTR(I,J)/EPJ
          END DO
      END DO
      
      DO K=1,3
          XNV(K)=XNDIR(K,K)
      END DO
      XNV(4)=XNDIR(1,2)
      XNV(5)=XNDIR(3,2)
      XNV(6)=XNDIR(1,3)
      
      
      FLOW=EPJ-R-YIELD
      
      XDP=0.
      IF(FLOW.GT.0.)THEN
                
          
      R0=R
      DO I=1,3
          DO J=1,3
              BSTR1T(I,J)=BSTR1(I,J)
              BSTR2T(I,J)=BSTR2(I,J)
          END DO
      END DO
      DO K=1,10
          XFLOW=(EPJ-(THREE*EG*XDP)-R0-YIELD)
          XPHI=(XFLOW/KK)**vn
          XPHIDP=(-(vn*THREE*EG)/KK)*((XFLOW/KK)**(vn-ONE))
          XPHIEPJ=(vn/XK)*((XFLOW/KK)**(vn-ONE))
          XPHIR=-XPHIEPJ
          XRES=XPHI-(XDP/DTIME)

          XNDIRMBSTR1 = 0.
          DO I=1,3
              DO J=1,3
                  XNDIRMBSTR1=XNDIRMBSTR1+XNDIR(I,J)*BSTR1(I,J)*
     1            (THREE/TWO)
              END DO
          END DO
          XNDIRMBSTR2 = 0.
          DO I=1,3
              DO J=1,3
                  XNDIRMBSTR2=XNDIRMBSTR2+XNDIR(I,J)*BSTR2(I,J)*
     1            (THREE/TWO)
              END DO
          END DO
          DEQPL=XRES/((ONE/DTIME)-XPHIDP+XPHIEPJ*
     1    (C1*a1+C2*a2)-XPHIEPJ*C1*XNDIRMBSTR1-XPHIEPJ*C2*XNDIRMBSTR2-
     2    XPHIR*b*(Q-R))
          XDP=XDP+DEQPL
          R=R0+b*(Q-R0)*XDP
          DO I=1,3
              DO J=1,3
                  BSTR1(I,J)=BSTR1T(I,J)+(C1*a1*(TWO/THREE))*
     1            XDP*XNDIR(I,J)*(THREE/TWO)-C1*XDP*BSTR1T(I,J)
              END DO
          END DO
          DO I=1,3
              DO J=1,3
                  BSTR2(I,J)=BSTR2T(I,J)+(C2*a2*(TWO/THREE))*
     1            XDP*XNDIR(I,J)*(THREE/TWO)-C2*XDP*BSTR2T(I,J)
              END DO
          END DO
          
          
          
          IF(ABS(XRES).LT.1.E-12)GOTO 10
      END DO
10    CONTINUE
                 
          
      DO I=1,3
          DO J=1,3
              DPSTRN(I,J)=(THREE/TWO)*XDP*EDSTR(I,J)/EPJ
          END DO
      END DO
      
      
      DO K=1,3
          DPSTRAN(K)=DPSTRN(K,K)
      END DO
      DPSTRAN(4) = 2*DPSTRN(1,2)
      DPSTRAN(5) = 2*DPSTRN(3,2)
      DPSTRAN(6) = 2*DPSTRN(1,3)

      
      DO K=1,6
          DESTRAN(K)=DSTRAN(K)-DPSTRAN(K)
      END DO

      
      CALL KMLT1(DDSDDE,DESTRAN,DSTRESS,NTENS)
      
      
      DO K=1,6
          STRESS(K)=STRESSOLD(K)+DSTRESS(K)
      END DO
      
      p=p+XDP
      
      
      DO K=1,3
          STATEV(K)=BSTR1(K,K)
      END DO
      STATEV(4)=BSTR1(1,2)
      STATEV(5)=BSTR1(3,2)
      STATEV(6)=BSTR1(1,3)

      
      DO K=1,3
          STATEV(K+6)=BSTR2(K,K)
      END DO
      STATEV(4+6)=BSTR2(1,2)
      STATEV(5+6)=BSTR2(3,2)
      STATEV(6+6)=BSTR2(1,3)

      
      STATEV(13) = p
      STATEV(14) = R
      
      END IF
      
      RETURN
      END
      
      
*********************************************** 
**           UTILITY    SUBROUTINES          **            
***********************************************
**    MULTIPLY 6X6 MATRIX WITH 6X1 VECTOR    **
*USER SUBROUTINE
      SUBROUTINE KMLT1(DM1,DM2,DM,NTENS)
      
      INCLUDE 'ABA_PARAM.INC'
      
      PARAMETER (M=6)
      
      DIMENSION DM1(M,M),DM2(M),DM(M)
      
      DO 10 I=1,NTENS
          X=0.0
          DO 20 K=1,NTENS
              Y=DM1(I,K)*DM2(K)
              X=X+Y
20        CONTINUE
          DM(I)=X
10     CONTINUE
      
      RETURN
      END
          
***********************************************
**              EFFECTIVE STRESS             **
**      (CONTRACTED MATRIX CALCULATION)      **
***********************************************
*USER SUBROUTINE
      SUBROUTINE KEFFP(EFF1,VAL1)
      
      INCLUDE 'ABA_PARAM.INC'
      
      PARAMETER (M=3,N=3)
      
      DIMENSION EFF1(M,N)
      
      X=0.0
      DO 10 I=1,M
          DO 10 J=1,N
              X=X+EFF1(I,J)*EFF1(I,J)
10    CONTINUE
      IF(X .LE. 0.0) GO TO 20
          VAL1=DSQRT((3.0/2.0)*X)
20    RETURN
      END
      
***********************************************
**      DOT PRODUCT OF TWO VECTORS           **
***********************************************
*USER SUBROUTINE
      SUBROUTINE DOTPROD(DM1,DM2,DM,NTENS)
      
      INCLUDE 'ABA_PARAM.INC'
      
C      PARAMETER (M=6)
      
      DIMENSION DM1(6),DM2(6)
      
      Y=0.0
      DO 20 K=1,NTENS
          X=DM1(K)*DM2(K)
          Y=X+Y
20    CONTINUE
      DM=Y
      RETURN
      END
      
**
      SUBROUTINE DYADICPROD(DM1,DM2,DM3,NTENS)
      
      INCLUDE 'ABA_PARAM.INC'
      
C      PARAMETER (M=6)
      
      DIMENSION DM1(6),DM2(6),DM3(6,6)
      
      DO I=1,6
          DO J=1,6
              DM3(I,J)=DM1(I)*DM2(J)
          END DO
      END DO
      
      RETURN
      END
      
***********************************************
**      DEVIATORIC STRESS CALCULATION        **
***********************************************
*USER SUBROUTINE
      SUBROUTINE KDEVIA(STRSS,XIDENTY,DEVITO)
      
      INCLUDE 'ABA_PARAM.INC'
      
      PARAMETER (M=3,N=3)
      
      DIMENSION STRSS(M,N),XIDENTY(M,N),DEVITO(M,N)
      
      X=0.0
      DO 10 I=1,M
          DO 10 J=1,N
              IF(I .EQ. J)THEN
                  X=X+STRSS(I,J)
              ELSE
              END IF
10    CONTINUE
      
      DO 20 I=1,M
          DO 20 J=1,N
              IF(I .EQ. J)THEN
                  DEVITO(I,J)=STRSS(I,J)-((1./3.)*X*XIDENTY(I,J))
              ELSE
                  DEVITO(I,J)=STRSS(I,J)
              END IF
20    CONTINUE
      RETURN
      END
C*****************************BOTTOM*****************************