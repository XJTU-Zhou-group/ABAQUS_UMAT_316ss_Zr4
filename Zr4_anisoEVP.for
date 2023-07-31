C********************************TOP***************************
C     UMAT FOR ABAQUS/STANDARD
C     AN ANISOTROPIC ELASTO-VISCOPLASTIC MODEL FOR ZIRCALOY-4
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
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),
     5 PSTRESS(NTENS),SDEV(6),DSTRESS(NTENS)
      
      PARAMETER (M=6,N=6,ID=3,ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 SIX=6.D0, NINE=9.D0, TOLER=1.D-5)
      
      DIMENSION BSTRESS(6),BSTRESS1(6),BSTRESS2(6),STRESSOLD(6),
     1 VSTRESS(6),EBSTRESS(6),XMM(M,N),XMN(M,N),XMQ(M,N),XMR(M,N),
     2 PFDIR(6),BSTRESST(6),BSTRESS1T(6),BSTRESS2T(6),XMNMPXNV(6),
     3 XMNMPXMR(6,6),EQBSTRESSV(6),EQDBSTRESSV(6),PFDIR1(6),
     4 DBSTRESS2(6),DBSTRESS1(6),DBSTRESS(6),
     5 DPSTRAN(6),DESTRAN(6),XHILSTRT(6),EQBSTRESST(6),
     6 BSTRMPXMR(6),XNVMPXMN(6),XNVMPXMQ(6),
     7 PFDIRMPXMN(6),PFDIRMPXMQ(6),XMNMPPFDIR(6),XMQMPBSTR2(6),
     8 XMQMPBSTR1(6),XMQMPBSTR(6),PFDIRMPXMNMPXMR(6),XMNMPXMRMPBSTR(6)
     9
      
      
      REAL*8 E,XNUE,XLN,XSN,EPSRZ,XQ,XQ1,XQ2,XRM,XRM1,ALPHAZ,XMZ,XM1
      REAL*8 YIELDZ,YSATISO,YSATSITA,B,BSITA,XM11,XM22,XM33,XM13,XM23
      REAL*8 XM12,XM44,XN11,XN22,XN33,XN13,XN23,XN12,XN44,VMSTR,xdp
      REAL*8 XNSITA,SITA,XFSITA,XMM,XMN,XMQ,XMR
      

      DO K1 = 1,6
          DO K2 = 1,6
              DDSDDE(K2,K1) = 0.D0
          END DO
      END DO
      XMM = XMM * 0.D0
      XMN = XMN * 0.D0
      XMQ = XMQ * 0.D0
      XMR = XMR * 0.D0
      
      
      E = PROPS(1)
      XNUE = PROPS(2)
      XLN = PROPS(3)
      XSN = PROPS(4)
      EPSRZ = PROPS(5)
      XQ = PROPS(6)
      XQ1 = PROPS(7)
      XQ2 = PROPS(8)
      XRM = PROPS(9)
      XRM1 = PROPS(10)
      ALPHAZ = PROPS(11)
      XMZ = PROPS(12)
      XM1 = PROPS(13)
      YIELDZ = PROPS(14)
      YSATISO = PROPS(15)
      YSATSITA = PROPS(16)
      B = PROPS(17)
      BSITA = PROPS(18)
      XM11 = PROPS(19)
      XM22 = PROPS(20)
      XM33 = PROPS(21)
      XM13 = PROPS(22)
      XM23 = PROPS(23)
      XM12 = PROPS(24)
      XM44 = PROPS(25)
      XN11 = PROPS(26)
      XN22 = PROPS(27)
      XN33 = PROPS(28)
      XN13 = PROPS(29)
      XN23 = PROPS(30)
      XN12 = PROPS(31)
      XN44 = PROPS(32)
      Q11 = PROPS(33)
      Q22 = PROPS(34)
      Q33 = PROPS(35)
      Q44 = PROPS(36)
      Q13 = PROPS(37)
      Q23 = PROPS(38)
      Q12 = PROPS(39)
      R11 = PROPS(40)
      R22 = PROPS(41)
      R33 = PROPS(42)
      R13 = PROPS(43)
      R23 = PROPS(44)
      R12 = PROPS(45)
      R44 = PROPS(46)
      
      EBULK3 = E/(ONE-TWO*XNUE)
      XK = EBULK3/THREE
      EG2 = E/(ONE+XNUE)
      EG = EG2/TWO
      ELAM = (EBULK3-EG2)/THREE
       
      
      !XM11 = TWO/THREE
      !XM22 = TWO/THREE
      !XM33 = TWO/THREE
      !XM13 = -ONE/THREE
      !XM23 = -ONE/THREE
      !XM12 = -ONE/THREE
      !XM44 = TWO

      
      !XN11 = TWO/THREE
      !XN22 = TWO/THREE
      !XN33 = TWO/THREE
      !XN13 = -ONE/THREE
      !XN23 = -ONE/THREE
      !XN12 = -ONE/THREE
      !XN44 = ONE/TWO      
      
      
      !Q11 = TWO/THREE
      !Q22 = TWO/THREE
      !Q33 = TWO/THREE
      !Q13 = -ONE/THREE
      !Q23 = -ONE/THREE
      !Q12 = -ONE/THREE
      !Q44 = ONE
           
      
      !R11 = TWO/THREE
      !R22 = TWO/THREE
      !R33 = TWO/THREE
      !R13 = -ONE/THREE
      !R23 = -ONE/THREE
      !R12 = -ONE/THREE
      !R44 = TWO

      DO K1 = 1,3
          DO K2 = 1,3
              DDSDDE(K2,K1) = ELAM
          END DO
          DDSDDE(K1,K1) = EG2+ELAM
      END DO
      
      DDSDDE(4,4)=EG
      DDSDDE(5,5)=EG
      DDSDDE(6,6)=EG
      
      DO K=1,6
          BSTRESS(K)=STATEV(K)
      END DO
      
      DO K=1,6
          BSTRESS1(K)=STATEV(K+6)
      END DO
      
      DO K=1,6
          BSTRESS2(K)=STATEV(K+12)
      END DO
      
      YISO = STATEV(19)
      YSITA = STATEV(20)
      p = STATEV(21)
      
      STATEV(61)=STRESS(1)
      STATEV(62)=STRESS(2)
      STATEV(63)=STRESS(3)
      STATEV(64)=STRESS(4)
      STATEV(65)=STRESS(5)
      STATEV(66)=STRESS(6)
      
      DO K=1,6
          STRESSOLD(K)=STRESS(K)
      END DO
      
      DO K = 1,6
          DSTRESS(K)=0.0D0
      END DO
      
      DO K1 = 1,NTENS
         DO K2 = 1,NTENS
              DSTRESS(K1)=DSTRESS(K1)+DDSDDE(K1,K2)*DSTRAN(K2)
         END DO
      END DO
      
      DO K=1,NTENS
          STRESS(K)=STRESS(K)+DSTRESS(K)
      END DO 
      
      DO K = 1,3
          VSTRESS(K)=STRESS(K)-((ONE/THREE)*
     1    (STRESS(1)+STRESS(2)+STRESS(3)))
      END DO
      VSTRESS(4)=STRESS(4)
      VSTRESS(5)=STRESS(5)
      VSTRESS(6)=STRESS(6)
      
      DO K=1,NTENS
          EBSTRESS(K)=VSTRESS(K)-BSTRESS(K)
      END DO
           
      
      XHILSTR1=XM11*EBSTRESS(1)*EBSTRESS(1)+XM22*EBSTRESS(2)*EBSTRESS(2)
     1+XM33*EBSTRESS(3)*EBSTRESS(3)+TWO*XM12*EBSTRESS(1)*EBSTRESS(2)+
     2TWO*XM23*EBSTRESS(2)*EBSTRESS(3)+TWO*XM13*EBSTRESS(1)*EBSTRESS(3)+
     3XM44*EBSTRESS(4)*EBSTRESS(4)+XM44*EBSTRESS(5)*EBSTRESS(5)+
     4XM44*EBSTRESS(6)*EBSTRESS(6)
      
      
      XHILSTR2=(THREE/TWO)*XHILSTR1
      XHILSTR3 = DSQRT(XHILSTR2)
      

      STATEV(23)=XMM(1,1)
      STATEV(24)=EBSTRESS(1)
      STATEV(25)=EBSTRESS(2)
      STATEV(26)=EBSTRESS(3)
      STATEV(27)=XHILSTRT(1)
      STATEV(28)=XHILSTRT(2)
      STATEV(29)=XHILSTRT(3)          
      STATEV(30)=XHILSTR1
      STATEV(31)=XHILSTR2     
      STATEV(32)=XHILSTR3         
      STATEV(33)=VSTRESS(1)      
      STATEV(34)=VSTRESS(2)      
      STATEV(35)=VSTRESS(3)      
      STATEV(36)=EBSTRESS(4)
      STATEV(37)=EBSTRESS(5)
      STATEV(38)=EBSTRESS(6)
      STATEV(39)=XHILSTRT(4)
      STATEV(40)=XHILSTRT(5)
      STATEV(41)=XHILSTRT(6)


      DO K = 1,NTENS
          PFDIR1(K)=0.D0
      END DO

      PFDIR1(1)=XM11*EBSTRESS(1)+XM12*EBSTRESS(2)+XM13*EBSTRESS(3)
      PFDIR1(2)=XM12*EBSTRESS(1)+XM22*EBSTRESS(2)+XM23*EBSTRESS(3)
      PFDIR1(3)=XM13*EBSTRESS(1)+XM23*EBSTRESS(2)+XM33*EBSTRESS(3)      
      PFDIR1(4)=XM44*EBSTRESS(4)
      PFDIR1(5)=XM44*EBSTRESS(5)
      PFDIR1(6)=XM44*EBSTRESS(6)
      
      DO K = 1,NTENS
          PFDIR(K)=0.D0
      END DO
      
      DO K = 1,NTENS
          PFDIR(K)=(THREE/TWO)*PFDIR1(K)/XHILSTR3
      END DO       
      
      YF=XHILSTR3

      IF(YF.GT.0.)THEN
      
      YISO0=YISO
      YSITA0=YSITA
      
      DO K=1,6
          BSTRESST(K)=BSTRESS(K)
          BSTRESS1T(K)=BSTRESS1(K)
          BSTRESS2T(K)=BSTRESS2(K)
      END DO
      xdp=1.D-10
      DO KNEWT =1,10
          XFLOW=XHILSTR3-(THREE*EG*xdp)
          XPHI=EPSRZ*((sinh(XFLOW/XLN))**XSN)
          XPHISTR=EPSRZ*XSN*(ONE/XLN)*((sinh(XFLOW/XLN))**(XSN-ONE))*
     1    cosh(XFLOW/XLN)
!          XPHIY=-EPSRZ*XSN*(ONE/XLN)*((sinh(XFLOW/XLN))**(XSN-ONE))*
!     1    cosh(XFLOW/XLN)
          XPHIDP=-EPSRZ*XSN*(ONE/XLN)*((sinh(XFLOW/XLN))**(XSN-ONE))*
     1    (cosh(XFLOW/XLN))*THREE*EG
          XRES=XPHI-(xdp/DTIME)
          
          STATEV(67)=XRES
          STATEV(82)=XFLOW
          STATEV(83)=XPHI          
          STATEV(84)=XPHISTR               
          STATEV(85)=XPHIDP          
          
          PFDIRMPXMNMPPFDIR=0.0D0
          PFDIRMPXMNMPPFDIR=XN11*PFDIR(1)*PFDIR(1)+XN22*PFDIR(2)*
     1    PFDIR(2)+XN33*PFDIR(3)*PFDIR(3)+TWO*XN12*PFDIR(1)*PFDIR(2)+
     2    TWO*XN23*PFDIR(2)*PFDIR(3)+TWO*XN13*PFDIR(1)*PFDIR(3)+
     3    XN44*PFDIR(4)*PFDIR(4)+XN44*PFDIR(5)*PFDIR(5)+
     4    XN44*PFDIR(6)*PFDIR(6)          
          
          STATEV(68)=PFDIRMPXMNMPPFDIR
          
          
          PFDIRMPXMQMPBSTR1=0.0D0
          PFDIRMPXMQMPBSTR1=Q11*PFDIR(1)*(BSTRESS(1)-BSTRESS1(1))+
     1    Q22*PFDIR(2)*(BSTRESS(2)-BSTRESS1(2))+Q33*PFDIR(3)*
     1    (BSTRESS(3)-BSTRESS1(3))+Q12*PFDIR(1)*(BSTRESS(2)-BSTRESS1(2))
     1    +Q12*PFDIR(2)*(BSTRESS(1)-BSTRESS1(1))+Q23*PFDIR(2)*
     1    (BSTRESS(3)-BSTRESS1(3))+Q23*PFDIR(3)*(BSTRESS(2)-BSTRESS1(2))
     1    +Q13*PFDIR(1)*(BSTRESS(3)-BSTRESS1(3))+Q13*PFDIR(3)*
     1    (BSTRESS(1)-BSTRESS1(1))+Q44*PFDIR(4)*(BSTRESS(4)-BSTRESS1(4))
     1    +Q44*PFDIR(5)*(BSTRESS(5)-BSTRESS1(5))+Q44*PFDIR(6)
     1    *(BSTRESS(6)-BSTRESS1(6))
        
          STATEV(69)=PFDIRMPXMQMPBSTR1
            
          
          XNR11=0.D0
          XNR12=0.D0
          XNR13=0.D0
          XNR21=0.D0
          XNR22=0.D0
          XNR23=0.D0
          XNR31=0.D0
          XNR32=0.D0
          XNR33=0.D0
          XNR44=0.D0
          XNR55=0.D0
          XNR66=0.D0
          
          XNR11=XN11*R11+XN12*R12+XN13*R13
          XNR12=XN11*R12+XN12*R22+XN13*R23
          XNR13=XN11*R13+XN12*R23+XN13*R33 
          XNR21=XN12*R11+XN22*R12+XN23*R13
          XNR22=XN12*R12+XN22*R22+XN23*R23
          XNR23=XN12*R13+XN22*R23+XN23*R33
          XNR31=XN13*R11+XN23*R12+XN33*R13
          XNR32=XN13*R12+XN23*R22+XN33*R23
          XNR33=XN13*R13+XN23*R23+XN33*R33
          XNR44=XN44*R44
          XNR55=XN44*R44
          XNR66=XN44*R44          
          
          PFDIRMPXMNMPXMRMPBSTR=XNR11*PFDIR(1)*BSTRESS(1)+
     1    XNR22*PFDIR(2)*BSTRESS(2)+XNR33*PFDIR(3)*
     1    BSTRESS(3)+XNR12*PFDIR(1)*BSTRESS(2)
     1    +XNR21*PFDIR(2)*BSTRESS(1)+XNR23*PFDIR(2)*
     1    BSTRESS(3)+XNR32*PFDIR(3)*BSTRESS(2)
     1    +XNR13*PFDIR(1)*BSTRESS(3)+XNR31*PFDIR(3)*
     1    BSTRESS(1)+XNR44*PFDIR(4)*BSTRESS(4)
     1    +XNR55*PFDIR(5)*BSTRESS(5)+XNR66*PFDIR(6)
     1    *BSTRESS(6)
          
          STATEV(70)=PFDIRMPXMNMPXMRMPBSTR          
         
          EQBSTRESS = 0.0D0          
          EQBSTRESS=R11*BSTRESS(1)*BSTRESS(1)+R22*BSTRESS(2)*BSTRESS(2)
     1    +R33*BSTRESS(3)*BSTRESS(3)+TWO*R12*BSTRESS(1)*BSTRESS(2)+
     2    TWO*R23*BSTRESS(2)*BSTRESS(3)+TWO*R13*BSTRESS(1)*BSTRESS(3)+
     3    R44*BSTRESS(4)*BSTRESS(4)+R44*BSTRESS(5)*BSTRESS(5)+
     4    R44*BSTRESS(6)*BSTRESS(6)          
          EQBSTRESS = DSQRT((THREE/TWO)*EQBSTRESS)
          
          STATEV(71)=EQBSTRESS           
        

          IF(EQBSTRESS.LT.1.D-5)THEN
          EQBSTRESS=1.D-5
          END IF          
          
          STATEV(72)=EQBSTRESS
          
          
          DEQPL=(XRES+XPHISTR*((XRM*((EQBSTRESS/ALPHAZ)**XMZ))+(XRM1*
     1    ((EQBSTRESS/ALPHAZ)**XM1)))*DTIME*(PFDIRMPXMNMPXMRMPBSTR
     2    /EQBSTRESS))/
     3    ((ONE/DTIME)-XPHIDP+((TWO/THREE)*(YIELDZ+YISO+YSITA)*
     4    PFDIRMPXMNMPPFDIR-PFDIRMPXMQMPBSTR1)*XPHISTR*XQ)
          
          
          STATEV(73)=DEQPL
          
          xdp=xdp+DEQPL
          
          DO K = 1,6
              XMNMPPFDIR(K)=0.0D0
          END DO          
          XMNMPPFDIR(1)=XN11*PFDIR(1)+XN12*PFDIR(2)+XN13*PFDIR(3)
          XMNMPPFDIR(2)=XN12*PFDIR(1)+XN22*PFDIR(2)+XN23*PFDIR(3)
          XMNMPPFDIR(3)=XN13*PFDIR(1)+XN23*PFDIR(2)+XN33*PFDIR(3)      
          XMNMPPFDIR(4)=XN44*PFDIR(4)
          XMNMPPFDIR(5)=XN44*PFDIR(5)
          XMNMPPFDIR(6)=XN44*PFDIR(6)          
          
          DO K = 1,6
              XMQMPBSTR2(K)=0.0D0
          END DO          
          XMQMPBSTR2(1)=Q11*BSTRESS2T(1)+Q12*BSTRESS2T(2)+Q13*
     1    BSTRESS2T(3)
          XMQMPBSTR2(2)=Q12*BSTRESS2T(1)+Q22*BSTRESS2T(2)+Q23*
     1    BSTRESS2T(3)
          XMQMPBSTR2(3)=Q13*BSTRESS2T(1)+Q23*BSTRESS2T(2)+Q33*
     1    BSTRESS2T(3)      
          XMQMPBSTR2(4)=Q44*BSTRESS2T(4)
          XMQMPBSTR2(5)=Q44*BSTRESS2T(5)
          XMQMPBSTR2(6)=Q44*BSTRESS2T(6)          
          
          
          DO K = 1,6
              XMQMPBSTR1(K)=0.0D0
          END DO          
          XMQMPBSTR1(1)=Q11*(BSTRESS1T(1)-BSTRESS2T(1))+Q12*
     1    (BSTRESS1T(2)-BSTRESS2T(2))+Q13*(BSTRESS1T(3)-BSTRESS2T(3))
          XMQMPBSTR1(2)=Q12*(BSTRESS1T(1)-BSTRESS2T(1))+Q22*
     1    (BSTRESS1T(2)-BSTRESS2T(2))+Q23*(BSTRESS1T(3)-BSTRESS2T(3))
          XMQMPBSTR1(3)=Q13*(BSTRESS1T(1)-BSTRESS2T(1))+Q23*
     1    (BSTRESS1T(2)-BSTRESS2T(2))+Q33*(BSTRESS1T(3)-BSTRESS2T(3))      
          XMQMPBSTR1(4)=Q44*(BSTRESS1T(4)-BSTRESS2T(4))
          XMQMPBSTR1(5)=Q44*(BSTRESS1T(5)-BSTRESS2T(5))
          XMQMPBSTR1(6)=Q44*(BSTRESS1T(6)-BSTRESS2T(6))            
          
        
          DO K = 1,6
              XMQMPBSTR(K)=0.0D0
          END DO          
          XMQMPBSTR(1)=Q11*(BSTRESST(1)-BSTRESS1T(1))+Q12*
     1    (BSTRESST(2)-BSTRESS1T(2))+Q13*(BSTRESST(3)-BSTRESS1T(3))
          XMQMPBSTR(2)=Q12*(BSTRESST(1)-BSTRESS1T(1))+Q22*
     1    (BSTRESST(2)-BSTRESS1T(2))+Q23*(BSTRESST(3)-BSTRESS1T(3))
          XMQMPBSTR(3)=Q13*(BSTRESST(1)-BSTRESS1T(1))+Q23*
     1    (BSTRESST(2)-BSTRESS1T(2))+Q33*(BSTRESST(3)-BSTRESS1T(3))      
          XMQMPBSTR(4)=Q44*(BSTRESST(4)-BSTRESS1T(4))
          XMQMPBSTR(5)=Q44*(BSTRESST(5)-BSTRESS1T(5))
          XMQMPBSTR(6)=Q44*(BSTRESST(6)-BSTRESS1T(6))           
          
                  
          
          DO K = 1,6
              XMNMPXMRMPBSTR(K)=0.0D0
          END DO
          XMNMPXMRMPBSTR(1)=XNR11*BSTRESST(1)+XNR12*BSTRESST(2)
     1    +XNR13*BSTRESST(3)
          XMNMPXMRMPBSTR(2)=XNR21*BSTRESST(1)+XNR22*BSTRESST(2)
     1    +XNR23*BSTRESST(3)
          XMNMPXMRMPBSTR(3)=XNR31*BSTRESST(1)+XNR32*BSTRESST(2)
     1    +XNR33*BSTRESST(3)      
          XMNMPXMRMPBSTR(4)=XNR44*BSTRESST(4)
          XMNMPXMRMPBSTR(5)=XNR55*BSTRESST(5)
          XMNMPXMRMPBSTR(6)=XNR66*BSTRESST(6)          
          
          
          EQBSTRESS = 0.0D0          
          EQBSTRESS=R11*BSTRESST(1)*BSTRESST(1)+R22*BSTRESST(2)*
     1    BSTRESST(2)+R33*BSTRESST(3)*BSTRESST(3)+TWO*R12*BSTRESST(1)*
     2    BSTRESST(2)+TWO*R23*BSTRESST(2)*BSTRESST(3)+TWO*R13*
     3    BSTRESST(1)*BSTRESST(3)+R44*BSTRESST(4)*BSTRESST(4)+R44*
     4    BSTRESST(5)*BSTRESST(5)+R44*BSTRESST(6)*BSTRESST(6)          
          EQBSTRESS = DSQRT((THREE/TWO)*EQBSTRESS)          
          
          STATEV(74)=EQBSTRESS

          
          IF(EQBSTRESS.LT.1.D-5)THEN
          EQBSTRESS=1.D-5
          END IF           
       
          STATEV(75)=EQBSTRESS
          
          DO K=1,6
              DBSTRESS2(K)=XQ2*((TWO/THREE)*(YIELDZ+YISO0+YSITA0)*
     1        XMNMPPFDIR(K)-XMQMPBSTR2(K))*xdp
          END DO
          
          DO K=1,6
              DBSTRESS1(K)=XQ1*((TWO/THREE)*(YIELDZ+YISO0+YSITA0)*
     1        XMNMPPFDIR(K)-XMQMPBSTR1(K))*xdp
          END DO
          
          DO K=1,6
              DBSTRESS(K)=XQ*((TWO/THREE)*(YIELDZ+YISO0+YSITA0)*
     1        XMNMPPFDIR(K)-XMQMPBSTR(K))*xdp-
     2        (((XRM*((EQBSTRESS/ALPHAZ)**XMZ))+(XRM1*
     3        ((EQBSTRESS/ALPHAZ)**XM1)))*DTIME*(XMNMPXMRMPBSTR(K)
     4        /EQBSTRESS))
          END DO
          
          DO K=1,6
              BSTRESS2(K)=BSTRESS2T(K)+DBSTRESS2(K)
              BSTRESS1(K)=BSTRESS1T(K)+DBSTRESS1(K)
              BSTRESS(K)=BSTRESST(K)+DBSTRESS(K)
          END DO
          
   
          EQDBSTRESS = 0.0D0          
          EQDBSTRESS=R11*DBSTRESS(1)*DBSTRESS(1)+R22*DBSTRESS(2)
     1    *DBSTRESS(2)+R33*DBSTRESS(3)*DBSTRESS(3)+TWO*R12*DBSTRESS(1)
     2    *DBSTRESS(2)+TWO*R23*DBSTRESS(2)*DBSTRESS(3)+TWO*R13*
     3    DBSTRESS(1)*DBSTRESS(3)+R44*DBSTRESS(4)*DBSTRESS(4)+
     4    R44*DBSTRESS(5)*DBSTRESS(5)+R44*DBSTRESS(6)*DBSTRESS(6)
          EQDBSTRESS = DSQRT((THREE/TWO)*EQDBSTRESS) 
          
          STATEV(76)=EQDBSTRESS
          
          IF(EQDBSTRESS.LT.1.D-15)THEN
          EQDBSTRESS=1.D-15
          END IF
          
          STATEV(77)=EQDBSTRESS
          
          
          BSTRMPXMRMPDBSTR=0.0D0
          BSTRMPXMRMPDBSTR=R11*BSTRESST(1)*DBSTRESS(1)+
     1    R22*BSTRESST(2)*DBSTRESS(2)+R33*BSTRESST(3)*
     1    DBSTRESS(3)+R12*BSTRESST(1)*DBSTRESS(2)
     1    +R12*BSTRESST(2)*DBSTRESS(1)+R23*BSTRESST(2)*
     1    DBSTRESS(3)+R23*BSTRESST(3)*DBSTRESS(2)
     1    +R13*BSTRESST(1)*DBSTRESS(3)+R13*BSTRESST(3)*
     1    DBSTRESS(1)+R44*BSTRESST(4)*DBSTRESS(4)
     1    +R44*BSTRESST(5)*DBSTRESS(5)+R44*BSTRESST(6)
     1    *DBSTRESS(6)          
          
          STATEV(78)=BSTRMPXMRMPDBSTR
               
          
          IF(EQDBSTRESS.LT.1.D-5)THEN
          EQDBSTRESS=1.D-5
          END IF
          
          IF(EQBSTRESS.LT.1.D-5)THEN
          EQBSTRESS=1.D-5
          END IF
          
          XNSITA=((THREE/TWO)*BSTRMPXMRMPDBSTR)/
     1    (EQDBSTRESS*EQBSTRESS)
          
          STATEV(79)=XNSITA          
          
          XFSITA=ONE-ABS(XNSITA)
          
          STATEV(80)=XFSITA        
                    
          YISO=YISO0+B*(YSATISO-YISO0)*xdp
          YSITA=YSITA0+BSITA*((YSATSITA*XFSITA)-YSITA0)*xdp
          
          IF(ABS(XRES).LT.1.E-12)GOTO 10
      END DO
10    CONTINUE          

      STATEV(81)=KNEWT

      DO K = 1,3
          DPSTRAN(K)=PFDIR(K)*xdp
      END DO
      DPSTRAN(4)=TWO*PFDIR(4)*xdp
      DPSTRAN(5)=TWO*PFDIR(5)*xdp
      DPSTRAN(6)=TWO*PFDIR(6)*xdp     
      
      DO K=1,6
          DESTRAN(K)=DSTRAN(K)-DPSTRAN(K)
      END DO
      
      DO K = 1,6
          DSTRESS(K)=0.0D0
      END DO
      
      DO K1 = 1,NTENS
         DO K2 = 1,NTENS
              DSTRESS(K1)=DSTRESS(K1)+DDSDDE(K1,K2)*DESTRAN(K2)
         END DO
      END DO
      
      STATEV(42)=PFDIR1(1)      
      STATEV(43)=PFDIR1(2)      
      STATEV(44)=PFDIR1(3)      
      STATEV(45)=PFDIR1(4)
      STATEV(46)=PFDIR1(5)
      STATEV(47)=PFDIR1(6)
      STATEV(48)=DPSTRAN(1)
      STATEV(49)=DPSTRAN(2)
      STATEV(50)=DPSTRAN(3)
      STATEV(51)=DPSTRAN(4)
      STATEV(52)=DPSTRAN(5)
      STATEV(53)=DPSTRAN(6)      
      STATEV(54)=PFDIR(1)      
      STATEV(55)=PFDIR(2)      
      STATEV(56)=PFDIR(3)      
      STATEV(57)=PFDIR(4)
      STATEV(58)=PFDIR(5)
      STATEV(59)=PFDIR(6)      
      STATEV(60)=xdp
      
      DO K=1,NTENS
          STRESS(K)=STRESSOLD(K)+DSTRESS(K)
      END DO
      
      p=p+xdp
      
      DO K=1,6
          STATEV(K)=BSTRESS(K)
      END DO
      
      DO K=1,6
          STATEV(K+6)=BSTRESS1(K)
      END DO
      
      DO K=1,6
          STATEV(K+12)=BSTRESS2(K)
      END DO
      
      STATEV(19) = YISO
      STATEV(20) = YSITA
      STATEV(21) = p

      END IF
      
      STATEV(22)=XFLOW
      RETURN
      END
C*****************************BOTTOM*****************************