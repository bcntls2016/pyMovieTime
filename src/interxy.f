      SUBROUTINE INTERxy(Nx,Ny,Xms,Yms,Phi,Nx0,Ny0,Hx0,Hy0,Phi0,K,KM1)
C
C     Subrutina que serveix per interpolar les funcions d'ona i el
C     Coulomb, desde un grid (Nx0-Hx0,Ny0-Hy0) cap un grid
C     (Nx-Hx,Ny-Hy)
C
C     Nx0  : Punts en la direccio x on coneixem les funcions
C     Ny0  : Al mateix en la direccio y
C     Hx0   : Pas en la direccio x
C     Hy0   :  "   "  "     "    y
C
C     K  : Nombre de punts que fem servir per interpolar la funcio en
C          cadescuna de les 2 direccions
C     KM1: Nombre maxim de derivades que fem servir en el
C          desenvolopament de Taylor
C
C      Use Grid, only : Xms,Yms, Nx, Ny
C      Use Onas, only : Phi
C      Use Nivells, only : Nphi
C      Use Derint
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(KM=25,KM2=KM*KM,KM21=KM2+1,KT=(KM+1)*KM2)
      Dimension Xms(Nx), Yms(Ny)
      Complex (Kind=8):: Phi(Nx,Ny),Phi0(Nx0,Ny0),Csto
      DIMENSION CDC(0:KM,KM2),W(KM21),C1OFAC(0:KM),A(KM,KM)
      DIMENSION SX(0:KM),SY(0:KM)
      LOGICAl*1 INTRP
      DATA CDC/KT*0.0D0/
C
C     Control del posibles errors greus abans de comensar
C
      IF(K.GT.Nx0.OR.K.GT.Ny0.OR.KM1.GT.(K-1))THEN
        WRITE(6,*)' From Interxy: Error en els parametres...'
     &  ,K,KM1,Nx0,Ny0
        RETURN
      ENDIF
      Xmax0=Nx0/2*Hx0
      Xmin0=-Xmax0
      Ymax0=Ny0/2*Hy0
      Ymin0=-Ymax0
      XI=Xmin0
      YI=Ymin0
      XF=Xmax0-hx0
      YF=Ymax0-hy0
C
C      Write(6,'(" Interxy: Xi,Yi.....",1p,2E15.6)')Xi,Yi
C      Write(6,'(" Interxy: Xf,Yf.....",1p,2E15.6)')Xf,Yf
C      Call Flush(6)
C      Write(6,'(" Interxy: Xms(1),Yms(1).......",1p,2E15.6)')
C     &Xms(1),Yms(1)
C      Write(6,'(" Interxy: Xms(Nx),Yms(Nx).....",1p,2E15.6)')
C     &Xms(Nx),Yms(Ny)
C      Call Flush(6)
C
      KO2=K/2
      KO2P1=KO2+1
      KK=K*K
      KKP1=KK+1
      NXMKO2=Nx0-KO2
      NYMKO2=Ny0-KO2
      C1OFAC(0)=1.0D0
      C1OFAC(1)=1.0D0
C
C     Aqui defineixo els "coeficients" de la derivada d'ordre 0
C     o sia la funcio (ha d'esser la matriu unitat)
C
      DO I=1,K
        II=I+(I-1)*K
        CDC(0,II)=1.0D0
      ENDDO
      DO IK=1,KM1
        C1OFAC(IK+1)=C1OFAC(IK)/DFLOAT(IK+1)
        CALL NEWDER(IK,K,100,H,A,W,1)
        STO=1.0D0/W(KKP1)*C1OFAC(IK)
        DO I=1,KK
          CDC(IK,I)=W(I)*STO
        ENDDO
      ENDDO
      DO Ix=1,Nx
        DO Iy=1,Ny
C
C     Aqui inicilitzem la matriu d'interpolacio
C
          DO IXx=1,K
            DO IYy=1,K
              A(IXx,IYy)=0.0D0
            ENDDO
          ENDDO
          X=Xms(Ix)
          Y=Yms(Iy)
C
C     Aqui mirem quins son els index del punt mes proper a (X,Y,Z)
C     i tambe calculem les coordenades d'aquest punt
C
          IX0=(X-XI)/Hx0+1.5D0
          IY0=(Y-YI)/Hy0+1.5D0
          INTRP=.TRUE.
C
C     Aqui comprobem que no ens sortim de mare
C
          IF(IX0.LT.1)THEN
            IX0=1
            INTRP=.FALSE.
          ENDIF
          IF(IX0.GT.Nx0)THEN
            IX0=Nx0
            INTRP=.FALSE.
          ENDIF
          IF(IY0.LT.1)THEN
            IY0=1
            INTRP=.FALSE.
          ENDIF
          IF(IY0.GT.Ny0)THEN
            IY0=Ny0
            INTRP=.FALSE.
          ENDIF
          X0=(IX0-1)*Hx0+XI
          Y0=(IY0-1)*Hy0+YI
          IF(.NOT.INTRP)THEN
            GO TO 100
          ENDIF
C
C     Index auxiliars per poder calcular les derivades de F en el punt
C     (X0,Y0)
C
          IF(IX0.LE.KO2P1)THEN
            IXX0=(IX0-1)*K
            JXX0=IX0
          ELSEIF(IX0.LE.NXMKO2)THEN
            IXX0=KO2*K
            JXX0=KO2P1
          ELSE
            IXX0=K*(IX0+K-Nx0-1)
            JXX0=K+IX0-Nx0
          ENDIF
          IXX1=IX0-JXX0
          IF(IY0.LE.KO2P1)THEN
            IYY0=(IY0-1)*K
            JYY0=IY0
          ELSEIF(IY0.LE.NYMKO2)THEN
            IYY0=KO2*K
           JYY0=KO2P1
          ELSE
            IYY0=K*(IY0+K-Ny0-1)
            JYY0=K+IY0-Ny0
          ENDIF
          IYY1=IY0-JYY0
          SSX=(X-X0)/Hx0
          SSY=(Y-Y0)/Hy0
          SX(0)=1.0D0
          SY(0)=1.0D0
          DO I=1,KM1
            SX(I)=SX(I-1)*SSX
            SY(I)=SY(I-1)*SSY
          ENDDO
          DO JX=1,K
            JJX=IXX0+JX
            DO JY=1,K
              JJY=IYY0+JY
              DO L=1,KM1
                DO M=0,L
                  LMM=L-M
                  AUXX=SX(LMM)*CDC(LMM,JJX)
                  AUXY=SY(M)*CDC(M,JJY)
                  A(JX,JY)=A(JX,JY)
     &                    +AUXX*AUXY
                ENDDO
              ENDDO
            ENDDO
          ENDDO
100       CONTINUE
            CSTO=PHI0(IX0,IY0)
            IF(INTRP)THEN
              DO JX=1,K
                IIX=IXX1+JX
                DO JY=1,K
                  IIY=IYY1+JY
                  CSTO=CSTO+A(JX,JY)*PHI0(IIX,IIY)
                ENDDO
              ENDDO
            ENDIF
            PHI(Ix,Iy)=CSTO
          ENDDO
        ENDDO
      RETURN
      END
