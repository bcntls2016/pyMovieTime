      SUBROUTINE SIMPS(K,IC,H,Y,S,N)
C
C **************************************************************
C       Integracio per el metode de Simpson per funcions
C       uniformement espaides
C
C       H ... Longitud del interval
C       Y ... Vector amb els valors de la funcio
C       S ... Vector amb els valors de la funcio integrada
C       N ... Numero de punts
C       K ... Numero de de la funcio que agafem per fer l'integral
C             sobre un interval
C       IC .. Parametre de control
C               -1 calcula la funcio integral
C                0 calcula integral total
C                1 torna els coeficients que fa servir per calcular
C                  les integrals
C                2 torna els pesos que fa servir per calcular
C                  les integrals
C
C       K     Te que esser parell i amb les condicions seguents
C                      2.ge.K.and.26.le.K
C                              i
C                           N.ge.K
C
C **************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SINTA/A(1547)
      DIMENSION Y(*),S(*)
C            MF  Utima posicio dels coeficients a fer servir
      MF=(2*K*K+3*K-2)*K/24
C            M   Nombre total de coeficients
      M=K*(K-1)/2
C            MI  Primera posicio dels coeficients a fer servir
      MI=MF-M+1
      IF(IC.EQ.1)THEN
        KKM1=K*(K-1)
        DO 5 I=MI,MF
          J=I-MI+1
          JJ=KKM1-J+1
          S(J)=A(I)
          S(JJ)=A(I)
5       CONTINUE
        RETURN
      ENDIF
      IF(MOD(K,2).NE.0.OR.K.GT.26.OR.K.GT.N)THEN
        WRITE(*,*)' Utilitzacio dolenta de la subrutina SIMPS...',K,N
        RETURN
      ENDIF
      KM1=K-1
      KO2=K/2
      KO2P1=KO2+1
      KO2M1=KO2-1
      KO2M2=KO2-2
      IF(IC.EQ.2)THEN
        DO I=1,N
          S(I)=0.0D0
        ENDDO
        DO I=2,KO2
          L0=MI+(I-2)*K
          DO J=1,K
            L=L0+J-1
            S(J)=S(J)+A(L)
          ENDDO
        ENDDO
        L0=MF-KO2M1
        DO  I=KO2P1,N-KO2M1
          J0=I-KO2
          DO J=J0,I-1
            J1=J-J0
            JJ=I+KO2M1-J1
            L=L0+J1
            S(J)=S(J)+A(L)
            S(JJ)=S(JJ)+A(L)
          ENDDO
        ENDDO
        I0=N-KO2M2
        J0=N-K+1
        DO I=I0,N
          L0=MI+(N-I+1)*K-1
          DO  J=J0,N
            L=L0-(J-J0)
            S(J)=S(J)+A(L)
          ENDDO
        ENDDO
        RETURN
      ENDIF
      S(1)=0.0D0
      AUX=0.0D0
      DO 10 I=2,KO2
        L0=MI+(I-2)*K
        DO 20 J=1,K
          L=L0+J-1
          AUX=AUX+A(L)*Y(J)
20      CONTINUE
        IF(IC.EQ.-1)THEN
          S(I)=S(I-1)+AUX*H
          AUX=0.0D0
        ENDIF
10    CONTINUE
      L0=MF-KO2M1
      DO 30 I=KO2P1,N-KO2M1
        J0=I-KO2
        DO 40 J=J0,I-1
          J1=J-J0
          JJ=I+KO2M1-J1
          L=L0+J1
          AUX=AUX+A(L)*(Y(J)+Y(JJ))
40      CONTINUE
        IF(IC.EQ.-1)THEN
          S(I)=S(I-1)+AUX*H
          AUX=0.0D0
        ENDIF
30    CONTINUE
      I0=N-KO2M2
      J0=N-K+1
      DO 50 I=I0,N
        L0=MI+(N-I+1)*K-1
        DO 60 J=J0,N
          L=L0-(J-J0)
          AUX=AUX+A(L)*Y(J)
60      CONTINUE
        IF(IC.EQ.-1)THEN
          S(I)=S(I-1)+AUX*H
          AUX=0.0D0
        ENDIF
50    CONTINUE
      IF(IC.EQ.0)THEN
        S(1)=AUX*H
      ENDIF
      RETURN
      END
