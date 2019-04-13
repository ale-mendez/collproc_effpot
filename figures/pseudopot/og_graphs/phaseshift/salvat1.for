c**==  DEMORAD.FOR                                                       ADBP0030
C                                                                       ADBP0031
C                                                                       ADBP0032
C         ********************************************                  ADBP0033
C         **   DEMO FOR SUBROUTINE PACKAGE RADIAL   **                  ADBP0034
C         ********************************************                  ADBP0035
C                                                                       ADBP0036
C             EXPONENTIALLY SCREENED COULOMB FIELDS.                    ADBP0037
C                                                                       ADBP0038
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0039
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP0040
      PARAMETER (SL=137.036D0,PI=3.1415926535897932D0)                  ADBP0041
      DIMENSION R0(NDIM),RV0(NDIM)                                      ADBP0042
C  ****  COULOMB WAVE FUNCTION PARAMETERS.                              ADBP0043
      COMMON/OCOUL/WAVNUM,ETA,DELTA                                     ADBP0044
C  ****  OUTPUT RADIAL FUNCTIONS.                                       ADBP0045
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER  

	parameter (NPOT=500)
	DIMENSION A(NPOT),B(NPOT),C(NPOT),D(NPOT)
	dimension r0pot(NPOT),zPOT(NPOT)
	dimension Psq(ndim),rmed(ndim),ovrmed(ndim)

C                	                                                       ADBP0047
C  ****  READ FIELD PARAMETERS.                                         ADBP0048
C                                                                       ADBP0049
   10 CONTINUE                                                          ADBP0050

C                                                                       ADBP0055
C  ****  POTENTIAL GRID.                                                ADBP0056
C                                                                       ADBP0057
      RATIO=1.15D0                                                      ADBP0058
      RNN=100.0D0/DMAX1(ALPHA,1.0D0)                                     ADBP0059
      NV=500                                                            ADBP0060
      STEP=RNN/(NV-100.0D0)                                             ADBP0061
      CALL GRID(R0,RATIO,STEP,RNN,NV)
      rmax = r0(nv)
      
c.... read potential
	call aropm(r0pot,zPOT,A,B,C,D,npts,rmax)
	do i=1,NV
        rr = r0(i)
	  RV0(i) = POTE(rr,r0pot,zPOT,A,B,C,D,npts)
	enddo

      CALL VINT(R0,RV0,NV)                                              ADBP0077
C                                                                       ADBP0078
   30 CONTINUE                                                          ADBP0079
      WRITE(6,*) '  '                                                   ADBP0080
      WRITE(6,*) '  SELECT ONE OPTION ...'                              ADBP0081
      WRITE(6,*) '    1: SCHRODINGER EQUATION. BOUND STATE.'            ADBP0082
      WRITE(6,*) '    2: SCHRODINGER EQUATION. FREE STATE.'             ADBP0083
      WRITE(6,*) '    3: DIRAC EQUATION. BOUND STATE.'                  ADBP0084
      WRITE(6,*) '    4: DIRAC EQUATION. FREE STATE.'                   ADBP0085
      READ(5,*) IOPT                                                    ADBP0086
      WRITE(6,*) '  '                                                   ADBP0087
C                                                                       ADBP0088
C  ****  SCHRODINGER EQUATION. BOUND STATE.                             ADBP0089
C                                                                       ADBP0090
      IF(IOPT.EQ.1) THEN                                                ADBP0091
        WRITE(6,*) '  ENTER N, L AND EPS ...'                           ADBP0092
        READ(5,*) N,L,EPS                                               ADBP0093
        EPS=DMAX1(EPS,1.0D-15)                                          ADBP0094
        IF(N.LT.1.OR.L.GE.N) GO TO 30                                   ADBP0095
        RN=2.0D3                                                        ADBP0096
        NGP=800                                                         ADBP0097
        STEP=RN/(NGP-100.0D0)                                           ADBP0098
        CALL GRID(RAD,RATIO,STEP,RN,NGP)                                ADBP0099
        E=-Z**2/(2.0D0*N*N)                                             ADBP0100
        CALL SBOUND(E,EPS,EPS,N,L)                                      ADBP0101
        IF(IER.NE.0) GO TO 30                                           ADBP0102
        WRITE(6,1100) Z,ZS,ALPHA,N,L,EPS,E                              ADBP0103
 1100 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           ADBP0104
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    ADBP0105
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      ADBP0106
     3  /7X,'BOUND STATE: N=',I4,', L=',I4,'   (EPS=',E8.1,')'          ADBP0107
     4  /7X,'BINDING ENERGY =',E22.15)                                  ADBP0108


c..... mean values
	do i=1,npts
		ovrmed(i) = 0.0d0
		Psq(i)=P(i)*P(i)
		rmed(i) = Psq(i)*rad(i)
		if (rad(i).gt.1.d-14) ovrmed(i) = Psq(i)/rad(i)
	write(600,*) rad(i),' ' ,P(i)
	enddo
	ilast = npts-1
	rlast = Rad(ilast)
      CALL SPLINE(Rad,Psq,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP1028
      CALL INTEG(Rad,A,B,C,D,0.0D0,RLAST,SUMP,ILAST) 
	write(6,*) 'norm=',sump
      CALL SPLINE(Rad,rmed,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP1028
      CALL INTEG(Rad,A,B,C,D,0.0D0,RLAST,SUMP,ILAST) 
	write(6,*) '<r>=',sump
      CALL SPLINE(Rad,ovrmed,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP1028
      CALL INTEG(Rad,A,B,C,D,0.0D0,RLAST,SUMP,ILAST) 
	write(6,*) '<1/r>=',sump


C                                                                       ADBP0109
C  ****  SCHRODINGER EQUATION. FREE STATE.                              ADBP0110
C                                                                       ADBP0111
      ELSE IF(IOPT.EQ.2) THEN                                           ADBP0112
        WRITE(6,*) '  ENTER E, L AND EPS ...'                           ADBP0113
        READ(5,*) E,L,EPS                                               ADBP0114
        EPS=DMAX1(EPS,1.0D-15)                                          ADBP0115
        IF(E.LT.0.0D0.OR.L.LT.0) GO TO 30                               ADBP0116
        NGP=800                                                         ADBP0117
        WAVEL=2.0D0*PI/DSQRT(E+E)                                       ADBP0118
        STEP=0.05D0*WAVEL                                               ADBP0119
        RN=STEP*(NGP-100)                                               ADBP0120
        CALL GRID(RAD,RATIO,STEP,RN,NGP)                                ADBP0121
        CALL SFREE(E,EPS,PHASE,L)                                       ADBP0122
        IF(IER.NE.0) GO TO 30                                           ADBP0123
        WRITE(6,1200) Z,ZS,ALPHA,E,L,EPS,PHASE,DELTA,ETA                ADBP0124
 1200 FORMAT(/1X,1P,'****  SCHRODINGER EQ. ',                           ADBP0125
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    ADBP0126
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      ADBP0127
     3  /7X,'FREE STATE: E=',E13.6,', L=',I4,'   (EPS=',E8.1,')'        ADBP0128
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              ADBP0129
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')         ADBP0130
C                                                                       ADBP0131
C  ****  DIRAC EQUATION. BOUND STATE.                                   ADBP0132
C                                                                       ADBP0133
      ELSE IF(IOPT.EQ.3) THEN                                           ADBP0134
        WRITE(6,*) '  ENTER N, K AND EPS ...'                           ADBP0135
        READ(5,*) N,K,EPS                                               ADBP0136
        EPS=DMAX1(EPS,1.0D-15)                                          ADBP0137
        IF(N.LT.1.OR.K.EQ.0.OR.K.GE.N.OR.K.LT.-N) GO TO 30              ADBP0138
        RN=2.0D3                                                        ADBP0139
        NGP=800                                                         ADBP0140
        STEP=RN/(NGP-100.0D0)                                           ADBP0141
        CALL GRID(RAD,RATIO,STEP,RN,NGP)                                ADBP0142
        E=-Z**2/(2.0D0*N*N)                                             ADBP0143
        CALL DBOUND(E,EPS,EPS,N,K)                                      ADBP0144
        IF(IER.NE.0) GO TO 30                                           ADBP0145
        WRITE(6,1300) Z,ZS,ALPHA,N,K,EPS,E                              ADBP0146
 1300 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           ADBP0147
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    ADBP0148
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      ADBP0149
     3  /7X,'BOUND STATE: N=',I4,', K=',I4,'   (EPS=',E8.1,')'          ADBP0150
     4  /7X,'BINDING ENERGY =',E22.15)                                  ADBP0151
C                                                                       ADBP0152
C  ****  DIRAC EQUATION. BOUND STATE.                                   ADBP0153
C                                                                       ADBP0154
      ELSE IF(IOPT.EQ.4) THEN                                           ADBP0155
        WRITE(6,*) '  ENTER E, K AND EPS ...'                           ADBP0156
        READ(5,*) E,K,EPS                                               ADBP0157
        EPS=DMAX1(EPS,1.0D-15)                                          ADBP0158
        IF(E.LT.0.0D0.OR.K.EQ.0) GO TO 30                               ADBP0159
        IF(K.LT.0) THEN                                                 ADBP0160
          L=-K-1                                                        ADBP0161
        ELSE                                                            ADBP0162
          L=K                                                           ADBP0163
        ENDIF                                                           ADBP0164
        NGP=800                                                         ADBP0165
        WAVEL=2.0D0*PI/DSQRT(E*(2.0D0+E/SL**2))                         ADBP0166
        STEP=0.05D0*WAVEL                                               ADBP0167
        RN=STEP*(NGP-100)                                               ADBP0168
        CALL GRID(RAD,RATIO,STEP,RN,NGP)                                ADBP0169
        CALL DFREE(E,EPS,PHASE,K)                                       ADBP0170
        IF(IER.NE.0) GO TO 30                                           ADBP0171
        WRITE(6,1400) Z,ZS,ALPHA,E,K,EPS,PHASE,DELTA,ETA                ADBP0172
 1400 FORMAT(/1X,1P,'****  DIRAC EQUATION. ',                           ADBP0173
     1  'POTENTIAL FUNCTION: R*V(R)=Z+ZS*DEXP(-A*R)'                    ADBP0174
     2  /7X,'Z=',E13.6,', ZS =',E13.6,', A=',E13.6                      ADBP0175
     3  /7X,'FREE STATE: E=',E13.6,', K=',I4,'   (EPS=',E8.1,')'        ADBP0176
     4  /7X,'  INNER PHASE SHIFT=',E22.15,                              ADBP0177
     5  /7X,'COULOMB PHASE SHIFT=',E22.15,'   (ETA=',E13.6,')')         ADBP0178
      ELSE                                                              ADBP0179
        GO TO 10                                                        ADBP0180
      ENDIF                                                             ADBP0181
C                                                                       ADBP0182
C  ****  RADIAL WAVE FUNCTIONS WRITTEN ON FILE 'WAVES.DAT'.             ADBP0183
C                                                                       ADBP0184
      OPEN(7,FILE='wavesalvat.dat')                                          ADBP0185
      DO 40 I=1,NGP                                                     ADBP0186
      IF(DABS(P(I)).LT.1.0D-35) P(I)=1.0D-35                            ADBP0187
      IF(DABS(Q(I)).LT.1.0D-35) Q(I)=1.0D-35                            ADBP0188
      WRITE(7,'(1X,1P,4E13.5)') RAD(I),P(I),Q(I)                        ADBP0189
   40 CONTINUE                                                          ADBP0190
      CLOSE(UNIT=7)                                                     ADBP0191
      GO TO 30                                                          ADBP0192
      END                                                               ADBP0193
C  **************************************************************       ADBP0194
C                        SUBROUTINE GRID                                ADBP0195
C  **************************************************************       ADBP0196
      SUBROUTINE GRID(R,RATIO,STEP,RN,NP)                               ADBP0197
C                                                                       ADBP0198
C     THIS SUBROUTINE SETS UP A RADIAL GRID R(I) (I=1, ..., NP)         ADBP0199
C  SUCH THAT                                                            ADBP0200
C     1) R(1)=0, R(NP)=RN,                                              ADBP0201
C     2) A*R(I)+B*DLOG(R(I))-C=I  (I.GT.0), WITH                        ADBP0202
C        A=1.0/STEP AND B=1.0/DLOG(RATIO).                              ADBP0203
C                                                                       ADBP0204
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0205
      PARAMETER (NDIM=2400)                                              ADBP0206
      DIMENSION R(NDIM)                                                 ADBP0207
C                                                                       ADBP0208
      IF(NP.LT.50) THEN                                                 ADBP0209
        WRITE(6,1001) NP                                                ADBP0210
 1001 FORMAT(1X,'*** ERROR IN GRID: NP =',I3,                           ADBP0211
     1  ' MUST BE LARGER THAN 50.')                                     ADBP0212
        STOP                                                            ADBP0213
      ENDIF                                                             ADBP0214
      IF(NP.GT.NDIM) THEN                                               ADBP0215
        WRITE(6,1002) NP,NDIM                                           ADBP0216
 1002 FORMAT(1X,'*** ERROR IN GRID: NP =',I5,                           ADBP0217
     1  ' IS LARGER THAN NDIM =',I5,'.')                                ADBP0218
        STOP                                                            ADBP0219
      ENDIF                                                             ADBP0220
      IF(STEP.LT.1.0D-10.OR.(RATIO-1.0D0).LT.1.0D-3) THEN               ADBP0221
        WRITE(6,1003) STEP, RATIO                                       ADBP0222
 1003 FORMAT(1X,'*** ERROR IN GRID: STEP =',1P,E10.3,' OR RATIO',       ADBP0223
     1  ' =',E10.3,' ARE TOO SMALL.')                                   ADBP0224
        STOP                                                            ADBP0225
      ENDIF                                                             ADBP0226
C                                                                       ADBP0227
      A=1.0D0/STEP                                                      ADBP0228
      B=1.0D0/DLOG(RATIO)                                               ADBP0229
      C=NP-A*RN-B*DLOG(RN)                                              ADBP0230
C                                                                       ADBP0231
      R(1)=0.0D0                                                        ADBP0232
      RR=1.0D-35                                                        ADBP0233
      FR=A*RR+B*DLOG(RR)+C-1                                            ADBP0234
      IF(FR.GT.0.0D0) THEN                                              ADBP0235
        WRITE(6,1004)                                                   ADBP0236
 1004 FORMAT(1X,'*** ERROR IN GRID: R(2) IS TOO SMALL.')                ADBP0237
        STOP                                                            ADBP0238
      ENDIF                                                             ADBP0239
      DO 3 I=2,NP                                                       ADBP0240
      CI=C-I                                                            ADBP0241
      RL=RR                                                             ADBP0242
      RU=RL                                                             ADBP0243
    1 RU=RU+RU                                                          ADBP0244
      FU=A*RU+B*DLOG(RU)+CI                                             ADBP0245
      IF(FU.LT.0.0D0) GO TO 1                                           ADBP0246
    2 RR=0.5D0*(RU+RL)                                                  ADBP0247
      FR=A*RR+B*DLOG(RR)+CI                                             ADBP0248
      IF(FR.GT.0.0D0) THEN                                              ADBP0249
        RU=RR                                                           ADBP0250
      ELSE                                                              ADBP0251
        RL=RR                                                           ADBP0252
      ENDIF                                                             ADBP0253
      IF(RU-RL.GT.1.0D-15*RR) GO TO 2                                   ADBP0254
      R(I)=RR                                                           ADBP0255
    3 CONTINUE                                                          ADBP0256
      RETURN                                                            ADBP0257
      END                                                               ADBP0258
C  **************************************************************       ADBP0259
C                        SUBROUTINE ERRSPL                              ADBP0260
C  **************************************************************       ADBP0261
      SUBROUTINE ERRSPL(ERR,X,Y,N)                                      ADBP0262
C                                                                       ADBP0263
C     THIS SUBROUTINE ESTIMATES THE ERROR INTRODUCED BY NATURAL         ADBP0264
C  CUBIC SPLINE INTERPOLATION IN A TABLE X(I),Y(I) (I=1,...,N).         ADBP0265
C  THE INTERPOLATION ERROR IN THE VICINITY OF X(K) IS APPROXIMA-        ADBP0266
C  TED BY THE DIFFERENCE BETWEEN Y(K) AND THE VALUE OBTAINED FROM       ADBP0267
C  THE SPLINE THAT INTERPOLATES THE TABLE WITH THE K-TH POINT RE-       ADBP0268
C  MOVED. ERR IS THE LARGEST RELATIVE ERROR ALONG THE TABLE.            ADBP0269
C                                                                       ADBP0270
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0271
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP0272
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP0273
      COMMON/STORE/F(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)              ADBP0274
      DIMENSION X(NDIM),Y(NDIM)                                         ADBP0275
      ERR=0.0D0                                                         ADBP0276
      N1=N-1                                                            ADBP0277
      DO 2 I=2,N1                                                       ADBP0278
      DO 1 J=1,N1                                                       ADBP0279
      IF(J.LT.I) THEN                                                   ADBP0280
        R(J)=X(J)                                                       ADBP0281
        F(J)=Y(J)                                                       ADBP0282
      ELSE                                                              ADBP0283
        R(J)=X(J+1)                                                     ADBP0284
        F(J)=Y(J+1)                                                     ADBP0285
      ENDIF                                                             ADBP0286
    1 CONTINUE                                                          ADBP0287
      CALL SPLINE(R,F,A,B,C,D,0.0D0,0.0D0,N1)                           ADBP0288
      RC=X(I)                                                           ADBP0289
      YI=A(I-1)+RC*(B(I-1)+RC*(C(I-1)+RC*D(I-1)))                       ADBP0290
      IF(DABS(Y(I)).GT.1.0D-3) THEN                                     ADBP0291
        ERRP=1.0D0-YI/Y(I)                                              ADBP0292
      ELSE                                                              ADBP0293
        ERRP=YI-Y(I)                                                    ADBP0294
      ENDIF                                                             ADBP0295
      ERR=DMAX1(ERR,DABS(ERRP))                                         ADBP0296
C     WRITE(6,'(1X,I3,1P,3E18.10,2E9.1)') I,X(I),Y(I),YI,ERRP,ERR       ADBP0297
    2 CONTINUE                                                          ADBP0298
      RETURN                                                            ADBP0299
      END                                                               ADBP0300

c**==    RADIAL.FOR                                                      ADBP0301
C                                                                       ADBP0302
C                                                                       ADBP0303
C           *****************************************                   ADBP0304
C           ***     SUBROUTINE PACKAGE RADIAL     ***                   ADBP0305
C           *****************************************                   ADBP0306
C                                                                       ADBP0307
C                                       BARCELONA, DECEMBER 1994.       ADBP0308
C                                                                       ADBP0309
C  NUMERICAL SOLUTION OF THE SCHRODINGER (S) AND DIRAC (D) RADIAL       ADBP0310
C  WAVE EQUATIONS. CUBIC SPLINE FIELD + POWER SERIES METHOD.            ADBP0311
C                                                                       ADBP0312
C  IT IS ASSUMED THAT THE POTENTIAL ENERGY V(R) IS SUCH THAT THE        ADBP0313
C  FUNCTION R*V(R) IS FINITE FOR ALL R AND TENDS TO CONSTANT VA-        ADBP0314
C  LUES WHEN R TENDS TO ZERO AND TO INFINITY.                           ADBP0315
C                                                                       ADBP0316
C****   ALL QUANTITIES ARE IN ATOMIC HARTREE UNITS.                     ADBP0317
C  FOR ELECTRONS AND POSITRONS                                          ADBP0318
C    UNIT OF LENGTH = A0 = 5.29177D-11 METRES (= BOHR RADIUS),          ADBP0319
C    UNIT OF ENERGY = E0 = 27.2114 EV (= HARTREE ENERGY).               ADBP0320
C  FOR PARTICLES OF MASS 'M' (IN UNITS OF THE ELECTRON MASS) THE        ADBP0321
C  CORRESPONDING UNITS ARE                                              ADBP0322
C    UNIT OF LENGTH = A0/M,                                             ADBP0323
C    UNIT OF ENERGY = M*E0.                                             ADBP0324
C                                                                       ADBP0325
C                                                                       ADBP0326
C     THE CALLING SEQUENCE FROM THE MAIN PROGRAM IS:                    ADBP0327
C                                                                       ADBP0328
C****   CALL VINT(R,RV,NV)                                              ADBP0329
C                                                                       ADBP0330
C  THIS IS AN INITIALIZATION ROUTINE. IT DETERMINES THE NATURAL         ADBP0331
C  CUBIC SPLINE THAT INTERPOLATES THE TABLE OF VALUES OF THE            ADBP0332
C  FUNCTION R*V(R) PROVIDED BY THE USER.                                ADBP0333
C   INPUT ARGUMENTS:                                                    ADBP0334
C     R(I) ..... INPUT POTENTIAL GRID POINTS (REPEATED VALUES ARE       ADBP0335
C                INTERPRETED AS DISCONTINUITIES).                       ADBP0336
C     RV(I) .... R(I) TIMES THE POTENTIAL ENERGY AT R=R(I).             ADBP0337
C     NV ....... NUMBER OF POINTS IN THE TABLE (.LE.NDIM).              ADBP0338
C  THE R(I) GRID _MUST_ INCLUDE THE ORIGIN (R=0), AND EXTEND UP         ADBP0339
C  TO RADIAL DISTANCES FOR WHICH THE FUNCTION R*V(R) REACHES ITS        ADBP0340
C  (CONSTANT) ASYMPTOTIC VALUE. THE INPUT GRID POINTS MUST BE IN        ADBP0341
C  NON-DECREASING ORDER, I.E. R(I+1).GE.R(I).                           ADBP0342
C                                                                       ADBP0343
C****   CALL SBOUND(E,EPS,DELL,N,L) OR DBOUND(E,EPS,DELL,N,K)           ADBP0344
C                                                                       ADBP0345
C  THESE SUBROUTINES SOLVE THE RADIAL WAVE EQUATIONS FOR BOUND          ADBP0346
C  STATES.                                                              ADBP0347
C   INPUT ARGUMENTS:                                                    ADBP0348
C     E ........ ESTIMATED BINDING ENERGY (A GOOD INITIAL ESTI-         ADBP0349
C                MATE SPEEDS UP THE CALCULATION).                       ADBP0350
C     EPS ...... GLOBAL TOLERANCE, I.E. ALLOWED RELATIVE ERROR IN       ADBP0351
C                THE SUMMATION OF THE RADIAL FUNCTION SERIES.           ADBP0352
C     DELL ..... EIGENVALUE TOLERANCE, I.E. THE RELATIVE ERROR OF       ADBP0353
C                THE COMPUTED EIGENVALUE IS REQUIRED TO BE LESS         ADBP0354
C                THAN DELL (A CONVENIENT VALUE IS DELL=EPS).            ADBP0355
C     N ........ PRINCIPAL QUANTUM NUMBER.                              ADBP0356
C     L ........ ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER.               ADBP0357
C     K ........ RELATIVISTIC ANGULAR MOMENTUM QUANTUM NUMBER.          ADBP0358
C                (NOTE: 0.LE.L.LE.N-1, -N.LE.K.LE.N-1)                  ADBP0359
C   OUTPUT ARGUMENT:                                                    ADBP0360
C     E ........ BINDING ENERGY.                                        ADBP0361
C                                                                       ADBP0362
C****   CALL SFREE(E,EPS,PHASE,L) OR DFREE(E,EPS,PHASE,K)               ADBP0363
C                                                                       ADBP0364
C  THESE SUBROUTINES SOLVE THE RADIAL WAVE EQUATIONS FOR FREE           ADBP0365
C  STATES.                                                              ADBP0366
C   INPUT ARGUMENTS:                                                    ADBP0367
C     E ........ KINETIC ENERGY.                                        ADBP0368
C     EPS ...... GLOBAL TOLERANCE, I.E. ALLOWED RELATIVE ERROR IN       ADBP0369
C                THE SUMMATION OF THE RADIAL FUNCTION SERIES.           ADBP0370
C     L ........ ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER.               ADBP0371
C     K ........ RELATIVISTIC ANGULAR MOMENTUM QUANTUM NUMBER.          ADBP0372
C                (NOTE: L.GE.0, K.NE.0)                                 ADBP0373
C   OUTPUT ARGUMENTS:                                                   ADBP0374
C     PHASE .... INNER PHASE SHIFT (IN RADIANS), CAUSED BY THE          ADBP0375
C                SHORT RANGE COMPONENT OF THE FIELD. FOR MODIFIED       ADBP0376
C                COULOMB FIELDS, WE HAVE                                ADBP0377
C                       TOTAL PHASE SHIFT = PHASE + DELTA               ADBP0378
C                WHERE DELTA IS THE COULOMB PHASE SHIFT (DELIVER-       ADBP0379
C                ED THROUGH THE COMMON BLOCK                            ADBP0380
C                       COMMON/OCOUL/RK,ETA,DELTA).                     ADBP0381
C                                                                       ADBP0382
C                                                                       ADBP0383
C  THE VALUES OF THE RADIAL FUNCTIONS ARE DELIVERED THROUGH THE         ADBP0384
C  COMMON BLOCK                                                         ADBP0385
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER              ADBP0386
C  WITH NDIM=2400 (CHANGE THE VALUE OF THIS PARAMETER IN _ALL_           ADBP0387
C  SUBROUTINES IF A LARGER NUMBER OF GRID POINTS IS NEEDED). THE        ADBP0388
C  GRID OF POINTS RAD(I), WHERE THE RADIAL WAVE FUNCTIONS ARE           ADBP0389
C  TABULATED, CAN BE ARBITRARILY SELECTED BY THE USER. THE QUAN-        ADBP0390
C  TITIES                                                               ADBP0391
C     RAD(I) ... USER'S RADIAL GRID,                                    ADBP0392
C     NGP ...... NUMBER OF GRID POINTS (.LE.NDIM),                      ADBP0393
C  MUST BE DEFINED BEFORE CALLING THE SOLUTION SUBROUTINES. THE         ADBP0394
C  RADIAL POINTS ARE SORTED IN INCREASING ORDER BY THE SOLUTION         ADBP0395
C  SUBROUTINES (TO AVOID CONFUSION, IT IS ADVISABLE TO ORDER THE        ADBP0396
C  INPUT GRID IN THIS WAY).                                             ADBP0397
C  THE OUTPUT QUANTITIES ARE                                            ADBP0398
C     P(I) ..... VALUE OF THE RADIAL FUNCTION P(R) AT THE I-TH          ADBP0399
C                GRID POINT, R=RAD(I).                                  ADBP0400
C     Q(I) ..... VALUE OF THE RADIAL FUNCTION Q(R) AT THE I-TH          ADBP0401
C                GRID POINT (= P'(R) FOR SCHRODINGER PARTICLES).        ADBP0402
C     ILAST .... *** BOUND STATES: FOR R.GT.RAD(ILAST), P(R) AND        ADBP0403
C                Q(R) ARE SET EQUAL TO 0.0D0.                           ADBP0404
C                *** FREE STATES: FOR R.GT.RAD(ILAST), P(R) AND         ADBP0405
C                Q(R) ARE OBTAINED IN TERMS OF THE REGULAR AND          ADBP0406
C                IRREGULAR ASYMPTOTIC COULOMB FUNCTIONS AS              ADBP0407
C                  P(R)=DCOS(PHASE)*FU(R)+DSIN(PHASE)*GU(R)             ADBP0408
C                  Q(R)=DCOS(PHASE)*FL(R)+DSIN(PHASE)*GL(R)             ADBP0409
C                WHERE FU, GU AND FL, GL ARE GIVEN BY SUBROUTINES       ADBP0410
C                SCOUL AND DCOUL WITH Z=RV(NV). WHEN THE ABSOLUTE       ADBP0411
C                VALUE OF RV(NV) IS LESS THAN EPS, Z IS SET EQUAL       ADBP0412
C                TO ZERO, SO THAT THE FUNCTIONS FU, GU, FL AND GL       ADBP0413
C                REDUCE TO SPHERICAL BESSEL FUNCTIONS OF INTEGER        ADBP0414
C                ORDER (GIVEN BY FUNCTION BESJN).                       ADBP0415
C     IER ...... ERROR CODE. A VALUE LARGER THAN ZERO INDICATES         ADBP0416
C                THAT SOME FATAL ERROR HAS BEEN FOUND DURING THE        ADBP0417
C                CALCULATION.                                           ADBP0418
C                                                                       ADBP0419
C                                                                       ADBP0420
C  BOUND STATE WAVE FUNCTIONS ARE NORMALIZED TO UNITY. THE ADOPT-       ADBP0421
C  ED NORMALIZATION FOR FREE STATES IS SUCH THAT P(R) OSCILLATES        ADBP0422
C  WITH UNIT AMPLITUDE IN THE ASYMPTOTIC REGION.                        ADBP0423
C                                                                       ADBP0424
C                                                                       ADBP0425
C****   ERROR CODES (AND TENTATIVE SOLUTIONS...):                       ADBP0426
C    0 .... EVERYTHING IS O.K.                                          ADBP0427
C    1 .... EMIN.GE.0 IN 'BOUND' (USE A DENSER GRID. IF THE ERROR       ADBP0428
C           PERSISTS THEN PROBABLY SUCH A BOUND STATE DOES NOT          ADBP0429
C           EXIST).                                                     ADBP0430
C    2 .... E=0 IN 'BOUND' (PROBABLY THIS BOUND STATE DOES NOT          ADBP0431
C           EXIST).                                                     ADBP0432
C    3 .... RAD(NGP) TOO SMALL IN 'BOUND' (EXTEND THE GRID TO           ADBP0433
C           LARGER RADII).                                              ADBP0434
C    4 .... SEVERAL ZEROS OF P(R) IN A SINGLE INTERVAL IN 'BOUND'       ADBP0435
C           (USE A DENSER GRID).                                        ADBP0436
C    5 .... E OUT OF RANGE IN 'BOUND' (ACCUMULATED ROUND-OFF            ADBP0437
C           ERRORS?).                                                   ADBP0438
C    6 .... RV(NGP)<<0 OR E>0 IN 'BOUND' (CHECK THE INPUT POTEN-        ADBP0439
C           TIAL VALUES).                                               ADBP0440
C    7 .... E.LE.0 IN 'FREE'.                                           ADBP0441
C    8 .... RAD(NGP) TOO SMALL IN 'FREE' (EXTEND THE GRID TO            ADBP0442
C           LARGER RADII).                                              ADBP0443
C                                                                       ADBP0444
C  THE PROGRAM STOPS WHEN THE INPUT QUANTUM NUMBERS ARE OUT OF          ADBP0445
C  RANGE.                                                               ADBP0446
C                                                                       ADBP0447
C                                                                       ADBP0448
C****  YOU CAN GET A PRINTED MANUAL OF THE SUBROUTINE PACKAGE BY        ADBP0449
C  RUNNING THE INPUT FILE 'RADIAL.TEX' THROUGH LaTeX, AND PRINT-        ADBP0450
C  ING THE OUTPUT FILE. BRACKETED NUMBERS IN COMMENT LINES OF THE       ADBP0451
C  PRESENT FILE INDICATE THE CORRESPONDING EQUATIONS IN THAT MA-        ADBP0452
C  NUAL.                                                                ADBP0453
C                                                                       ADBP0454
C  **************************************************************       ADBP0455
C                        SUBROUTINE VINT                                ADBP0456
C  **************************************************************       ADBP0457
      SUBROUTINE VINT(R,RV,NV)                                          ADBP0458
C                                                                       ADBP0459
C     NATURAL CUBIC SPLINE INTERPOLATION FOR R*V(R) FROM THE            ADBP0460
C  INPUT RADII AND POTENTIAL VALUES (128).                              ADBP0461
C                                                                       ADBP0462
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0463
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP0464
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),       ADBP0465
     1             VD(NPPG),NVT                                         ADBP0466
      COMMON/RGRID/X(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP0467
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)              ADBP0468
      DIMENSION R(NDIM),RV(NDIM)                                        ADBP0469
C                                                                       ADBP0470
      IF(R(1).LT.0.0D0) THEN                                            ADBP0471
      WRITE(6,2101)                                                     ADBP0472
 2101 FORMAT(1X,'*** ERROR IN VINT: R(1).LT.0.')                        ADBP0473
      STOP                                                              ADBP0474
      ENDIF                                                             ADBP0475
      IF(NV.GT.NDIM) THEN                                               ADBP0476
      WRITE(6,2102) NDIM                                                ADBP0477
 2102 FORMAT(1X,'*** ERROR IN VINT: INPUT POTENTIAL GRID WITH ',        ADBP0478
     1  'MORE THAN ',I5,' DATA POINTS.')                                ADBP0479
      STOP                                                              ADBP0480
      ENDIF                                                             ADBP0481
      R(1)=0.0D0                                                        ADBP0482
C                                                                       ADBP0483
      IO=0                                                              ADBP0484
      I=0                                                               ADBP0485
      K=0                                                               ADBP0486
    1 I=I+1                                                             ADBP0487
      K=K+1                                                             ADBP0488
      X(K)=R(I)                                                         ADBP0489
      Y(K)=RV(I)                                                        ADBP0490
      IF(I.EQ.NV) GO TO 2                                               ADBP0491
      IF(R(I).LT.R(I+1)-1.0D-12) GO TO 1                                ADBP0492
    2 CONTINUE                                                          ADBP0493
C                                                                       ADBP0494
      CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)                            ADBP0495
C                                                                       ADBP0496
      K=K-1                                                             ADBP0497
      DO 3 J=1,K                                                        ADBP0498
      IO=IO+1                                                           ADBP0499
      RG(IO)=X(J)                                                       ADBP0500
      RVG(IO)=Y(J)                                                      ADBP0501
      VA(IO)=A(J)                                                       ADBP0502
      VB(IO)=B(J)                                                       ADBP0503
      VC(IO)=C(J)                                                       ADBP0504
      VD(IO)=D(J)                                                       ADBP0505
    3 CONTINUE                                                          ADBP0506
      IF(I.LT.NV) THEN                                                  ADBP0507
        K=0                                                             ADBP0508
        GO TO 1                                                         ADBP0509
      ENDIF                                                             ADBP0510
C  ****  AN EXTRA POINT IS ADDED TO THE GRID, AND R*V(R) IS SET         ADBP0511
C        EQUAL TO RV(NV) FOR R.GE.R(NV)                                 ADBP0512
      IO=IO+1                                                           ADBP0513
      RG(IO)=X(K+1)                                                     ADBP0514
      RVG(IO)=Y(K+1)                                                    ADBP0515
      VA(IO)=RVG(IO)                                                    ADBP0516
      VB(IO)=0.0D0                                                      ADBP0517
      VC(IO)=0.0D0                                                      ADBP0518
      VD(IO)=0.0D0                                                      ADBP0519
      NVT=IO+1                                                          ADBP0520
      RG(NVT)=2.0D0*RG(IO)                                              ADBP0521
      RVG(NVT)=RVG(IO)                                                  ADBP0522
      VA(NVT)=RVG(IO)                                                   ADBP0523
      VB(NVT)=0.0D0                                                     ADBP0524
      VC(NVT)=0.0D0                                                     ADBP0525
      VD(NVT)=0.0D0                                                     ADBP0526
      RETURN                                                            ADBP0527
      END                                                               ADBP0528
C  **************************************************************       ADBP0529
C                         SUBROUTINE SBOUND                             ADBP0530
C  **************************************************************       ADBP0531
      SUBROUTINE SBOUND(E,EPS,DELL,N,L)                                 ADBP0532
C                                                                       ADBP0533
C     THIS SUBROUTINE SOLVES THE SCHRODINGER EQUATION FOR BOUND         ADBP0534
C  STATES.                                                              ADBP0535
C                                                                       ADBP0536
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0537
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP0538
C  ****  SET IWR=0 TO AVOID PRINTING PARTIAL RESULTS.                   ADBP0539
      PARAMETER (IWR=1)                                                 ADBP0540
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP0541
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP0542
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP0543
     1             VD(NPPG),NVT                                         ADBP0544
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)              ADBP0545
      COMMON/NZT/NZMAX                                                  ADBP0546
C                                                                       ADBP0547
      IER=0                                                             ADBP0548
      IF(N.LT.1) THEN                                                   ADBP0549
      WRITE(6,2101)                                                     ADBP0550
 2101 FORMAT(1X,'*** ERROR IN SBOUND: N.LT.1.')                         ADBP0551
      STOP                                                              ADBP0552
      ENDIF                                                             ADBP0553
      IF(L.LT.0) THEN                                                   ADBP0554
      WRITE(6,2102)                                                     ADBP0555
 2102 FORMAT(1X,'*** ERROR IN SBOUND: L.LT.0.')                         ADBP0556
      STOP                                                              ADBP0557
      ENDIF                                                             ADBP0558
      IF(L.GE.N) THEN                                                   ADBP0559
      WRITE(6,2103)                                                     ADBP0560
 2103 FORMAT(1X,'*** ERROR IN SBOUND: L.GE.N.')                         ADBP0561
      STOP                                                              ADBP0562
      ENDIF                                                             ADBP0563
C  ****  RADIAL QUANTUM NUMBER.                                         ADBP0564
      NR=N-L-1                                                          ADBP0565
C                                                                       ADBP0566
      DELL=DABS(DELL)                                                   ADBP0567
      IF(E.GT.-1.0D-1) E=-1.0D-1                                        ADBP0568
      FL1=0.5D0*L*(L+1)                                                 ADBP0569
C                                                                       ADBP0570
C  ****  MERGE THE 'RG' AND 'RAD' GRIDS,                                ADBP0571
C                                                                       ADBP0572
      T=DMAX1(0.5D0*EPS,1.0D-10)                                        ADBP0573
      DO 100 I=1,NVT                                                    ADBP0574
      R(I)=RG(I)                                                        ADBP0575
      IND(I)=I                                                          ADBP0576
  100 CONTINUE                                                          ADBP0577
      NRT=NVT                                                           ADBP0578
      DO 200 I=1,NGP                                                    ADBP0579
      RLOC=RAD(I)                                                       ADBP0580
      CALL FINDI(RG,RLOC,NVT,J)                                         ADBP0581
      TST=DMIN1(DABS(RLOC-RG(J)),DABS(RLOC-RG(J+1)))                    ADBP0582
      IF(TST.LT.T) GO TO 200                                            ADBP0583
      NRT=NRT+1                                                         ADBP0584
      R(NRT)=RLOC                                                       ADBP0585
      IND(NRT)=J                                                        ADBP0586
  200 CONTINUE                                                          ADBP0587
C  ****  ... AND SORT THE RESULTING R-GRID IN INCREASING ORDER.         ADBP0588
      DO 400 I=1,NRT-1                                                  ADBP0589
      RMIN=1.0D35                                                       ADBP0590
      IMIN=I                                                            ADBP0591
      DO 300 J=I,NRT                                                    ADBP0592
      IF(R(J).LT.RMIN) THEN                                             ADBP0593
        RMIN=R(J)                                                       ADBP0594
        IMIN=J                                                          ADBP0595
      ENDIF                                                             ADBP0596
  300 CONTINUE                                                          ADBP0597
      IF(IMIN.NE.I) THEN                                                ADBP0598
        RMIN=R(I)                                                       ADBP0599
        R(I)=R(IMIN)                                                    ADBP0600
        R(IMIN)=RMIN                                                    ADBP0601
        INDMIN=IND(I)                                                   ADBP0602
        IND(I)=IND(IMIN)                                                ADBP0603
        IND(IMIN)=INDMIN                                                ADBP0604
      ENDIF                                                             ADBP0605
  400 CONTINUE                                                          ADBP0606
C                                                                       ADBP0607
C  ****  MINIMUM OF THE EFFECTIVE RADIAL POTENTIAL. (165)               ADBP0608
C                                                                       ADBP0609
      EMIN=1.0D0                                                        ADBP0610
      DO 1 I=2,NRT                                                      ADBP0611
      RN=R(I)                                                           ADBP0612
      J=IND(I)                                                          ADBP0613
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))                          ADBP0614
      EMIN=DMIN1(EMIN,(RVN+FL1/RN)/RN)                                  ADBP0615
    1 CONTINUE                                                          ADBP0616
      IF(EMIN.GT.-1.0D-35) THEN                                         ADBP0617
      IER=1                                                             ADBP0618
      WRITE(6,1001)                                                     ADBP0619
 1001 FORMAT(1X,'*** ERROR 1 (IN SBOUND): EMIN.GE.0.'/5X,               ADBP0620
     1'(USE A DENSER GRID. IF THE ERROR PERSISTS THEN PROBABLY',        ADBP0621
     2/6X'SUCH A BOUND STATE DOES NOT EXIST).')                         ADBP0622
      RETURN                                                            ADBP0623
      ENDIF                                                             ADBP0624
      RN=R(NRT)                                                         ADBP0625
      J=IND(NRT)                                                        ADBP0626
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))                          ADBP0627
      ESUP=(RVN+FL1/RN)/RN                                              ADBP0628
      IF(ESUP.GT.0.0D0) ESUP=0.0D0                                      ADBP0629
      IF(L.EQ.0) THEN                                                   ADBP0630
        TEHYDR=-(R(1)/N)**2                                             ADBP0631
        EMIN=DMIN1(E,TEHYDR,-10.0D0)                                    ADBP0632
      ENDIF                                                             ADBP0633
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)                    ADBP0634
      EMAX=ESUP                                                         ADBP0635
      ICMIN=0                                                           ADBP0636
      ICMAX=0                                                           ADBP0637
C                                                                       ADBP0638
    2 CONTINUE                                                          ADBP0639
      IF(E.GT.-1.0D-16) THEN                                            ADBP0640
      IER=2                                                             ADBP0641
      WRITE(6,1002)                                                     ADBP0642
 1002 FORMAT(1X,'*** ERROR 2 (IN SBOUND): E=0.'/5X,                     ADBP0643
     1'(PROBABLY THIS BOUND STATE DOES NOT EXIST).')                    ADBP0644
      RETURN                                                            ADBP0645
      ENDIF                                                             ADBP0646
C  ****  OUTER TURNING POINT. (165)                                     ADBP0647
      DO 3 K=2,NRT                                                      ADBP0648
      IOTP=NRT+2-K                                                      ADBP0649
      RN=R(IOTP)                                                        ADBP0650
      J=IND(IOTP)                                                       ADBP0651
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))                          ADBP0652
      EKIN=E-(RVN+FL1/RN)/RN                                            ADBP0653
      IF(EKIN.GT.0.0D0) GO TO 4                                         ADBP0654
    3 CONTINUE                                                          ADBP0655
    4 IOTP=IOTP+1                                                       ADBP0656
      IF(IOTP.GT.NRT-1) THEN                                            ADBP0657
      IER=3                                                             ADBP0658
      WRITE(6,1003)                                                     ADBP0659
 1003 FORMAT(1X,'*** ERROR 3 (IN SBOUND): RAD(NGP) TOO SMALL.'          ADBP0660
     1/5X,'(EXTEND THE GRID TO LARGER RADII).')                         ADBP0661
      RETURN                                                            ADBP0662
      ENDIF                                                             ADBP0663
C                                                                       ADBP0664
C  ****  OUTWARD SOLUTION.                                              ADBP0665
C                                                                       ADBP0666
      CALL SOUTW(E,EPS,L,NR,NZERO,IOTP)                                 ADBP0667
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN                               ADBP0668
      IER=4                                                             ADBP0669
      WRITE(6,1004)                                                     ADBP0670
 1004 FORMAT(1X,'*** ERROR 4 (IN SBOUND): SEVERAL ZEROS OF P(R)',       ADBP0671
     1/5X,'IN A SINGLE INTERVAL (USE A DENSER GRID).')                  ADBP0672
      RETURN                                                            ADBP0673
      ENDIF                                                             ADBP0674
C                                                                       ADBP0675
C  ****  TOO MANY NODES.                                                ADBP0676
      IF(NZERO.GT.NR) THEN                                              ADBP0677
      IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)                        ADBP0678
      EMAX=E                                                            ADBP0679
      ICMAX=1                                                           ADBP0680
      E=0.5D0*(EMIN+E)                                                  ADBP0681
      IF(IWR.NE.0) THEN                                                 ADBP0682
      WRITE(6,2000) N,L                                                 ADBP0683
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP0684
      WRITE(6,2002) E                                                   ADBP0685
      WRITE(6,2004) EMIN,EMAX                                           ADBP0686
      ENDIF                                                             ADBP0687
C                                                                       ADBP0688
      IF(EMAX-EMIN.LT.DELL*DABS(EMIN)) THEN                             ADBP0689
      IER=5                                                             ADBP0690
      WRITE(6,1005)                                                     ADBP0691
      RETURN                                                            ADBP0692
      ENDIF                                                             ADBP0693
C                                                                       ADBP0694
      GO TO 2                                                           ADBP0695
      ENDIF                                                             ADBP0696
C  ****  TOO FEW NODES.                                                 ADBP0697
      IF(NZERO.LT.NR) THEN                                              ADBP0698
      ICMIN=1                                                           ADBP0699
      EMIN=E                                                            ADBP0700
      E=0.5D0*(E+EMAX)                                                  ADBP0701
      IF(IWR.NE.0) THEN                                                 ADBP0702
      WRITE(6,2000) N,L                                                 ADBP0703
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP0704
      WRITE(6,2002) E                                                   ADBP0705
      WRITE(6,2004) EMIN,EMAX                                           ADBP0706
      ENDIF                                                             ADBP0707
C                                                                       ADBP0708
      IF(EMAX-EMIN.LT.DELL*DABS(EMIN)) THEN                             ADBP0709
      IER=5                                                             ADBP0710
      WRITE(6,1005)                                                     ADBP0711
      RETURN                                                            ADBP0712
      ENDIF                                                             ADBP0713
C                                                                       ADBP0714
      GO TO 2                                                           ADBP0715
      ENDIF                                                             ADBP0716
C  ****  THE CORRECT NUMBER OF NODES HAS BEEN FOUND.                    ADBP0717
      PO=P(IOTP)                                                        ADBP0718
      QO=Q(IOTP)                                                        ADBP0719
C                                                                       ADBP0720
C  ****  INWARD SOLUTION                                                ADBP0721
C                                                                       ADBP0722
      CALL SINW(E,EPS,L,IOTP)                                           ADBP0723
      IF(IER.GT.0) RETURN                                               ADBP0724
C  ****  MATCHING OF THE OUTWARD AND INWARD SOLUTIONS.                  ADBP0725
      FACT=PO/P(IOTP)                                                   ADBP0726
      DO 5 I=IOTP,ILAST                                                 ADBP0727
      P(I)=P(I)*FACT                                                    ADBP0728
    5 Q(I)=Q(I)*FACT                                                    ADBP0729
      QI=Q(IOTP)                                                        ADBP0730
      RLAST=R(ILAST)                                                    ADBP0731
C  ****  NORMALIZATION                                                  ADBP0732
      CALL SPLINE(R,P,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP0733
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUM,ILAST)                      ADBP0734
C  ****  EIGENVALUE CORRECTION. (175)                                   ADBP0735
      IF(SUM.LT.1.0D-15) SUM=1.0D0                                      ADBP0736
      DE=PO*(QO-QI)/(SUM+SUM)                                           ADBP0737
      EP=E+DE                                                           ADBP0738
C                                                                       ADBP0739
      IF(DE.LT.0.0D0) THEN                                              ADBP0740
      ICMAX=1                                                           ADBP0741
      EMAX=E                                                            ADBP0742
      ENDIF                                                             ADBP0743
C                                                                       ADBP0744
      IF(DE.GT.0.0D0) THEN                                              ADBP0745
      ICMIN=1                                                           ADBP0746
      EMIN=E                                                            ADBP0747
      ENDIF                                                             ADBP0748
C                                                                       ADBP0749
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN                                ADBP0750
      EMIN=1.1D0*EMIN                                                   ADBP0751
      IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)                                  ADBP0752
      ENDIF                                                             ADBP0753
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)                   ADBP0754
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)                   ADBP0755
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)                                  ADBP0756
C                                                                       ADBP0757
      IF(IWR.NE.0) THEN                                                 ADBP0758
      WRITE(6,2000) N,L                                                 ADBP0759
 2000 FORMAT(/2X,'SUBROUTINE SBOUND.   N =',I3,'   L =',I3)             ADBP0760
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP0761
 2001 FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,              ADBP0762
     1'   NGP =',I5)                                                    ADBP0763
      WRITE(6,2002) EP                                                  ADBP0764
 2002 FORMAT(2X,'E NEW = ',1P,D22.15)                                   ADBP0765
      WRITE(6,2003) E,DE                                                ADBP0766
 2003 FORMAT(2X,'E OLD = ',1P,D22.15,'   DE = ',D11.4)                  ADBP0767
      WRITE(6,2004) EMIN,EMAX                                           ADBP0768
 2004 FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)                  ADBP0769
      ENDIF                                                             ADBP0770
C                                                                       ADBP0771
      IF(EP.GE.ESUP.AND.DABS(E-ESUP).LT.DELL*DABS(ESUP)) THEN           ADBP0772
      IER=5                                                             ADBP0773
      WRITE(6,1005)                                                     ADBP0774
 1005 FORMAT(1X,'*** ERROR 5 (IN SBOUND): E OUT OF RANGE.'/5X,          ADBP0775
     1'(ACCUMULATED ROUND-OFF ERRORS?).')                               ADBP0776
      RETURN                                                            ADBP0777
      ENDIF                                                             ADBP0778
      EO=E                                                              ADBP0779
      E=EP                                                              ADBP0780
      IF(DMIN1(DABS(DE),DABS(E-EO)).GT.DABS(E*DELL)) GO TO 2            ADBP0781
C  ****  NORMALIZATION                                                  ADBP0782
      FACT=1.0D0/DSQRT(SUM)                                             ADBP0783
      DO 6 I=1,ILAST                                                    ADBP0784
      P(I)=P(I)*FACT                                                    ADBP0785
    6 Q(I)=Q(I)*FACT                                                    ADBP0786
C                                                                       ADBP0787
C  ****  EXTRACT THE 'RAD' GRID...                                      ADBP0788
C                                                                       ADBP0789
      DO 500 I=1,NGP                                                    ADBP0790
      RLOC=RAD(I)                                                       ADBP0791
      CALL FINDI(R,RLOC,NRT,J)                                          ADBP0792
      IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN                                 ADBP0793
        PIO(I)=P(J)                                                     ADBP0794
        QIO(I)=Q(J)                                                     ADBP0795
      ELSE                                                              ADBP0796
        PIO(I)=P(J+1)                                                   ADBP0797
        QIO(I)=Q(J+1)                                                   ADBP0798
      ENDIF                                                             ADBP0799
  500 CONTINUE                                                          ADBP0800
C                                                                       ADBP0801
      IF(DABS(P(NRT)).GT.1.0D-5*DABS(P(IOTP))) THEN                     ADBP0802
      IER=3                                                             ADBP0803
      WRITE(6,1003)                                                     ADBP0804
      ENDIF                                                             ADBP0805
C                                                                       ADBP0806
      RLOC=R(ILAST)                                                     ADBP0807
      CALL FINDI(RAD,RLOC,NGP,ILAST)                                    ADBP0808
      ILAST=ILAST+1                                                     ADBP0809
      RETURN                                                            ADBP0810
      END                                                               ADBP0811
C  **************************************************************       ADBP0812
C                         SUBROUTINE DBOUND                             ADBP0813
C  **************************************************************       ADBP0814
      SUBROUTINE DBOUND(E,EPS,DELL,N,K)                                 ADBP0815
C                                                                       ADBP0816
C     THIS SUBROUTINE SOLVES THE DIRAC EQUATION FOR BOUND STATES.       ADBP0817
C                                                                       ADBP0818
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP0819
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG,                   ADBP0820
     1  SL=137.036D0)                                                   ADBP0821
C  ****  SET IWR=0 TO AVOID PRINTING PARTIAL RESULTS.                   ADBP0822
      PARAMETER (IWR=1)                                                 ADBP0823
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP0824
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP0825
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP0826
     1             VD(NPPG),NVT                                         ADBP0827
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)              ADBP0828
      COMMON/NZT/NZMAX                                                  ADBP0829
C                                                                       ADBP0830
      IER=0                                                             ADBP0831
      IF(N.LT.1) THEN                                                   ADBP0832
      WRITE(6,2101)                                                     ADBP0833
 2101 FORMAT(1X,'*** ERROR IN DBOUND: N.LT.1.')                         ADBP0834
      STOP                                                              ADBP0835
      ENDIF                                                             ADBP0836
      IF(K.EQ.0) THEN                                                   ADBP0837
      WRITE(6,2102)                                                     ADBP0838
 2102 FORMAT(1X,'*** ERROR IN DBOUND: K.EQ.0.')                         ADBP0839
      STOP                                                              ADBP0840
      ENDIF                                                             ADBP0841
      IF(K.LT.-N) THEN                                                  ADBP0842
      WRITE(6,2103)                                                     ADBP0843
 2103 FORMAT(1X,'*** ERROR IN DBOUND: K.LT.-N.')                        ADBP0844
      STOP                                                              ADBP0845
      ENDIF                                                             ADBP0846
      IF(K.GE.N) THEN                                                   ADBP0847
      WRITE(6,2104)                                                     ADBP0848
 2104 FORMAT(1X,'*** ERROR IN DBOUND: K.GE.N.')                         ADBP0849
      STOP                                                              ADBP0850
      ENDIF                                                             ADBP0851
C  ****  ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER. (15)                  ADBP0852
      IF(K.LT.0) THEN                                                   ADBP0853
      L=-K-1                                                            ADBP0854
      ELSE                                                              ADBP0855
      L=K                                                               ADBP0856
      ENDIF                                                             ADBP0857
C  ****  RADIAL QUANTUM NUMBER.                                         ADBP0858
      NR=N-L-1                                                          ADBP0859
C                                                                       ADBP0860
      DELL=DABS(DELL)                                                   ADBP0861
      IF(E.GT.-1.0D-1) E=-1.0D-1                                        ADBP0862
      FL1=0.5D0*L*(L+1)                                                 ADBP0863
C                                                                       ADBP0864
C  ****  MERGE THE 'RG' AND 'RAD' GRIDS,                                ADBP0865
C                                                                       ADBP0866
      T=DMAX1(0.5D0*EPS,1.0D-10)                                        ADBP0867
      DO 100 I=1,NVT                                                    ADBP0868
      R(I)=RG(I)                                                        ADBP0869
      IND(I)=I                                                          ADBP0870
  100 CONTINUE                                                          ADBP0871
      NRT=NVT                                                           ADBP0872
      DO 200 I=1,NGP                                                    ADBP0873
      RLOC=RAD(I)                                                       ADBP0874
      CALL FINDI(RG,RLOC,NVT,J)                                         ADBP0875
      TST=DMIN1(DABS(RLOC-RG(J)),DABS(RLOC-RG(J+1)))                    ADBP0876
      IF(TST.LT.T) GO TO 200                                            ADBP0877
      NRT=NRT+1                                                         ADBP0878
      R(NRT)=RLOC                                                       ADBP0879
      IND(NRT)=J                                                        ADBP0880
  200 CONTINUE                                                          ADBP0881
C  ****  ... AND SORT THE RESULTING R-GRID IN INCREASING ORDER.         ADBP0882
      DO 400 I=1,NRT-1                                                  ADBP0883
      RMIN=1.0D35                                                       ADBP0884
      IMIN=I                                                            ADBP0885
      DO 300 J=I,NRT                                                    ADBP0886
      IF(R(J).LT.RMIN) THEN                                             ADBP0887
        RMIN=R(J)                                                       ADBP0888
        IMIN=J                                                          ADBP0889
      ENDIF                                                             ADBP0890
  300 CONTINUE                                                          ADBP0891
      IF(IMIN.NE.I) THEN                                                ADBP0892
        RMIN=R(I)                                                       ADBP0893
        R(I)=R(IMIN)                                                    ADBP0894
        R(IMIN)=RMIN                                                    ADBP0895
        INDMIN=IND(I)                                                   ADBP0896
        IND(I)=IND(IMIN)                                                ADBP0897
        IND(IMIN)=INDMIN                                                ADBP0898
      ENDIF                                                             ADBP0899
  400 CONTINUE                                                          ADBP0900
C                                                                       ADBP0901
C  ****  MINIMUM OF THE EFFECTIVE RADIAL POTENTIAL. (165)               ADBP0902
C                                                                       ADBP0903
      EMIN=1.0D0                                                        ADBP0904
      DO 1 I=2,NRT                                                      ADBP0905
      RN=R(I)                                                           ADBP0906
      J=IND(I)                                                          ADBP0907
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))                          ADBP0908
      EMIN=DMIN1(EMIN,(RVN+FL1/RN)/RN)                                  ADBP0909
    1 CONTINUE                                                          ADBP0910
      IF(EMIN.GT.-1.0D-35) THEN                                         ADBP0911
      IER=1                                                             ADBP0912
      WRITE(6,1001)                                                     ADBP0913
 1001 FORMAT(1X,'*** ERROR 1 (IN DBOUND): EMIN.GE.0.'/5X,               ADBP0914
     1'(USE A DENSER GRID. IF THE ERROR PERSISTS THEN PROBABLY',        ADBP0915
     2/6X'SUCH A BOUND STATE DOES NOT EXIST).')                         ADBP0916
      RETURN                                                            ADBP0917
      ENDIF                                                             ADBP0918
      RN=R(NRT)                                                         ADBP0919
      J=IND(NRT)                                                        ADBP0920
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))                          ADBP0921
      ESUP=(RVN+FL1/RN)/RN                                              ADBP0922
      IF(ESUP.GT.0.0D0) ESUP=0.0D0                                      ADBP0923
      IF(L.EQ.0) THEN                                                   ADBP0924
        TEHYDR=-(R(1)/N)**2                                             ADBP0925
        EMIN=DMIN1(E,TEHYDR,-10.0D0)                                    ADBP0926
      ENDIF                                                             ADBP0927
      IF(EMIN.LT.-SL*SL) EMIN=-SL*SL                                    ADBP0928
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)                    ADBP0929
      EMAX=ESUP                                                         ADBP0930
      ICMIN=0                                                           ADBP0931
      ICMAX=0                                                           ADBP0932
C                                                                       ADBP0933
    2 CONTINUE                                                          ADBP0934
      IF(E.GT.-1.0D-16) THEN                                            ADBP0935
      IER=2                                                             ADBP0936
      WRITE(6,1002)                                                     ADBP0937
 1002 FORMAT(1X,'*** ERROR 2 (IN DBOUND): E=0.'/5X,                     ADBP0938
     1'(PROBABLY THIS BOUND STATE DOES NOT EXIST).')                    ADBP0939
      RETURN                                                            ADBP0940
      ENDIF                                                             ADBP0941
C  ****  OUTER TURNING POINT. (165)                                     ADBP0942
      DO 3 J=2,NRT                                                      ADBP0943
      IOTP=NRT+2-J                                                      ADBP0944
      RN=R(IOTP)                                                        ADBP0945
      JJ=IND(IOTP)                                                      ADBP0946
      RVN=VA(JJ)+RN*(VB(JJ)+RN*(VC(JJ)+RN*VD(JJ)))                      ADBP0947
      EKIN=E-(RVN+FL1/RN)/RN                                            ADBP0948
      IF(EKIN.GT.0.0D0) GO TO 4                                         ADBP0949
    3 CONTINUE                                                          ADBP0950
    4 IOTP=IOTP+1                                                       ADBP0951
      IF(IOTP.GT.NRT-1) THEN                                            ADBP0952
      IER=3                                                             ADBP0953
      WRITE(6,1003)                                                     ADBP0954
 1003 FORMAT(1X,'*** ERROR 3 (IN DBOUND): RAD(NGP) TOO SMALL.'          ADBP0955
     1/5X,'(EXTEND THE GRID TO LARGER RADII).')                         ADBP0956
      RETURN                                                            ADBP0957
      ENDIF                                                             ADBP0958
C                                                                       ADBP0959
C  ****  OUTWARD SOLUTION.                                              ADBP0960
C                                                                       ADBP0961
      CALL DOUTW(E,EPS,K,NR,NZERO,IOTP)                                 ADBP0962
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN                               ADBP0963
      IER=4                                                             ADBP0964
      WRITE(6,1004)                                                     ADBP0965
 1004 FORMAT(1X,'*** ERROR 4 (IN DBOUND): SEVERAL ZEROS OF P(R)',       ADBP0966
     1/5X,'IN A SINGLE INTERVAL (USE A DENSER GRID).')                  ADBP0967
      RETURN                                                            ADBP0968
      ENDIF                                                             ADBP0969
C                                                                       ADBP0970
C  ****  TOO MANY NODES.                                                ADBP0971
      IF(NZERO.GT.NR) THEN                                              ADBP0972
      IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)                        ADBP0973
      EMAX=E                                                            ADBP0974
      ICMAX=1                                                           ADBP0975
      E=0.5D0*(EMIN+E)                                                  ADBP0976
      IF(IWR.NE.0) THEN                                                 ADBP0977
      WRITE(6,2000) N,K                                                 ADBP0978
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP0979
      WRITE(6,2002) E                                                   ADBP0980
      WRITE(6,2004) EMIN,EMAX                                           ADBP0981
      ENDIF                                                             ADBP0982
C                                                                       ADBP0983
      IF(EMAX-EMIN.LT.DELL*DABS(EMIN)) THEN                             ADBP0984
      IER=5                                                             ADBP0985
      WRITE(6,1005)                                                     ADBP0986
      RETURN                                                            ADBP0987
      ENDIF                                                             ADBP0988
C                                                                       ADBP0989
      GO TO 2                                                           ADBP0990
      ENDIF                                                             ADBP0991
C  ****  TOO FEW NODES.                                                 ADBP0992
      IF(NZERO.LT.NR) THEN                                              ADBP0993
      ICMIN=1                                                           ADBP0994
      EMIN=E                                                            ADBP0995
      E=0.5D0*(E+EMAX)                                                  ADBP0996
      IF(IWR.NE.0) THEN                                                 ADBP0997
      WRITE(6,2000) N,K                                                 ADBP0998
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP0999
      WRITE(6,2002) E                                                   ADBP1000
      WRITE(6,2004) EMIN,EMAX                                           ADBP1001
      ENDIF                                                             ADBP1002
C                                                                       ADBP1003
      IF(EMAX-EMIN.LT.DELL*DABS(EMIN)) THEN                             ADBP1004
      IER=5                                                             ADBP1005
      WRITE(6,1005)                                                     ADBP1006
      RETURN                                                            ADBP1007
      ENDIF                                                             ADBP1008
C                                                                       ADBP1009
      GO TO 2                                                           ADBP1010
      ENDIF                                                             ADBP1011
C  ****  THE CORRECT NUMBER OF NODES HAS BEEN FOUND.                    ADBP1012
      PO=P(IOTP)                                                        ADBP1013
      QO=Q(IOTP)                                                        ADBP1014
C                                                                       ADBP1015
C  ****  INWARD SOLUTION                                                ADBP1016
C                                                                       ADBP1017
      CALL DINW(E,EPS,K,IOTP)                                           ADBP1018
      IF(IER.GT.0) RETURN                                               ADBP1019
C  ****  MATCHING OF THE OUTWARD AND INWARD SOLUTIONS.                  ADBP1020
      FACT=PO/P(IOTP)                                                   ADBP1021
      DO 5 I=IOTP,NRT                                                   ADBP1022
      P(I)=P(I)*FACT                                                    ADBP1023
    5 Q(I)=Q(I)*FACT                                                    ADBP1024
      QI=Q(IOTP)                                                        ADBP1025
      RLAST=R(ILAST)                                                    ADBP1026
C  ****  NORMALIZATION                                                  ADBP1027
      CALL SPLINE(R,P,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP1028
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUMP,ILAST)                     ADBP1029
      CALL SPLINE(R,Q,A,B,C,D,0.0D0,0.0D0,ILAST)                        ADBP1030
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUMQ,ILAST)                     ADBP1031
      SUM=SUMP+SUMQ                                                     ADBP1032
C  ****  EIGENVALUE CORRECTION. (179)                                   ADBP1033
      IF(SUM.LT.1.0D-15) SUM=1.0D0                                      ADBP1034
      DE=SL*PO*(QI-QO)/SUM                                              ADBP1035
      EP=E+DE                                                           ADBP1036
C                                                                       ADBP1037
      IF(DE.LT.0.0D0) THEN                                              ADBP1038
      ICMAX=1                                                           ADBP1039
      EMAX=E                                                            ADBP1040
      ENDIF                                                             ADBP1041
C                                                                       ADBP1042
      IF(DE.GT.0.0D0) THEN                                              ADBP1043
      ICMIN=1                                                           ADBP1044
      EMIN=E                                                            ADBP1045
      ENDIF                                                             ADBP1046
C                                                                       ADBP1047
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN                                ADBP1048
      EMIN=1.1D0*EMIN                                                   ADBP1049
      IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)                                  ADBP1050
      ENDIF                                                             ADBP1051
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)                   ADBP1052
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)                   ADBP1053
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)                                  ADBP1054
C                                                                       ADBP1055
      IF(IWR.NE.0) THEN                                                 ADBP1056
      WRITE(6,2000) N,K                                                 ADBP1057
 2000 FORMAT(/2X,'SUBROUTINE DBOUND.   N =',I3,'   K =',I3)             ADBP1058
      WRITE(6,2001) NR,NZERO,IOTP,NRT                                   ADBP1059
 2001 FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,              ADBP1060
     1'   NGP =',I5)                                                    ADBP1061
      WRITE(6,2002) EP                                                  ADBP1062
 2002 FORMAT(2X,'E NEW = ',1P,D22.15)                                   ADBP1063
      WRITE(6,2003) E,DE                                                ADBP1064
 2003 FORMAT(2X,'E OLD = ',1P,D22.15,'   DE = ',D11.4)                  ADBP1065
      WRITE(6,2004) EMIN,EMAX                                           ADBP1066
 2004 FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)                  ADBP1067
      ENDIF                                                             ADBP1068
C                                                                       ADBP1069
      IF(EP.GT.ESUP.AND.DABS(E-ESUP).LT.DELL*DABS(ESUP)) THEN           ADBP1070
      IER=5                                                             ADBP1071
      WRITE(6,1005)                                                     ADBP1072
 1005 FORMAT(1X,'*** ERROR 5 (IN DBOUND): E OUT OF RANGE.'/5X,          ADBP1073
     1'(ACCUMULATED ROUND-OFF ERRORS?).')                               ADBP1074
      RETURN                                                            ADBP1075
      ENDIF                                                             ADBP1076
      EO=E                                                              ADBP1077
      E=EP                                                              ADBP1078
      IF(DMIN1(DABS(DE),DABS(E-EO)).GT.DABS(E*DELL)) GO TO 2            ADBP1079
C  ****  NORMALIZATION                                                  ADBP1080
      FACT=1.0D0/DSQRT(SUM)                                             ADBP1081
      DO 6 I=1,ILAST                                                    ADBP1082
      P(I)=P(I)*FACT                                                    ADBP1083
    6 Q(I)=Q(I)*FACT                                                    ADBP1084
C                                                                       ADBP1085
C  ****  EXTRACT THE 'RAD' GRID...                                      ADBP1086
C                                                                       ADBP1087
      DO 500 I=1,NGP                                                    ADBP1088
      RLOC=RAD(I)                                                       ADBP1089
      CALL FINDI(R,RLOC,NRT,J)                                          ADBP1090
      IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN                                 ADBP1091
        PIO(I)=P(J)                                                     ADBP1092
        QIO(I)=Q(J)                                                     ADBP1093
      ELSE                                                              ADBP1094
        PIO(I)=P(J+1)                                                   ADBP1095
        QIO(I)=Q(J+1)                                                   ADBP1096
      ENDIF                                                             ADBP1097
  500 CONTINUE                                                          ADBP1098
C                                                                       ADBP1099
      IF(DABS(P(NRT)).GT.1.0D-5*DABS(P(IOTP))) THEN                     ADBP1100
      IER=3                                                             ADBP1101
      WRITE(6,1003)                                                     ADBP1102
      ENDIF                                                             ADBP1103
C                                                                       ADBP1104
      RLOC=R(ILAST)                                                     ADBP1105
      CALL FINDI(RAD,RLOC,NGP,ILAST)                                    ADBP1106
      ILAST=ILAST+1                                                     ADBP1107
      RETURN                                                            ADBP1108
      END                                                               ADBP1109
C  **************************************************************       ADBP1110
C                         SUBROUTINE SFREE                              ADBP1111
C  **************************************************************       ADBP1112
      SUBROUTINE SFREE(E,EPS,PHASE,L)                                   ADBP1113
C                                                                       ADBP1114
C     THIS SUBROUTINE SOLVES THE SCHRODINGER EQUATION FOR FREE          ADBP1115
C  STATES.                                                              ADBP1116
C                                                                       ADBP1117
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1118
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG,                   ADBP1119
     1  PI=3.1415926535897932D0,TPI=PI+PI)                              ADBP1120
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1121
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1122
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1123
     1             VD(NPPG),NVT                                         ADBP1124
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)          ADBP1125
      COMMON/NZT/NZMAX                                                  ADBP1126
      COMMON/OCOUL/RK,ETA,DELTA                                         ADBP1127
	external besjn

      ETA=0.0D0                                                         ADBP1128
      DELTA=0.0D0                                                       ADBP1129
C                                                                       ADBP1130
      IER=0                                                             ADBP1131
      IF(L.LT.0) THEN                                                   ADBP1132
      WRITE(6,2101)                                                     ADBP1133
 2101 FORMAT(1X,'*** ERROR IN SFREE: L.LT.0.')                          ADBP1134
      STOP                                                              ADBP1135
      ENDIF                                                             ADBP1136
      FL1=0.5D0*L*(L+1)                                                 ADBP1137
C                                                                       ADBP1138
      IF(E.LE.0.0D0) THEN                                               ADBP1139
      IER=7                                                             ADBP1140
      WRITE(6,1007)                                                     ADBP1141
 1007 FORMAT(1X,'*** ERROR 7  (IN SFREE): E.LE.0.')                     ADBP1142
      RETURN                                                            ADBP1143
      ENDIF                                                             ADBP1144
      RK=DSQRT(E+E)                                                     ADBP1145
C                                                                       ADBP1146
C  ****  MERGE THE 'RG' AND 'RAD' GRIDS,                                ADBP1147
C                                                                       ADBP1148
      T=DMAX1(0.5D0*EPS,1.0D-10)                                        ADBP1149
      DO 100 I=1,NVT                                                    ADBP1150
      R(I)=RG(I)                                                        ADBP1151
      IND(I)=I                                                          ADBP1152
  100 CONTINUE                                                          ADBP1153
      NRT=NVT                                                           ADBP1154
      DO 200 I=1,NGP                                                    ADBP1155
      RLOC=RAD(I)                                                       ADBP1156
      CALL FINDI(RG,RLOC,NVT,J)                                         ADBP1157
      TST=DMIN1(DABS(RLOC-RG(J)),DABS(RLOC-RG(J+1)))                    ADBP1158
      IF(TST.LT.T) GO TO 200                                            ADBP1159
      NRT=NRT+1                                                         ADBP1160
      R(NRT)=RLOC                                                       ADBP1161
      IND(NRT)=J                                                        ADBP1162
  200 CONTINUE                                                          ADBP1163
C  ****  ... AND SORT THE RESULTING R-GRID IN INCREASING ORDER.         ADBP1164
      DO 400 I=1,NRT-1                                                  ADBP1165
      RMIN=1.0D35                                                       ADBP1166
      IMIN=I                                                            ADBP1167
      DO 300 J=I,NRT                                                    ADBP1168
      IF(R(J).LT.RMIN) THEN                                             ADBP1169
        RMIN=R(J)                                                       ADBP1170
        IMIN=J                                                          ADBP1171
      ENDIF                                                             ADBP1172
  300 CONTINUE                                                          ADBP1173
      IF(IMIN.NE.I) THEN                                                ADBP1174
        RMIN=R(I)                                                       ADBP1175
        R(I)=R(IMIN)                                                    ADBP1176
        R(IMIN)=RMIN                                                    ADBP1177
        INDMIN=IND(I)                                                   ADBP1178
        IND(I)=IND(IMIN)                                                ADBP1179
        IND(IMIN)=INDMIN                                                ADBP1180
      ENDIF                                                             ADBP1181
  400 CONTINUE                                                          ADBP1182
C                                                                       ADBP1183
C  ****  ASYMPTOTIC SOLUTION.                                           ADBP1184
C                                                                       ADBP1185
      ZINF=RV(NVT)                                                      ADBP1186
      IF(DABS(ZINF).LT.EPS) THEN                                        ADBP1187
C  ****  FINITE RANGE FIELDS.                                           ADBP1188
      WN=RK                                                             ADBP1189
      ILAST=NRT+1                                                       ADBP1190
      DO 1 I=4,NRT                                                      ADBP1191
      IL=ILAST-1                                                        ADBP1192
      RN=R(IL)                                                          ADBP1193
      INJ=IND(IL)                                                       ADBP1194
      RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))                  ADBP1195
      T=EPS*DABS(E*RN-FL1/RN)                                           ADBP1196
      X=RK*RN                                                           ADBP1197
      IF(DABS(RVN).GT.T) GO TO 2                                        ADBP1198
      BNL1=BESJN(2,L+1,X)                                               ADBP1199
      IF(DABS(BNL1).GT.100.0D0) GO TO 2                                 ADBP1200
      BNL=BESJN(2,L,X)                                                  ADBP1201
      BJL=BESJN(1,L,X)                                                  ADBP1202
      BJL1=BESJN(1,L+1,X)                                               ADBP1203
      ILAST=IL                                                          ADBP1204
      PA(ILAST)=X*BJL                                                   ADBP1205
      PB(ILAST)=-X*BNL                                                  ADBP1206
      QA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)                               ADBP1207
      QB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)                              ADBP1208
    1 CONTINUE                                                          ADBP1209
    2 CONTINUE                                                          ADBP1210
      IF(ILAST.EQ.NRT+1) THEN                                           ADBP1211
      IER=8                                                             ADBP1212
      WRITE(6,1008)                                                     ADBP1213
 1008 FORMAT(1X,'*** ERROR 8  (IN SFREE): RAD(NGP) TOO SMALL.'          ADBP1214
     1/5X,'(EXTEND THE GRID TO LARGER RADII).')                         ADBP1215
      RETURN                                                            ADBP1216
      ENDIF                                                             ADBP1217
      ELSE                                                              ADBP1218
C  ****  COULOMB FIELDS.                                                ADBP1219
      TAS=EPS*DABS(ZINF)                                                ADBP1220
      ILAST=NRT+1                                                       ADBP1221
      DO 3 I=4,NRT                                                      ADBP1222
      IL=ILAST-1                                                        ADBP1223
      RN=R(IL)                                                          ADBP1224
      INJ=IND(IL)                                                       ADBP1225
      RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))                  ADBP1226
      IF(DABS(RVN-ZINF).GT.TAS) GO TO 4                                 ADBP1227
      CALL SCOUL(ZINF,E,L,RN,P0,Q0,P1,Q1,ERR)                           ADBP1228
      IF(ERR.GT.EPS) GO TO 4                                            ADBP1229
      ILAST=IL                                                          ADBP1230
      PA(ILAST)=P0                                                      ADBP1231
      PB(ILAST)=P1                                                      ADBP1232
      QA(ILAST)=Q0                                                      ADBP1233
      QB(ILAST)=Q1                                                      ADBP1234
    3 CONTINUE                                                          ADBP1235
    4 CONTINUE                                                          ADBP1236
      IF(ILAST.EQ.NRT+1) THEN                                           ADBP1237
      IER=8                                                             ADBP1238
      WRITE(6,1008)                                                     ADBP1239
      RETURN                                                            ADBP1240
      ENDIF                                                             ADBP1241
      ENDIF                                                             ADBP1242
C                                                                       ADBP1243
C  ****  OUTWARD SOLUTION.                                              ADBP1244
C                                                                       ADBP1245
      CALL SOUTW(E,EPS,L,1,NZERO,ILAST)                                 ADBP1246
C                                                                       ADBP1247
C  ****  PHASE SHIFT. (187)                                             ADBP1248
C                                                                       ADBP1249
      IF(DABS(P(ILAST)).GT.EPS) THEN                                    ADBP1250
      RATIO=Q(ILAST)/P(ILAST)                                           ADBP1251
      PHASE=DATAN2(RATIO*PA(ILAST)-QA(ILAST),                           ADBP1252
     1              QB(ILAST)-RATIO*PB(ILAST))                          ADBP1253
      CD=DCOS(PHASE)                                                    ADBP1254
      SD=DSIN(PHASE)                                                    ADBP1255
      RNORM=(CD*PA(ILAST)+SD*PB(ILAST))/P(ILAST)                        ADBP1256
      ELSE                                                              ADBP1257
      PHASE=DATAN2(-PA(ILAST),PB(ILAST))                                ADBP1258
      CD=DCOS(PHASE)                                                    ADBP1259
      SD=DSIN(PHASE)                                                    ADBP1260
      RNORM=(CD*QA(ILAST)+SD*QB(ILAST))/Q(ILAST)                        ADBP1261
      ENDIF                                                             ADBP1262
C                                                                       ADBP1263
      IF(RNORM.LT.0.0D0) THEN                                           ADBP1264
      RNORM=-RNORM                                                      ADBP1265
      CD=-CD                                                            ADBP1266
      SD=-SD                                                            ADBP1267
      PHASE=PHASE+PI                                                    ADBP1268
      ENDIF                                                             ADBP1269
      IF(PHASE.GE.0.0D0) THEN                                           ADBP1270
        PHASE=DMOD(PHASE,TPI)                                           ADBP1271
        IF(PHASE.GT.PI) PHASE=PHASE-TPI                                 ADBP1272
      ELSE                                                              ADBP1273
        PHASE=-DMOD(-PHASE,TPI)                                         ADBP1274
        IF(PHASE.LT.-PI) PHASE=PHASE+TPI                                ADBP1275
      ENDIF                                                             ADBP1276
C  ****  NORMALIZED WAVE FUNCTION. (188)                                ADBP1277
      DO 5 I=1,ILAST                                                    ADBP1278
      P(I)=RNORM*P(I)                                                   ADBP1279
      Q(I)=RNORM*Q(I)                                                   ADBP1280
    5 CONTINUE                                                          ADBP1281
      IF(ILAST.LT.NRT) THEN                                             ADBP1282
      DO 6 I=ILAST+1,NRT                                                ADBP1283
      P(I)=CD*PA(I)+SD*PB(I)                                            ADBP1284
      Q(I)=CD*QA(I)+SD*QB(I)                                            ADBP1285
    6 CONTINUE                                                          ADBP1286
      ENDIF                                                             ADBP1287
C                                                                       ADBP1288
C  ****  EXTRACT THE 'RAD' GRID...                                      ADBP1289
C                                                                       ADBP1290
      DO 500 I=1,NGP                                                    ADBP1291
      RLOC=RAD(I)                                                       ADBP1292
      CALL FINDI(R,RLOC,NRT,J)                                          ADBP1293
      IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN                                 ADBP1294
        PIO(I)=P(J)                                                     ADBP1295
        QIO(I)=Q(J)                                                     ADBP1296
      ELSE                                                              ADBP1297
        PIO(I)=P(J+1)                                                   ADBP1298
        QIO(I)=Q(J+1)                                                   ADBP1299
      ENDIF                                                             ADBP1300
  500 CONTINUE                                                          ADBP1301
C                                                                       ADBP1302
      RLOC=R(ILAST)                                                     ADBP1303
      CALL FINDI(RAD,RLOC,NGP,ILAST)                                    ADBP1304
      ILAST=ILAST+1                                                     ADBP1305
      RETURN                                                            ADBP1306
      END                                                               ADBP1307
C  **************************************************************       ADBP1308
C                         SUBROUTINE DFREE                              ADBP1309
C  **************************************************************       ADBP1310
      SUBROUTINE DFREE(E,EPS,PHASE,K)                                   ADBP1311
C                                                                       ADBP1312
C     THIS SUBROUTINE SOLVES THE DIRAC EQUATION FOR FREE STATES.        ADBP1313
C                                                                       ADBP1314
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1315
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG,                   ADBP1316
     1  SL=137.036D0,PI=3.1415926535897932D0,TPI=PI+PI)                 ADBP1317
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1318
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1319
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1320
     1             VD(NPPG),NVT                                         ADBP1321
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)          ADBP1322
      COMMON/NZT/NZMAX                                                  ADBP1323
      COMMON/OCOUL/RK,ETA,DELTA                                         ADBP1324

	external besjn

      ETA=0.0D0                                                         ADBP1325
      DELTA=0.0D0                                                       ADBP1326
C                                                                       ADBP1327
      IER=0                                                             ADBP1328
      IF(K.EQ.0) THEN                                                   ADBP1329
      WRITE(6,2101)                                                     ADBP1330
 2101 FORMAT(1X,'*** ERROR IN DFREE: K.EQ.0.')                          ADBP1331
      STOP                                                              ADBP1332
      ENDIF                                                             ADBP1333
C                                                                       ADBP1334
      IF(E.LE.0.0D0) THEN                                               ADBP1335
      IER=7                                                             ADBP1336
      WRITE(6,1007)                                                     ADBP1337
 1007 FORMAT(1X,'*** ERROR 7  (IN DFREE): E.LE.0.')                     ADBP1338
      RETURN                                                            ADBP1339
      ENDIF                                                             ADBP1340
C  ****  ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER. (15)                  ADBP1341
      IF(K.LT.0) THEN                                                   ADBP1342
      L=-K-1                                                            ADBP1343
      KSIGN=1                                                           ADBP1344
      ELSE                                                              ADBP1345
      L=K                                                               ADBP1346
      KSIGN=-1                                                          ADBP1347
      ENDIF                                                             ADBP1348
      FL1=0.5D0*L*(L+1)                                                 ADBP1349
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL                                    ADBP1350
C                                                                       ADBP1351
C  ****  MERGE THE 'RG' AND 'RAD' GRIDS,                                ADBP1352
C                                                                       ADBP1353
      IF(NGP.GT.NDIM) THEN                                              ADBP1354
      WRITE(6,2102) NDIM                                                ADBP1355
 2102 FORMAT(1X,'*** ERROR IN DFREE: INPUT POTENTIAL GRID WITH',        ADBP1356
     1  ' MORE THAN ',I5,' DATA POINTS.')                               ADBP1357
      STOP                                                              ADBP1358
      ENDIF                                                             ADBP1359
      T=DMAX1(0.5D0*EPS,1.0D-10)                                        ADBP1360
      DO 100 I=1,NVT                                                    ADBP1361
      R(I)=RG(I)                                                        ADBP1362
      IND(I)=I                                                          ADBP1363
  100 CONTINUE                                                          ADBP1364
      NRT=NVT                                                           ADBP1365
      DO 200 I=1,NGP                                                    ADBP1366
      RLOC=RAD(I)                                                       ADBP1367
      CALL FINDI(RG,RLOC,NVT,J)                                         ADBP1368
      TST=DMIN1(DABS(RLOC-RG(J)),DABS(RLOC-RG(J+1)))                    ADBP1369
      IF(TST.LT.T) GO TO 200                                            ADBP1370
      NRT=NRT+1                                                         ADBP1371
      R(NRT)=RLOC                                                       ADBP1372
      IND(NRT)=J                                                        ADBP1373
  200 CONTINUE                                                          ADBP1374
C  ****  ... AND SORT THE RESULTING R-GRID IN INCREASING ORDER.         ADBP1375
      DO 400 I=1,NRT-1                                                  ADBP1376
      RMIN=1.0D35                                                       ADBP1377
      IMIN=I                                                            ADBP1378
      DO 300 J=I,NRT                                                    ADBP1379
      IF(R(J).LT.RMIN) THEN                                             ADBP1380
        RMIN=R(J)                                                       ADBP1381
        IMIN=J                                                          ADBP1382
      ENDIF                                                             ADBP1383
  300 CONTINUE                                                          ADBP1384
      IF(IMIN.NE.I) THEN                                                ADBP1385
        RMIN=R(I)                                                       ADBP1386
        R(I)=R(IMIN)                                                    ADBP1387
        R(IMIN)=RMIN                                                    ADBP1388
        INDMIN=IND(I)                                                   ADBP1389
        IND(I)=IND(IMIN)                                                ADBP1390
        IND(IMIN)=INDMIN                                                ADBP1391
      ENDIF                                                             ADBP1392
  400 CONTINUE                                                          ADBP1393
C                                                                       ADBP1394
C  ****  ASYMPTOTIC SOLUTION.                                           ADBP1395
C                                                                       ADBP1396
      ZINF=RV(NVT)                                                      ADBP1397
      IF(DABS(ZINF).LT.EPS) THEN                                        ADBP1398
C  ****  FINITE RANGE FIELDS.                                           ADBP1399
      WN=RK                                                             ADBP1400
      FACTOR=DSQRT(E/(E+2.0D0*SL*SL))                                   ADBP1401
      ILAST=NRT+1                                                       ADBP1402
      DO 1 I=4,NRT                                                      ADBP1403
      IL=ILAST-1                                                        ADBP1404
      RN=R(IL)                                                          ADBP1405
      INJ=IND(IL)                                                       ADBP1406
      RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))                  ADBP1407
      T=EPS*RN*DABS(E*RN-FL1/RN)                                        ADBP1408
      X=RK*RN                                                           ADBP1409
      IF(DABS(RVN).GT.T) GO TO 2                                        ADBP1410
      BNL=BESJN(2,L,X)                                                  ADBP1411
      IF(DABS(BNL).GT.100.0D0) GO TO 2                                  ADBP1412
      BNL1=BESJN(2,L+KSIGN,X)                                           ADBP1413
      IF(DABS(BNL1).GT.100.0D0) GO TO 2                                 ADBP1414
      BJL=BESJN(1,L,X)                                                  ADBP1415
      BJL1=BESJN(1,L+KSIGN,X)                                           ADBP1416
      ILAST=IL                                                          ADBP1417
      PA(ILAST)=X*BJL                                                   ADBP1418
      PB(ILAST)=-X*BNL                                                  ADBP1419
      QA(ILAST)=FACTOR*KSIGN*X*BJL1                                     ADBP1420
      QB(ILAST)=-FACTOR*KSIGN*X*BNL1                                    ADBP1421
    1 CONTINUE                                                          ADBP1422
    2 CONTINUE                                                          ADBP1423
      IF(ILAST.EQ.NRT+1) THEN                                           ADBP1424
      IER=8                                                             ADBP1425
      WRITE(6,1008)                                                     ADBP1426
 1008 FORMAT(1X,'*** ERROR 8  (IN DFREE): RAD(NGP) TOO SMALL.'          ADBP1427
     1/5X,'(EXTEND THE GRID TO LARGER RADII).')                         ADBP1428
      RETURN                                                            ADBP1429
      ENDIF                                                             ADBP1430
      ELSE                                                              ADBP1431
C  ****  COULOMB FIELDS.                                                ADBP1432
      TAS=EPS*DABS(ZINF)                                                ADBP1433
      ILAST=NRT+1                                                       ADBP1434
      DO 3 I=4,NRT                                                      ADBP1435
      IL=ILAST-1                                                        ADBP1436
      RN=R(IL)                                                          ADBP1437
      INJ=IND(IL)                                                       ADBP1438
      RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))                  ADBP1439
      IF(DABS(RVN-ZINF).GT.TAS) GO TO 4                                 ADBP1440
      CALL DCOUL(ZINF,E,K,RN,P0,Q0,P1,Q1,ERR)                           ADBP1441
      IF(ERR.GT.EPS) GO TO 4                                            ADBP1442
      ILAST=IL                                                          ADBP1443
      PA(ILAST)=P0                                                      ADBP1444
      PB(ILAST)=P1                                                      ADBP1445
      QA(ILAST)=Q0                                                      ADBP1446
      QB(ILAST)=Q1                                                      ADBP1447
    3 CONTINUE                                                          ADBP1448
    4 CONTINUE                                                          ADBP1449
      IF(ILAST.EQ.NRT+1) THEN                                           ADBP1450
      IER=8                                                             ADBP1451
      WRITE(6,1008)                                                     ADBP1452
      RETURN                                                            ADBP1453
      ENDIF                                                             ADBP1454
      ENDIF                                                             ADBP1455
C                                                                       ADBP1456
C  ****  OUTWARD SOLUTION.                                              ADBP1457
C                                                                       ADBP1458
      CALL DOUTW(E,EPS,K,1,NZERO,ILAST)                                 ADBP1459
C                                                                       ADBP1460
C  ****  PHASE SHIFT. (187)                                             ADBP1461
C                                                                       ADBP1462
      RM=R(ILAST)                                                       ADBP1463
      IL=IND(ILAST-1)                                                   ADBP1464
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))                         ADBP1465
      FG=(E-VF+2.0D0*SL*SL)/SL                                          ADBP1466
      PO=P(ILAST)                                                       ADBP1467
      POP=-K*PO/RM-FG*Q(ILAST)                                          ADBP1468
      IL=IND(ILAST)                                                     ADBP1469
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))                         ADBP1470
      FG=(E-VF+2.0D0*SL*SL)/SL                                          ADBP1471
      PIA=PA(ILAST)                                                     ADBP1472
      PIAP=-K*PIA/RM-FG*QA(ILAST)                                       ADBP1473
      PIB=PB(ILAST)                                                     ADBP1474
      PIBP=-K*PIB/RM-FG*QB(ILAST)                                       ADBP1475
C                                                                       ADBP1476
      IF(DABS(PO).GT.EPS) THEN                                          ADBP1477
      RATIO=POP/PO                                                      ADBP1478
      PHASE=DATAN2(RATIO*PIA-PIAP,PIBP-RATIO*PIB)                       ADBP1479
      CD=DCOS(PHASE)                                                    ADBP1480
      SD=DSIN(PHASE)                                                    ADBP1481
      RNORM=(CD*PIA+SD*PIB)/PO                                          ADBP1482
      ELSE                                                              ADBP1483
      PHASE=DATAN2(-PIA,PIB)                                            ADBP1484
      CD=DCOS(PHASE)                                                    ADBP1485
      SD=DSIN(PHASE)                                                    ADBP1486
      RNORM=(CD*PIAP+SD*PIBP)/POP                                       ADBP1487
      ENDIF                                                             ADBP1488
C                                                                       ADBP1489
      IF(RNORM.LT.0.0D0) THEN                                           ADBP1490
      RNORM=-RNORM                                                      ADBP1491
      CD=-CD                                                            ADBP1492
      SD=-SD                                                            ADBP1493
      PHASE=PHASE+PI                                                    ADBP1494
      ENDIF                                                             ADBP1495
      IF(PHASE.GE.0.0D0) THEN                                           ADBP1496
        PHASE=DMOD(PHASE,TPI)                                           ADBP1497
        IF(PHASE.GT.PI) PHASE=PHASE-TPI                                 ADBP1498
      ELSE                                                              ADBP1499
        PHASE=-DMOD(-PHASE,TPI)                                         ADBP1500
        IF(PHASE.LT.-PI) PHASE=PHASE+TPI                                ADBP1501
      ENDIF                                                             ADBP1502
C  ****  NORMALIZED WAVE FUNCTION. (188)                                ADBP1503
      DO 5 I=1,ILAST                                                    ADBP1504
      P(I)=RNORM*P(I)                                                   ADBP1505
      Q(I)=RNORM*Q(I)                                                   ADBP1506
    5 CONTINUE                                                          ADBP1507
      IF(ILAST.LT.NRT) THEN                                             ADBP1508
      DO 6 I=ILAST+1,NRT                                                ADBP1509
      P(I)=CD*PA(I)+SD*PB(I)                                            ADBP1510
      Q(I)=CD*QA(I)+SD*QB(I)                                            ADBP1511
    6 CONTINUE                                                          ADBP1512
      ENDIF                                                             ADBP1513
C                                                                       ADBP1514
C  ****  EXTRACT THE 'RAD' GRID...                                      ADBP1515
C                                                                       ADBP1516
      DO 500 I=1,NGP                                                    ADBP1517
      RLOC=RAD(I)                                                       ADBP1518
      CALL FINDI(R,RLOC,NRT,J)                                          ADBP1519
      IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN                                 ADBP1520
        PIO(I)=P(J)                                                     ADBP1521
        QIO(I)=Q(J)                                                     ADBP1522
      ELSE                                                              ADBP1523
        PIO(I)=P(J+1)                                                   ADBP1524
        QIO(I)=Q(J+1)                                                   ADBP1525
      ENDIF                                                             ADBP1526
  500 CONTINUE                                                          ADBP1527
C                                                                       ADBP1528
      RLOC=R(ILAST)                                                     ADBP1529
      CALL FINDI(RAD,RLOC,NGP,ILAST)                                    ADBP1530
      ILAST=ILAST+1                                                     ADBP1531
      RETURN                                                            ADBP1532
      END                                                               ADBP1533
C  **************************************************************       ADBP1534
C                       SUBROUTINE SOUTW                                ADBP1535
C  **************************************************************       ADBP1536
      SUBROUTINE SOUTW(E,EPS,L,NR,NZERO,IOTP)                           ADBP1537
C                                                                       ADBP1538
C     OUTWARD SOLUTION OF THE SCHRODINGER RADIAL EQUATION FOR A         ADBP1539
C   PIECEWISE CUBIC FIELD. POWER SERIES METHOD.                         ADBP1540
C                                                                       ADBP1541
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1542
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP1543
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1544
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1545
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1546
     1             VD(NPPG),NVT                                         ADBP1547
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1548
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP1549
      COMMON/NZT/NZMAX                                                  ADBP1550
      NZERO=0                                                           ADBP1551
      NZMAX=0                                                           ADBP1552
      AL=L                                                              ADBP1553
      N1=IOTP-1                                                         ADBP1554
C                                                                       ADBP1555
      P(1)=0.0D0                                                        ADBP1556
      Q(1)=0.0D0                                                        ADBP1557
      DO 2 I=1,N1                                                       ADBP1558
      RA=R(I)                                                           ADBP1559
      RB=R(I+1)                                                         ADBP1560
      IN=IND(I)                                                         ADBP1561
      RV0=VA(IN)                                                        ADBP1562
      RV1=VB(IN)                                                        ADBP1563
      RV2=VC(IN)                                                        ADBP1564
      RV3=VD(IN)                                                        ADBP1565
      PI=P(I)                                                           ADBP1566
      QI=Q(I)                                                           ADBP1567
      CALL SCH(E,AL,EPS)                                                ADBP1568
      NZERO=NZERO+NCHS                                                  ADBP1569
      IF(NCHS.GT.NZMAX) NZMAX=NCHS                                      ADBP1570
      IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN                             ADBP1571
      P(I+1)=PF                                                         ADBP1572
      Q(I+1)=QF                                                         ADBP1573
      IF(I.EQ.1) GO TO 2                                                ADBP1574
C  ****  RENORMALIZATION.                                               ADBP1575
      IF(RLN.GT.0.0D0) THEN                                             ADBP1576
      FACT=DEXP(-RLN)                                                   ADBP1577
      DO 1 K=1,I                                                        ADBP1578
      P(K)=P(K)*FACT                                                    ADBP1579
      Q(K)=Q(K)*FACT                                                    ADBP1580
    1 CONTINUE                                                          ADBP1581
      ENDIF                                                             ADBP1582
    2 CONTINUE                                                          ADBP1583
      RETURN                                                            ADBP1584
      END                                                               ADBP1585
C  **************************************************************       ADBP1586
C                       SUBROUTINE SINW                                 ADBP1587
C  **************************************************************       ADBP1588
      SUBROUTINE SINW(E,EPS,L,IOTP)                                     ADBP1589
C                                                                       ADBP1590
C     INWARD SOLUTION OF THE SCHRODINGER RADIAL EQUATION FOR A          ADBP1591
C   PIECEWISE CUBIC FIELD. POWER SERIES METHOD.                         ADBP1592
C                                                                       ADBP1593
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1594
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG,                   ADBP1595
     1  TRINF=5625.0D0)                                                 ADBP1596
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1597
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1598
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1599
     1             VD(NPPG),NVT                                         ADBP1600
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1601
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP1602
      AL=L                                                              ADBP1603
C  ****  WKB SOLUTION AT THE OUTER GRID POINT. (168)                    ADBP1604
      N=NRT                                                             ADBP1605
    1 N1=IND(N-1)                                                       ADBP1606
      RN=R(N)                                                           ADBP1607
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))                      ADBP1608
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))                     ADBP1609
      CMU=2.0D0*RN*(RVN-E*RN)+AL*(AL+1)                                 ADBP1610
      IF(CMU.LE.0.0D0) THEN                                             ADBP1611
      IER=6                                                             ADBP1612
      WRITE(6,1006)                                                     ADBP1613
 1006 FORMAT(1X,'*** ERROR 6 (IN SBOUND): RV(NGP)<<0 OR E>0.',          ADBP1614
     1/5X,'(CHECK THE INPUT POTENTIAL VALUES).')                        ADBP1615
      RETURN                                                            ADBP1616
      ENDIF                                                             ADBP1617
C  ****  PRACTICAL INFINITY. (170)                                      ADBP1618
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN                              ADBP1619
      CRAT=1.0D0-RN*(RVN+RN*RVNP-2*E*RN)/(CMU*CMU)                      ADBP1620
      P(N)=1.0D0                                                        ADBP1621
      Q(N)=(0.5D0/RN)*CRAT-DSQRT(CMU)/RN                                ADBP1622
      ILAST=N                                                           ADBP1623
      ELSE                                                              ADBP1624
      P(N)=0.0D0                                                        ADBP1625
      Q(N)=0.0D0                                                        ADBP1626
      N=N-1                                                             ADBP1627
      GO TO 1                                                           ADBP1628
      ENDIF                                                             ADBP1629
C                                                                       ADBP1630
      N1=N-IOTP                                                         ADBP1631
      DO 3 J=1,N1                                                       ADBP1632
      I=N-J                                                             ADBP1633
      I1=I+1                                                            ADBP1634
      RA=R(I1)                                                          ADBP1635
      RB=R(I)                                                           ADBP1636
      IN=IND(I)                                                         ADBP1637
      RV0=VA(IN)                                                        ADBP1638
      RV1=VB(IN)                                                        ADBP1639
      RV2=VC(IN)                                                        ADBP1640
      RV3=VD(IN)                                                        ADBP1641
      PI=P(I1)                                                          ADBP1642
      QI=Q(I1)                                                          ADBP1643
      CALL SCH(E,AL,EPS)                                                ADBP1644
      P(I)=PF                                                           ADBP1645
      Q(I)=QF                                                           ADBP1646
C  ****  RENORMALIZATION.                                               ADBP1647
      IF(RLN.GT.0.0D0) THEN                                             ADBP1648
      FACT=DEXP(-RLN)                                                   ADBP1649
      DO 2 K=I1,N                                                       ADBP1650
      P(K)=P(K)*FACT                                                    ADBP1651
      Q(K)=Q(K)*FACT                                                    ADBP1652
    2 CONTINUE                                                          ADBP1653
      ENDIF                                                             ADBP1654
    3 CONTINUE                                                          ADBP1655
      RETURN                                                            ADBP1656
      END                                                               ADBP1657
C  **************************************************************       ADBP1658
C                       SUBROUTINE DOUTW                                ADBP1659
C  **************************************************************       ADBP1660
      SUBROUTINE DOUTW(E,EPS,K,NR,NZERO,IOTP)                           ADBP1661
C                                                                       ADBP1662
C     OUTWARD SOLUTION OF THE DIRAC RADIAL EQUATION FOR A PIECE-        ADBP1663
C   WISE CUBIC FIELD. POWER SERIES METHOD.                              ADBP1664
C                                                                       ADBP1665
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1666
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG)                   ADBP1667
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1668
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1669
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1670
     1             VD(NPPG),NVT                                         ADBP1671
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1672
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP1673
      COMMON/NZT/NZMAX                                                  ADBP1674
      NZERO=0                                                           ADBP1675
      NZMAX=0                                                           ADBP1676
      AK=K                                                              ADBP1677
      N1=IOTP-1                                                         ADBP1678
C                                                                       ADBP1679
      P(1)=0.0D0                                                        ADBP1680
      Q(1)=0.0D0                                                        ADBP1681
      DO 2 I=1,N1                                                       ADBP1682
      RA=R(I)                                                           ADBP1683
      RB=R(I+1)                                                         ADBP1684
      IN=IND(I)                                                         ADBP1685
      RV0=VA(IN)                                                        ADBP1686
      RV1=VB(IN)                                                        ADBP1687
      RV2=VC(IN)                                                        ADBP1688
      RV3=VD(IN)                                                        ADBP1689
      PI=P(I)                                                           ADBP1690
      QI=Q(I)                                                           ADBP1691
      CALL DIR(E,AK,EPS)                                                ADBP1692
      NZERO=NZERO+NCHS                                                  ADBP1693
      IF(NCHS.GT.NZMAX) NZMAX=NCHS                                      ADBP1694
      IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN                             ADBP1695
      P(I+1)=PF                                                         ADBP1696
      Q(I+1)=QF                                                         ADBP1697
      IF(I.EQ.1) GO TO 2                                                ADBP1698
C  ****  RENORMALIZATION.                                               ADBP1699
      IF(RLN.GT.0.0D0) THEN                                             ADBP1700
      FACT=DEXP(-RLN)                                                   ADBP1701
      DO 1 J=1,I                                                        ADBP1702
      P(J)=P(J)*FACT                                                    ADBP1703
      Q(J)=Q(J)*FACT                                                    ADBP1704
    1 CONTINUE                                                          ADBP1705
      ENDIF                                                             ADBP1706
    2 CONTINUE                                                          ADBP1707
      RETURN                                                            ADBP1708
      END                                                               ADBP1709
C  **************************************************************       ADBP1710
C                       SUBROUTINE DINW                                 ADBP1711
C  **************************************************************       ADBP1712
      SUBROUTINE DINW(E,EPS,K,IOTP)                                     ADBP1713
C                                                                       ADBP1714
C     INWARD SOLUTION OF THE DIRAC RADIAL EQUATION FOR A PIECE-         ADBP1715
C   WISE CUBIC FIELD. POWER SERIES METHOD.                              ADBP1716
C                                                                       ADBP1717
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1718
      PARAMETER (NDIM=2400,NPPG=NDIM+1,NPTG=NDIM+NPPG,                   ADBP1719
     1  TRINF=5625.0D0,SL=137.036D0)                                    ADBP1720
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER          ADBP1721
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT                ADBP1722
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),        ADBP1723
     1             VD(NPPG),NVT                                         ADBP1724
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1725
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP1726
C  ****  ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER. (15)                  ADBP1727
      IF(K.LT.0) THEN                                                   ADBP1728
      L=-K-1                                                            ADBP1729
      ELSE                                                              ADBP1730
      L=K                                                               ADBP1731
      ENDIF                                                             ADBP1732
      AK=K                                                              ADBP1733
      AL=L                                                              ADBP1734
C  ****  WKB SOLUTION AT THE OUTER GRID POINT. (177)                    ADBP1735
      N=NRT                                                             ADBP1736
      FACT=(E+2.0D0*SL*SL)/(SL*SL)                                      ADBP1737
    1 N1=IND(N-1)                                                       ADBP1738
      RN=R(N)                                                           ADBP1739
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))                      ADBP1740
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))                     ADBP1741
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1)                                  ADBP1742
      IF(CMU.LE.0.0D0) THEN                                             ADBP1743
      IER=6                                                             ADBP1744
      WRITE(6,1006)                                                     ADBP1745
 1006 FORMAT(1X,'*** ERROR 6 (IN DBOUND): RV(NGP)<<0 OR E>0.',          ADBP1746
     1/5X,'(CHECK THE INPUT POTENTIAL VALUES).')                        ADBP1747
      RETURN                                                            ADBP1748
      ENDIF                                                             ADBP1749
C  ****  PRACTICAL INFINITY. (170)                                      ADBP1750
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN                              ADBP1751
      CRAT=(0.5D0-DSQRT(CMU))/RN-0.25D0*FACT*(RVN+RN*RVNP               ADBP1752
     1    -2.0D0*RN*E)/CMU                                              ADBP1753
      P(N)=1.0D0                                                        ADBP1754
      Q(N)=-SL*(CRAT+AK/RN)/(E+2.0D0*SL*SL)                             ADBP1755
      ILAST=N                                                           ADBP1756
      ELSE                                                              ADBP1757
      P(N)=0.0D0                                                        ADBP1758
      Q(N)=0.0D0                                                        ADBP1759
      N=N-1                                                             ADBP1760
      GO TO 1                                                           ADBP1761
      ENDIF                                                             ADBP1762
C                                                                       ADBP1763
      N1=N-IOTP                                                         ADBP1764
      DO 3 J=1,N1                                                       ADBP1765
      I=N-J                                                             ADBP1766
      I1=I+1                                                            ADBP1767
      RA=R(I1)                                                          ADBP1768
      RB=R(I)                                                           ADBP1769
      IN=IND(I)                                                         ADBP1770
      RV0=VA(IN)                                                        ADBP1771
      RV1=VB(IN)                                                        ADBP1772
      RV2=VC(IN)                                                        ADBP1773
      RV3=VD(IN)                                                        ADBP1774
      PI=P(I1)                                                          ADBP1775
      QI=Q(I1)                                                          ADBP1776
      CALL DIR(E,AK,EPS)                                                ADBP1777
      P(I)=PF                                                           ADBP1778
      Q(I)=QF                                                           ADBP1779
C  ****  RENORMALIZATION.                                               ADBP1780
      IF(RLN.GT.0.0D0) THEN                                             ADBP1781
      FACT=DEXP(-RLN)                                                   ADBP1782
      DO 2 M=I1,N                                                       ADBP1783
      P(M)=P(M)*FACT                                                    ADBP1784
      Q(M)=Q(M)*FACT                                                    ADBP1785
    2 CONTINUE                                                          ADBP1786
      ENDIF                                                             ADBP1787
    3 CONTINUE                                                          ADBP1788
      RETURN                                                            ADBP1789
      END                                                               ADBP1790
C  **************************************************************       ADBP1791
C                       SUBROUTINE SCH                                  ADBP1792
C  **************************************************************       ADBP1793
      SUBROUTINE SCH(E,AL,EPS)                                          ADBP1794
C                                                                       ADBP1795
C      THIS SUBROUTINE SOLVES THE RADIAL SCHRODINGER EQUATION FOR       ADBP1796
C   A CENTRAL FIELD V(R) SUCH THAT                                      ADBP1797
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3                     ADBP1798
C      GIVEN THE BOUNDARY CONDITIONS (I.E. THE VALUE OF THE             ADBP1799
C   RADIAL FUNCTION AND ITS DERIVATIVE) AT RA, THE SOLUTION IN          ADBP1800
C   THE INTERVAL BETWEEN RA AND RB IS GENERATED BY USING A              ADBP1801
C   PIECEWISE POWER SERIES EXPANSION FOR A PARTITION OF THE             ADBP1802
C   INTERVAL, SUITABLY CHOSEN TO ALLOW FAST CONVERGENCE OF THE          ADBP1803
C   SERIES.                                                             ADBP1804
C                                                                       ADBP1805
C   INPUT ARGUMENTS:                                                    ADBP1806
C      E ..................... PARTICLE KINETIC ENERGY                  ADBP1807
C      AL .................... ORBITAL ANGULAR MOMENTUM QUANTUM         ADBP1808
C                              NUMBER                                   ADBP1809
C                                                                       ADBP1810
C   INPUT (COMMON POTEN):                                               ADBP1811
C      RV0, RV1, RV2, RV3 .... POTENTIAL PARAMETERS                     ADBP1812
C                                                                       ADBP1813
C   INPUT-OUTPUT (COMMON SINOUT):                                       ADBP1814
C      RA, RB ................ INTERVAL END POINTS (INPUT)              ADBP1815
C      PI, QI ................ VALUES OF THE RADIAL FUNCTION            ADBP1816
C                              AND ITS DERIVATIVE AT RA (INPUT)         ADBP1817
C      PF, QF ................ VALUES OF THE RADIAL FUNCTION            ADBP1818
C                              AND ITS DERIVATIVE AT RB (OUTPUT)        ADBP1819
C      RLN ................... DLOG OF THE RE-NORMALIZING FACTOR        ADBP1820
C      NSTEP ................. NUMBER OF STEPS                          ADBP1821
C      NCHS .................. NUMBER OF ZEROS OF P(R) IN (RA,RB)       ADBP1822
C                                                                       ADBP1823
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1824
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1825
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP1826
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,NSUM                        ADBP1827
      NCHS=0                                                            ADBP1828
      RLN=0.0D0                                                         ADBP1829
C                                                                       ADBP1830
      H=RB-RA                                                           ADBP1831
      IF(H.LT.0.0D0) THEN                                               ADBP1832
      DIRECT=-1.0D0                                                     ADBP1833
      ELSE                                                              ADBP1834
      DIRECT=1.0D0                                                      ADBP1835
      ENDIF                                                             ADBP1836
      K=-2                                                              ADBP1837
      NSTEP=0                                                           ADBP1838
C                                                                       ADBP1839
      R1=RA                                                             ADBP1840
      P1=PI                                                             ADBP1841
      Q1=QI                                                             ADBP1842
    1 R0=R1                                                             ADBP1843
      P0=P1                                                             ADBP1844
      Q0=Q1                                                             ADBP1845
    2 IOUT=0                                                            ADBP1846
      R1=R0+H                                                           ADBP1847
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN                        ADBP1848
      R1=RB                                                             ADBP1849
      H=RB-R0                                                           ADBP1850
      IOUT=1                                                            ADBP1851
      ENDIF                                                             ADBP1852
      CALL SCH0(E,AL,EPS)                                               ADBP1853
C                                                                       ADBP1854
      K=K+1                                                             ADBP1855
      IF(NSUM.GT.15) GO TO 3                                            ADBP1856
      IF(K.LT.0) GO TO 4                                                ADBP1857
      H=H+H                                                             ADBP1858
      K=0                                                               ADBP1859
      GO TO 4                                                           ADBP1860
    3 IF(NSUM.LT.60) GO TO 4                                            ADBP1861
      H=0.5D0*H                                                         ADBP1862
      K=-4                                                              ADBP1863
      GO TO 2                                                           ADBP1864
    4 NSTEP=NSTEP+1                                                     ADBP1865
      TST=DABS(P1)                                                      ADBP1866
      IF(TST.GT.1.0D2) THEN                                             ADBP1867
C  ****  RENORMALIZATION.                                               ADBP1868
      RLN=RLN+DLOG(TST)                                                 ADBP1869
      P1=P1/TST                                                         ADBP1870
      Q1=Q1/TST                                                         ADBP1871
      ENDIF                                                             ADBP1872
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1                    ADBP1873
      IF(IOUT.EQ.0) GO TO 1                                             ADBP1874
C  ****  OUTPUT.                                                        ADBP1875
      PF=P1                                                             ADBP1876
      QF=Q1                                                             ADBP1877
      RETURN                                                            ADBP1878
      END                                                               ADBP1879
C  **************************************************************       ADBP1880
C                       SUBROUTINE SCH0                                 ADBP1881
C  **************************************************************       ADBP1882
      SUBROUTINE SCH0(E,AL,EPS)                                         ADBP1883
C                                                                       ADBP1884
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP1885
C  ****  OVERFLOW LEVEL.                                                ADBP1886
      PARAMETER (OVER=1.0D15)                                           ADBP1887
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP1888
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,NSUM                        ADBP1889
C                                                                       ADBP1890
      RVE=RV1-E                                                         ADBP1891
      IF(R0.GT.1.0D-10) GO TO 3                                         ADBP1892
C                                                                       ADBP1893
C  **** FIRST INTERVAL. (134-143)                                       ADBP1894
C                                                                       ADBP1895
      S=AL+1                                                            ADBP1896
      U0=AL*S                                                           ADBP1897
      U1=2*RV0*R1                                                       ADBP1898
      U2=2*RVE*R1**2                                                    ADBP1899
      U3=2*RV2*R1**3                                                    ADBP1900
      U4=2*RV3*R1**4                                                    ADBP1901
      UT=U0+U1+U2+U3+U4                                                 ADBP1902
C                                                                       ADBP1903
      CA(1)=1.0D0                                                       ADBP1904
      CA(2)=U1*CA(1)/((S+1)*S-U0)                                       ADBP1905
      CA(3)=(U1*CA(2)+U2*CA(1))/((S+2)*(S+1)-U0)                        ADBP1906
      CA(4)=(U1*CA(3)+U2*CA(2)+U3*CA(1))                                ADBP1907
     1     /((S+3)*(S+2)-U0)                                            ADBP1908
C                                                                       ADBP1909
      P1=CA(1)+CA(2)+CA(3)+CA(4)                                        ADBP1910
      Q1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)                    ADBP1911
      P2P1=S*(S-1)*CA(1)+(S+1)*S*CA(2)+(S+2)*(S+1)*CA(3)                ADBP1912
     1     +(S+3)*(S+2)*CA(4)                                           ADBP1913
C                                                                       ADBP1914
      DO 1 I=5,60                                                       ADBP1915
      K=I-1                                                             ADBP1916
      CA(I)=(U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)+U4*CA(I-4))                 ADBP1917
     1     /((S+K)*(S+K-1)-U0)                                          ADBP1918
      P1=P1+CA(I)                                                       ADBP1919
      DQ1=(S+K)*CA(I)                                                   ADBP1920
      Q1=Q1+DQ1                                                         ADBP1921
      P2P1=P2P1+(S+K-1)*DQ1                                             ADBP1922
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP1923
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(P2P1))                           ADBP1924
      IF(TST.GT.OVER) THEN                                              ADBP1925
      NSUM=100                                                          ADBP1926
      RETURN                                                            ADBP1927
      ENDIF                                                             ADBP1928
      T1=DABS(CA(I))                                                    ADBP1929
      T2=DABS(R1*R1*(P2P1-UT*P1))                                       ADBP1930
      TST1=EPS*DMAX1(DABS(P1),DABS(Q1)/I)                               ADBP1931
      TST2=EPS*DMAX1(DABS(P1),DABS(Q1))                                 ADBP1932
      IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 2                             ADBP1933
    1 CONTINUE                                                          ADBP1934
C  ****  RENORMALIZATION. (143)                                         ADBP1935
    2 NSUM=K+1                                                          ADBP1936
      Q1=Q1/(P1*R1)                                                     ADBP1937
      P1=1.0D0                                                          ADBP1938
      RETURN                                                            ADBP1939
C                                                                       ADBP1940
C  **** MIDDLE REGION. (134-139)                                        ADBP1941
C                                                                       ADBP1942
    3 CONTINUE                                                          ADBP1943
      H=R1-R0                                                           ADBP1944
      H2=H*H                                                            ADBP1945
C                                                                       ADBP1946
      RHO=H/R0                                                          ADBP1947
      U0=AL*(AL+1)+2*R0*(RV0+R0*(RVE+R0*(RV2+R0*RV3)))                  ADBP1948
      U1=2*(RV0+R0*(2*RVE+R0*(3*RV2+R0*4*RV3)))*H                       ADBP1949
      U2=2*(RVE+R0*(3*RV2+R0*6*RV3))*H2                                 ADBP1950
      U3=2*(RV2+R0*4*RV3)*H2*H                                          ADBP1951
      U4=2*RV3*H2*H2                                                    ADBP1952
      UT=U0+U1+U2+U3+U4                                                 ADBP1953
C                                                                       ADBP1954
      CA(1)=P0                                                          ADBP1955
      CA(2)=Q0*H                                                        ADBP1956
      CA(3)=RHO*RHO*U0*CA(1)/2                                          ADBP1957
      CA(4)=RHO*(RHO*(U0*CA(2)+U1*CA(1))-4*CA(3))/6                     ADBP1958
      CAK=(U0-2)*CA(3)+U1*CA(2)+U2*CA(1)                                ADBP1959
      CA(5)=RHO*(RHO*CAK-12*CA(4))/12                                   ADBP1960
      CAK=(U0-6)*CA(4)+U1*CA(3)+U2*CA(2)+U3*CA(1)                       ADBP1961
      CA(6)=RHO*(RHO*CAK-24*CA(5))/20                                   ADBP1962
C                                                                       ADBP1963
      P1=CA(1)+CA(2)+CA(3)+CA(4)+CA(5)+CA(6)                            ADBP1964
      Q1=CA(2)+2*CA(3)+3*CA(4)+4*CA(5)+5*CA(6)                          ADBP1965
      P2P1=2*CA(3)+6*CA(4)+12*CA(5)+20*CA(6)                            ADBP1966
C                                                                       ADBP1967
      DO 4 I=7,60                                                       ADBP1968
      K=I-1                                                             ADBP1969
      CAK=(U0-(K-2)*(K-3))*CA(I-2)+U1*CA(I-3)+U2*CA(I-4)                ADBP1970
     1   +U3*CA(I-5)+U4*CA(I-6)                                         ADBP1971
      CA(I)=RHO*(RHO*CAK-2*(K-1)*(K-2)*CA(K))/(K*(K-1))                 ADBP1972
      P1=P1+CA(I)                                                       ADBP1973
      DQ1=K*CA(I)                                                       ADBP1974
      Q1=Q1+DQ1                                                         ADBP1975
      P2P1=P2P1+K*(K-1)*CA(I)                                           ADBP1976
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP1977
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(P2P1))                           ADBP1978
      IF(TST.GT.OVER) THEN                                              ADBP1979
      NSUM=100                                                          ADBP1980
      RETURN                                                            ADBP1981
      ENDIF                                                             ADBP1982
      T1=DABS(CA(I))                                                    ADBP1983
      T2=DABS(R1*R1*P2P1-H2*UT*P1)                                      ADBP1984
      TST1=EPS*DMAX1(DABS(P1),DABS(Q1)/I)                               ADBP1985
      TST2=EPS*DMAX1(DABS(P1),DABS(Q1))                                 ADBP1986
      IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 5                             ADBP1987
    4 CONTINUE                                                          ADBP1988
C                                                                       ADBP1989
    5 NSUM=K+1                                                          ADBP1990
      Q1=Q1/H                                                           ADBP1991
      RETURN                                                            ADBP1992
      END                                                               ADBP1993
C  **************************************************************       ADBP1994
C                       SUBROUTINE DIR                                  ADBP1995
C  **************************************************************       ADBP1996
      SUBROUTINE DIR(E,AK,EPS)                                          ADBP1997
C                                                                       ADBP1998
C      THIS SUBROUTINE SOLVES THE RADIAL DIRAC EQUATION FOR A           ADBP1999
C   CENTRAL FIELD V(R) SUCH THAT                                        ADBP2000
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3                     ADBP2001
C      GIVEN THE BOUNDARY CONDITIONS (I.E. THE VALUE OF THE             ADBP2002
C   LARGE AND SMALL RADIAL FUNCTIONS) AT RA, THE SOLUTION IN            ADBP2003
C   THE INTERVAL BETWEEN RA AND RB IS GENERATED BY USING A              ADBP2004
C   PIECEWISE POWER SERIES EXPANSION FOR A PARTITION OF THE             ADBP2005
C   INTERVAL, SUITABLY CHOSEN TO ALLOW FAST CONVERGENCE OF THE          ADBP2006
C   SERIES.                                                             ADBP2007
C                                                                       ADBP2008
C   INPUT ARGUMENTS:                                                    ADBP2009
C      E ..................... PARTICLE KINETIC ENERGY                  ADBP2010
C      AK .................... RELATIVISTIC ANGULAR MOMENTUM            ADBP2011
C                              QUANTUM NUMBER                           ADBP2012
C                                                                       ADBP2013
C   INPUT (COMMON POTEN):                                               ADBP2014
C      RV0, RV1, RV2, RV3 .... POTENTIAL PARAMETERS                     ADBP2015
C                                                                       ADBP2016
C   INPUT-OUTPUT (COMMON DINOUT):                                       ADBP2017
C      RA, RB ................ INTERVAL END POINTS (INPUT)              ADBP2018
C      PI, QI ................ VALUES OF THE LARGE AND SMALL            ADBP2019
C                              RADIAL FUNCTIONS AT RA (INPUT)           ADBP2020
C      PF, QF ................ VALUES OF THE LARGE AND SMALL            ADBP2021
C                              RADIAL FUNCTIONS AT RB (OUTPUT)          ADBP2022
C      RLN ................... DLOG OF THE RE-NORMALIZING FACTOR        ADBP2023
C      EPS ................... ESTIMATE OF THE GLOBAL ERROR IN          ADBP2024
C                              PF AND QF                                ADBP2025
C      NSTEP ................. NUMBER OF STEPS                          ADBP2026
C      NCHS .................. NUMBER OF ZEROS OF P(R) IN (RA,RB)       ADBP2027
C                                                                       ADBP2028
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP2029
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP2030
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS                    ADBP2031
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM                 ADBP2032
      NCHS=0                                                            ADBP2033
      RLN=0.0D0                                                         ADBP2034
C                                                                       ADBP2035
      H=RB-RA                                                           ADBP2036
      IF(H.LT.0.0D0) THEN                                               ADBP2037
      DIRECT=-1.0D0                                                     ADBP2038
      ELSE                                                              ADBP2039
      DIRECT=1.0D0                                                      ADBP2040
      ENDIF                                                             ADBP2041
      K=-2                                                              ADBP2042
      NSTEP=0                                                           ADBP2043
C                                                                       ADBP2044
      R1=RA                                                             ADBP2045
      P1=PI                                                             ADBP2046
      Q1=QI                                                             ADBP2047
    1 R0=R1                                                             ADBP2048
      P0=P1                                                             ADBP2049
      Q0=Q1                                                             ADBP2050
    2 IOUT=0                                                            ADBP2051
      R1=R0+H                                                           ADBP2052
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN                        ADBP2053
      R1=RB                                                             ADBP2054
      H=RB-R0                                                           ADBP2055
      IOUT=1                                                            ADBP2056
      ENDIF                                                             ADBP2057
      CALL DIR0(E,AK,EPS)                                               ADBP2058
C                                                                       ADBP2059
      K=K+1                                                             ADBP2060
      IF(NSUM.GT.15) GO TO 3                                            ADBP2061
      IF(K.LT.0) GO TO 4                                                ADBP2062
      H=H+H                                                             ADBP2063
      K=0                                                               ADBP2064
      GO TO 4                                                           ADBP2065
    3 IF(NSUM.LT.60) GO TO 4                                            ADBP2066
      H=0.5D0*H                                                         ADBP2067
      K=-4                                                              ADBP2068
      GO TO 2                                                           ADBP2069
    4 NSTEP=NSTEP+1                                                     ADBP2070
      TST=DABS(P1)                                                      ADBP2071
      IF(TST.GT.1.0D2) THEN                                             ADBP2072
C  ****  RENORMALIZATION.                                               ADBP2073
      RLN=RLN+DLOG(TST)                                                 ADBP2074
      P1=P1/TST                                                         ADBP2075
      Q1=Q1/TST                                                         ADBP2076
      ENDIF                                                             ADBP2077
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1                    ADBP2078
      IF(IOUT.EQ.0) GO TO 1                                             ADBP2079
C  ****  OUTPUT.                                                        ADBP2080
      PF=P1                                                             ADBP2081
      QF=Q1                                                             ADBP2082
      RETURN                                                            ADBP2083
      END                                                               ADBP2084
C  **************************************************************       ADBP2085
C                       SUBROUTINE DIR0                                 ADBP2086
C  **************************************************************       ADBP2087
      SUBROUTINE DIR0(E,AK,EPS)                                         ADBP2088
C                                                                       ADBP2089
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP2090
C  ****  SPEED OF LIGHT AND OVERFLOW LEVEL.                             ADBP2091
      PARAMETER (SL=137.036D0,OVER=1.0D15)                              ADBP2092
      COMMON/POTEN/RV0,RV1,RV2,RV3                                      ADBP2093
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM                 ADBP2094
C                                                                       ADBP2095
      ISIG=1                                                            ADBP2096
      IF(AK.GT.0.0D0) ISIG=-1                                           ADBP2097
      H=R1-R0                                                           ADBP2098
      H2=H*H                                                            ADBP2099
      RVE=RV1-E                                                         ADBP2100
C                                                                       ADBP2101
      IF(R0.GT.1.0D-10) GO TO 7                                         ADBP2102
C                                                                       ADBP2103
C  **** FIRST INTERVAL. (147-164)                                       ADBP2104
C                                                                       ADBP2105
      U0=RV0/SL                                                         ADBP2106
      U1=RVE*R1/SL                                                      ADBP2107
      U2=RV2*R1**2/SL                                                   ADBP2108
      U3=RV3*R1**3/SL                                                   ADBP2109
      UT=U0+U1+U2+U3                                                    ADBP2110
      UQ=UT-2*SL*R1                                                     ADBP2111
      UH=U1-2*SL*R1                                                     ADBP2112
      IF(DABS(U0).LT.1.0D-10) GO TO 2                                   ADBP2113
C                                                                       ADBP2114
C  ****  U0.NE.0. (155-159)                                             ADBP2115
      S=DSQRT(AK*AK-U0*U0)                                              ADBP2116
      DS=S+S                                                            ADBP2117
      CA(1)=1.0D0                                                       ADBP2118
      CB(1)=(S+AK)/U0                                                   ADBP2119
      CAI=U1*CA(1)                                                      ADBP2120
      CBI=UH*CB(1)                                                      ADBP2121
      CA(2)=(U0*CAI+(S+1-AK)*CBI)/(DS+1)                                ADBP2122
      CB(2)=(-(S+1+AK)*CAI+U0*CBI)/(DS+1)                               ADBP2123
      CAI=U1*CA(2)+U2*CA(1)                                             ADBP2124
      CBI=UH*CB(2)+U2*CB(1)                                             ADBP2125
      CA(3)=(U0*CAI+(S+2-AK)*CBI)/(2*(DS+2))                            ADBP2126
      CB(3)=(-(S+2+AK)*CAI+U0*CBI)/(2*(DS+2))                           ADBP2127
      P1=CA(1)+CA(2)+CA(3)                                              ADBP2128
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)                               ADBP2129
      Q1=CB(1)+CB(2)+CB(3)                                              ADBP2130
      QP1=S*CB(1)+(S+1)*CB(2)+(S+2)*CB(3)                               ADBP2131
C                                                                       ADBP2132
      DO 1 I=4,60                                                       ADBP2133
      K=I-1                                                             ADBP2134
      CAI=U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)                                ADBP2135
      CBI=UH*CB(K)+U2*CB(I-2)+U3*CB(I-3)                                ADBP2136
      CA(I)=(U0*CAI+(S+K-AK)*CBI)/(K*(DS+K))                            ADBP2137
      CB(I)=(-(S+K+AK)*CAI+U0*CBI)/(K*(DS+K))                           ADBP2138
      P1=P1+CA(I)                                                       ADBP2139
      PP1=PP1+(S+K)*CA(I)                                               ADBP2140
      Q1=Q1+CB(I)                                                       ADBP2141
      QP1=QP1+(S+K)*CB(I)                                               ADBP2142
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP2143
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))                  ADBP2144
      IF(TST.GT.OVER) THEN                                              ADBP2145
      NSUM=100                                                          ADBP2146
      RETURN                                                            ADBP2147
      ENDIF                                                             ADBP2148
      T1A=DABS(R1*PP1+H*(AK*P1-UQ*Q1))                                  ADBP2149
      T1B=DABS(R1*QP1-H*(AK*Q1-UT*P1))                                  ADBP2150
      T1=DMAX1(T1A,T1B)                                                 ADBP2151
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))                                 ADBP2152
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))                                  ADBP2153
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6                               ADBP2154
    1 CONTINUE                                                          ADBP2155
      GO TO 6                                                           ADBP2156
C                                                                       ADBP2157
C  ****  U0.EQ.0 AND SIG=1. (160,161)                                   ADBP2158
    2 IF(ISIG.LT.0) GO TO 4                                             ADBP2159
      S=DABS(AK)                                                        ADBP2160
      DS1=S+S+1                                                         ADBP2161
      CA(1)=1.0D0                                                       ADBP2162
      CB(1)=-U1*CA(1)/DS1                                               ADBP2163
      CA(2)=0.0D0                                                       ADBP2164
      CB(2)=-U2*CA(1)/(DS1+1)                                           ADBP2165
      CA(3)=UH*CB(1)/2                                                  ADBP2166
      CB(3)=-(U1*CA(3)+U3*CA(1))/(DS1+2)                                ADBP2167
      CA(4)=(UH*CB(2)+U2*CB(1))/3                                       ADBP2168
      CB(4)=-(U1*CA(4)+U2*CA(3))/(DS1+3)                                ADBP2169
      P1=CA(1)+CA(2)+CA(3)+CA(4)                                        ADBP2170
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)                   ADBP2171
      Q1=CB(1)+CB(2)+CB(3)+CB(4)                                        ADBP2172
      QP1=(S+1)*CB(1)+(S+2)*CB(2)+(S+3)*CB(3)                           ADBP2173
C                                                                       ADBP2174
      DO 3 I=5,60                                                       ADBP2175
      K=I-1                                                             ADBP2176
      CA(I)=(UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K                        ADBP2177
      CB(I)=-(U1*CA(I)+U2*CA(K)+U3*CA(I-2))/(DS1+K)                     ADBP2178
      P1=P1+CA(I)                                                       ADBP2179
      PP1=PP1+(S+K)*CA(I)                                               ADBP2180
      Q1=Q1+CB(I)                                                       ADBP2181
      QP1=QP1+(S+I)*CB(I)                                               ADBP2182
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP2183
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))                  ADBP2184
      IF(TST.GT.OVER) THEN                                              ADBP2185
      NSUM=100                                                          ADBP2186
      RETURN                                                            ADBP2187
      ENDIF                                                             ADBP2188
      T1A=DABS(R1*PP1+H*(AK*P1-UQ*Q1))                                  ADBP2189
      T1B=DABS(R1*QP1-H*(AK*Q1-UT*P1))                                  ADBP2190
      T1=DMAX1(T1A,T1B)                                                 ADBP2191
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))                                 ADBP2192
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))                                  ADBP2193
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6                               ADBP2194
    3 CONTINUE                                                          ADBP2195
      GO TO 6                                                           ADBP2196
C                                                                       ADBP2197
C  ****  U0.EQ.0 AND SIG=-1. (162,163)                                  ADBP2198
    4 S=DABS(AK)+1                                                      ADBP2199
      DS1=S+DABS(AK)                                                    ADBP2200
      CB(1)=1.0D0                                                       ADBP2201
      CA(1)=UH*CB(1)/DS1                                                ADBP2202
      CB(2)=0.0D0                                                       ADBP2203
      CA(2)=U2*CB(1)/(DS1+1)                                            ADBP2204
      CB(3)=-U1*CA(1)/2                                                 ADBP2205
      CA(3)=(UH*CB(3)+U3*CB(1))/(DS1+2)                                 ADBP2206
      CB(4)=-(U1*CA(2)+U2*CA(1))/3                                      ADBP2207
      CA(4)=(UH*CB(4)+U2*CB(3))/(DS1+3)                                 ADBP2208
      P1=CA(1)+CA(2)+CA(3)+CA(4)                                        ADBP2209
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)                   ADBP2210
      Q1=CB(1)+CB(2)+CB(3)+CB(4)                                        ADBP2211
      QP1=(S-1)*CB(1)+S*CB(2)+(S+1)*CB(3)                               ADBP2212
C                                                                       ADBP2213
      DO 5 I=5,60                                                       ADBP2214
      K=I-1                                                             ADBP2215
      CB(I)=-(U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K                       ADBP2216
      CA(I)=(UH*CB(I)+U2*CB(K)+U3*CB(I-2))/(DS1+K)                      ADBP2217
      P1=P1+CA(I)                                                       ADBP2218
      PP1=PP1+(S+K)*CA(I)                                               ADBP2219
      Q1=Q1+CB(I)                                                       ADBP2220
      QP1=QP1+(S+K-1)*CB(I)                                             ADBP2221
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP2222
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))                  ADBP2223
      IF(TST.GT.OVER) THEN                                              ADBP2224
      NSUM=100                                                          ADBP2225
      RETURN                                                            ADBP2226
      ENDIF                                                             ADBP2227
      T1A=DABS(R1*PP1+H*(AK*P1-UQ*Q1))                                  ADBP2228
      T1B=DABS(R1*QP1-H*(AK*Q1-UT*P1))                                  ADBP2229
      T1=DMAX1(T1A,T1B)                                                 ADBP2230
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))                                 ADBP2231
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))                                  ADBP2232
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6                               ADBP2233
    5 CONTINUE                                                          ADBP2234
C  ****  RENORMALIZATION. (164)                                         ADBP2235
    6 NSUM=K+1                                                          ADBP2236
      Q1=Q1/P1                                                          ADBP2237
      P1=1.0D0                                                          ADBP2238
      RETURN                                                            ADBP2239
C                                                                       ADBP2240
C  **** MIDDLE REGION. (148-152)                                        ADBP2241
C                                                                       ADBP2242
    7 CONTINUE                                                          ADBP2243
      RHO=H/R0                                                          ADBP2244
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/SL                              ADBP2245
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/SL                                 ADBP2246
      U2=(RV2+R0*3*RV3)*H2/SL                                           ADBP2247
      U3=RV3*H*H2/SL                                                    ADBP2248
      UB=U0-2*SL*R0                                                     ADBP2249
      UH=U1-2*SL*H                                                      ADBP2250
      UT=U0+U1+U2+U3                                                    ADBP2251
      UQ=UT-2*SL*R1                                                     ADBP2252
C                                                                       ADBP2253
      CA(1)=P0                                                          ADBP2254
      CB(1)=Q0                                                          ADBP2255
      CA(2)=RHO*(-AK*CA(1)+UB*CB(1))                                    ADBP2256
      CB(2)=-RHO*(U0*CA(1)-AK*CB(1))                                    ADBP2257
      CA(3)=RHO*(-(1+AK)*CA(2)+UB*CB(2)+UH*CB(1))/2                     ADBP2258
      CB(3)=-RHO*(U0*CA(2)+(1-AK)*CB(2)+U1*CA(1))/2                     ADBP2259
      CAI=-(2+AK)*CA(3)+UB*CB(3)+UH*CB(2)+U2*CB(1)                      ADBP2260
      CBI=U0*CA(3)+(2-AK)*CB(3)+U1*CA(2)+U2*CA(1)                       ADBP2261
      CA(4)=RHO*CAI/3                                                   ADBP2262
      CB(4)=-RHO*CBI/3                                                  ADBP2263
C                                                                       ADBP2264
      P1=CA(1)+CA(2)+CA(3)+CA(4)                                        ADBP2265
      PP1=CA(2)+2*CA(3)+3*CA(4)                                         ADBP2266
      Q1=CB(1)+CB(2)+CB(3)+CB(4)                                        ADBP2267
      QP1=CB(2)+2*CB(3)+3*CB(4)                                         ADBP2268
C                                                                       ADBP2269
      DO 9 I=5,60                                                       ADBP2270
      K=I-1                                                             ADBP2271
      CAI=-(K-1+AK)*CA(K)+UB*CB(K)+UH*CB(I-2)                           ADBP2272
     1   +U2*CB(I-3)+U3*CB(I-4)                                         ADBP2273
      CBI=U0*CA(K)+(K-1-AK)*CB(K)+U1*CA(I-2)                            ADBP2274
     1   +U2*CA(I-3)+U3*CA(I-4)                                         ADBP2275
      CA(I)=RHO*CAI/K                                                   ADBP2276
      CB(I)=-RHO*CBI/K                                                  ADBP2277
      P1=P1+CA(I)                                                       ADBP2278
      PP1=PP1+K*CA(I)                                                   ADBP2279
      Q1=Q1+CB(I)                                                       ADBP2280
      QP1=QP1+K*CB(I)                                                   ADBP2281
C  ****  CHECK OVERFLOW LIMIT.                                          ADBP2282
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))                  ADBP2283
      IF(TST.GT.OVER) THEN                                              ADBP2284
      NSUM=100                                                          ADBP2285
      RETURN                                                            ADBP2286
      ENDIF                                                             ADBP2287
      T1A=DABS(R1*PP1+H*(AK*P1-UQ*Q1))                                  ADBP2288
      T1B=DABS(R1*QP1-H*(AK*Q1-UT*P1))                                  ADBP2289
      T1=DMAX1(T1A,T1B)                                                 ADBP2290
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))                                 ADBP2291
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))                                  ADBP2292
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 10                              ADBP2293
    9 CONTINUE                                                          ADBP2294
C                                                                       ADBP2295
   10 NSUM=K+1                                                          ADBP2296
      RETURN                                                            ADBP2297
      END                                                               ADBP2298
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP2299
CCCCCCCCCC         COULOMB AND BESSEL FUNCTIONS        CCCCCCCCCC       ADBP2300
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP2301
C  **************************************************************       ADBP2302
C                       SUBROUTINE SCOUL                                ADBP2303
C  **************************************************************       ADBP2304
      SUBROUTINE SCOUL(Z,E,L,R,F,FP,G,GP,ERR)                           ADBP2305
C                                                                       ADBP2306
C     THIS SUBROUTINE COMPUTES RADIAL SCHRODINGER-COULOMB WAVE          ADBP2307
C  FUNCTIONS FOR FREE STATES.                                           ADBP2308
C                                                                       ADBP2309
C  **** ALL QUANTITIES IN ATOMIC UNITS.                                 ADBP2310
C                                                                       ADBP2311
C  INPUT ARGUMENTS:                                                     ADBP2312
C     Z ........ FIELD STRENGTH, I.E. VALUE OF R*V(R) (ASSUMED          ADBP2313
C                CONSTANT).                                             ADBP2314
C     E ........ PARTICLE KINETIC ENERGY (POSITIVE).                    ADBP2315
C     K ........ ANGULAR MOMENTUM QUANTUM NUMBER KAPPA (.NE.0).         ADBP2316
C     R ........ RADIAL DISTANCE (POSITIVE).                            ADBP2317
C                                                                       ADBP2318
C  OUTPUT ARGUMENTS:                                                    ADBP2319
C     F, FP .... REGULAR RADIAL SCHRODINGER-COULOMB FUNCTION AND        ADBP2320
C                ITS DERIVATIVE.                                        ADBP2321
C     G, GP .... IRREGULAR RADIAL SCHRODINGER-COULOMB FUNCTION          ADBP2322
C                AND ITS DERIVATIVE.                                    ADBP2323
C     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-       ADBP2324
C                CERTAINTY).                                            ADBP2325
C  OUTPUT THROUGH COMMON/OCOUL/:                                        ADBP2326
C     WAVNUM ... WAVE NUMBER.                                           ADBP2327
C     ETA ...... SOMMERFELD'S PARAMETER.                                ADBP2328
C     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).                    ADBP2329
C                                                                       ADBP2330
C     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THEY        ADBP2331
C  OSCILLATE WITH UNIT AMPLITUDE.                                       ADBP2332
C                                                                       ADBP2333
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,         ADBP2334
C                                 AND FUNCTION CLGAM.                   ADBP2335
C                                                                       ADBP2336
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2337
      COMMON/OFCOUL/DELTA0                                              ADBP2338
      COMMON/OCOUL/WAVNUM,ETA,DELTA                                     ADBP2339
C                                                                       ADBP2340

	external besjn

      IF(E.LT.0.0001D0.OR.L.LT.0) THEN                                  ADBP2341
        F=0.0D0                                                         ADBP2342
        FP=0.0D0                                                        ADBP2343
        G=1.0D35                                                        ADBP2344
        GP=-1.0D35                                                      ADBP2345
        ERR=1.0D0                                                       ADBP2346
        IF(E.LT.0.0001D0) WRITE(6,2101)                                 ADBP2347
 2101   FORMAT(1X,'*** ERROR IN SCOUL: E IS TOO SMALL.')                ADBP2348
        IF(L.LT.0) WRITE(6,2102)                                        ADBP2349
 2102   FORMAT(1X,'*** ERROR IN SCOUL: L.LT.0.')                        ADBP2350
        RETURN                                                          ADBP2351
      ENDIF                                                             ADBP2352
C                                                                       ADBP2353
C  ****  PARAMETERS.                                                    ADBP2354
C                                                                       ADBP2355
      WAVNUM=DSQRT(E+E)                                                 ADBP2356
      IF(DABS(Z).GT.0.00001D0) THEN                                     ADBP2357
        ETA=Z/WAVNUM                                                    ADBP2358
        ICAL=0                                                          ADBP2359
      ELSE                                                              ADBP2360
        ETA=0.0D0                                                       ADBP2361
        ICAL=1                                                          ADBP2362
      ENDIF                                                             ADBP2363
      RLAMB=L                                                           ADBP2364
      X=WAVNUM*R                                                        ADBP2365
      IF(ICAL.EQ.1) GO TO 1                                             ADBP2366
C                                                                       ADBP2367
C  ************  COULOMB FUNCTIONS.                                     ADBP2368
C                                                                       ADBP2369
      DELTA0=1.0D30                                                     ADBP2370
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                             ADBP2371
      FP=FP*WAVNUM                                                      ADBP2372
      GP=GP*WAVNUM                                                      ADBP2373
      IF(DELTA0.LT.1.0D29) THEN                                         ADBP2374
        DELTA=DELTA0                                                    ADBP2375
      ELSE                                                              ADBP2376
        DELTA=DELTAC(ETA,RLAMB)                                         ADBP2377
      ENDIF                                                             ADBP2378
      RETURN                                                            ADBP2379
C                                                                       ADBP2380
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                       ADBP2381
C                                                                       ADBP2382
    1 CONTINUE                                                          ADBP2383
      F=X*BESJN(1,L,X)                                                  ADBP2384
      G=-X*BESJN(2,L,X)                                                 ADBP2385
      FP=((L+1)*BESJN(1,L,X)-X*BESJN(1,L+1,X))*WAVNUM                   ADBP2386
      GP=-((L+1)*BESJN(2,L,X)-X*BESJN(2,L+1,X))*WAVNUM                  ADBP2387
      DELTA=0.0D0                                                       ADBP2388
      ERR=0.0D0                                                         ADBP2389
      RETURN                                                            ADBP2390
      END                                                               ADBP2391
C  **************************************************************       ADBP2392
C                       SUBROUTINE DCOUL                                ADBP2393
C  **************************************************************       ADBP2394
      SUBROUTINE DCOUL(Z,E,K,R,FU,FL,GU,GL,ERR)                         ADBP2395
C                                                                       ADBP2396
C     THIS SUBROUTINE COMPUTES RADIAL DIRAC-COULOMB WAVE FUNC-          ADBP2397
C  TIONS FOR FREE STATES.                                               ADBP2398
C                                                                       ADBP2399
C  **** ALL QUANTITIES IN ATOMIC UNITS.                                 ADBP2400
C                                                                       ADBP2401
C  INPUT ARGUMENTS:                                                     ADBP2402
C     Z ........ FIELD STRENGTH, I.E. VALUE OF R*V(R) (ASSUMED          ADBP2403
C                CONSTANT).                                             ADBP2404
C     E ........ PARTICLE KINETIC ENERGY (POSITIVE).                    ADBP2405
C     K ........ ANGULAR MOMENTUM QUANTUM NUMBER KAPPA (.NE.0).         ADBP2406
C     R ........ RADIAL DISTANCE (POSITIVE).                            ADBP2407
C                                                                       ADBP2408
C  OUTPUT ARGUMENTS:                                                    ADBP2409
C     FU, FL ... UPPER AND LOWER COMPONENTS OF THE REGULAR RADIAL       ADBP2410
C                COULOMB FUNCTION.                                      ADBP2411
C     GU, GL ... UPPER AND LOWER COMPONENTS OF THE IRREGULAR RA-        ADBP2412
C                DIAL COULOMB FUNCTION.                                 ADBP2413
C     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-       ADBP2414
C                CERTAINTY).                                            ADBP2415
C  OUTPUT THROUGH COMMON/OCOUL/:                                        ADBP2416
C     WAVNUM ... WAVE NUMBER.                                           ADBP2417
C     ETA ...... SOMMERFELD'S PARAMETER.                                ADBP2418
C     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).                    ADBP2419
C                                                                       ADBP2420
C     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THE         ADBP2421
C  UPPER COMPONENT OSCILLATES WITH UNIT AMPLITUDE.                      ADBP2422
C                                                                       ADBP2423
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,         ADBP2424
C                                 AND FUNCTION CLGAM.                   ADBP2425
C                                                                       ADBP2426
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2427
      PARAMETER (PI=3.1415926535897933D0,PIH=0.5D0*PI,TPI=PI+PI,        ADBP2428
     1  SL=137.036D0,SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)             ADBP2429
      COMMON/OFCOUL/DELTA0                                              ADBP2430
      COMMON/OCOUL/WAVNUM,ETA,DELTA                                     ADBP2431
C                                                                       ADBP2432
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN                                  ADBP2433
        FU=0.0D0                                                        ADBP2434
        FL=0.0D0                                                        ADBP2435
        GU=0.0D0                                                        ADBP2436
        GL=0.0D0                                                        ADBP2437
        ERR=1.0D0                                                       ADBP2438
        IF(E.LT.0.0001D0) WRITE(6,2101)                                 ADBP2439
 2101   FORMAT(1X,'*** ERROR IN DCOUL: E IS TOO SMALL.')                ADBP2440
        IF(K.EQ.0) WRITE(6,2102)                                        ADBP2441
 2102   FORMAT(1X,'*** ERROR IN DCOUL: K.EQ.0.')                        ADBP2442
        RETURN                                                          ADBP2443
      ENDIF                                                             ADBP2444
C                                                                       ADBP2445
C  ****  PARAMETERS.                                                    ADBP2446
C                                                                       ADBP2447
      IF(DABS(Z).GT.0.00001D0) THEN                                     ADBP2448
        ZETA=Z*ALPHA                                                    ADBP2449
        ICAL=0                                                          ADBP2450
      ELSE                                                              ADBP2451
        ZETA=0.0D0                                                      ADBP2452
        ICAL=1                                                          ADBP2453
      ENDIF                                                             ADBP2454
      RLAMBS=K*K-ZETA*ZETA                                              ADBP2455
      RLAMB=DSQRT(RLAMBS)                                               ADBP2456
      W=E+SL2                                                           ADBP2457
      PC=DSQRT(E*(E+TSL2))                                              ADBP2458
      WAVNUM=PC/SL                                                      ADBP2459
      ETA=ZETA*W/PC                                                     ADBP2460
      X=WAVNUM*R                                                        ADBP2461
      RLA=DSQRT(RLAMBS+ETA*ETA)                                         ADBP2462
      IF(ICAL.EQ.1) GO TO 1                                             ADBP2463
C                                                                       ADBP2464
C  ************  COULOMB FUNCTIONS.                                     ADBP2465
C                                                                       ADBP2466
      DELTA0=1.0D30                                                     ADBP2467
      RLAMB1=RLAMB-1.0D0                                                ADBP2468
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR)                    ADBP2469
      IF(ERR.GE.1.0D-6) RETURN                                          ADBP2470
      SLA=(RLAMB/X)+(ETA/RLAMB)                                         ADBP2471
      F=RLAMB*(SLA*FM1-FPM1)/RLA                                        ADBP2472
      G=RLAMB*(SLA*GM1-GPM1)/RLA                                        ADBP2473
C                                                                       ADBP2474
C  ****  DIRAC-COULOMB WAVE FUNCTIONS AND PHASE SHIFT.                  ADBP2475
C                                                                       ADBP2476
      P1=K+RLAMB                                                        ADBP2477
      P2=RLAMB*SL2-K*W                                                  ADBP2478
      RNUR=ZETA*(W+SL2)                                                 ADBP2479
      RNUI=-P1*PC                                                       ADBP2480
      RNORM=1.0D0/(DSQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)                    ADBP2481
      IF(K.GT.0) THEN                                                   ADBP2482
        L=K                                                             ADBP2483
      ELSE                                                              ADBP2484
        L=-K-1                                                          ADBP2485
      ENDIF                                                             ADBP2486
C                                                                       ADBP2487
      IF(DELTA0.GT.1.0D29) DELTA0=DELTAC(ETA,RLAMB)                     ADBP2488
      RNU=DATAN2(RNUI,RNUR)                                             ADBP2489
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0                                  ADBP2490
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN                                    ADBP2491
        RNORM=-RNORM                                                    ADBP2492
        DELTA=DELTA-PI                                                  ADBP2493
      ENDIF                                                             ADBP2494
      IF(DELTA.GE.0.0D0) THEN                                           ADBP2495
        DELTA=DMOD(DELTA,TPI)                                           ADBP2496
      ELSE                                                              ADBP2497
        DELTA=-DMOD(-DELTA,TPI)                                         ADBP2498
      ENDIF                                                             ADBP2499
      Q2=P1*P2*RNORM                                                    ADBP2500
      Q1=RLA*PC*RNORM                                                   ADBP2501
      P1=P1*Q1                                                          ADBP2502
      Q1=ZETA*Q1                                                        ADBP2503
      P2=ZETA*P2*RNORM                                                  ADBP2504
C                                                                       ADBP2505
      FU=P1*F+P2*FM1                                                    ADBP2506
      GU=P1*G+P2*GM1                                                    ADBP2507
      FL=Q1*F+Q2*FM1                                                    ADBP2508
      GL=Q1*G+Q2*GM1                                                    ADBP2509
      RETURN                                                            ADBP2510
C                                                                       ADBP2511
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                       ADBP2512
C                                                                       ADBP2513
    1 CONTINUE                                                          ADBP2514
      RLAMB=IABS(K)                                                     ADBP2515
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)                           ADBP2516
      IF(ERR.GE.1.0D-6) RETURN                                          ADBP2517
      FM1=(RLAMB*F/X)+FP                                                ADBP2518
      GM1=(RLAMB*G/X)+GP                                                ADBP2519
      FACT=DSQRT(E/(E+TSL2))                                            ADBP2520
      IF(K.LT.0) THEN                                                   ADBP2521
        FU=FM1                                                          ADBP2522
        GU=GM1                                                          ADBP2523
        FL=FACT*F                                                       ADBP2524
        GL=FACT*G                                                       ADBP2525
      ELSE                                                              ADBP2526
        FU=F                                                            ADBP2527
        GU=G                                                            ADBP2528
        FL=-FACT*FM1                                                    ADBP2529
        GL=-FACT*GM1                                                    ADBP2530
      ENDIF                                                             ADBP2531
      DELTA=0.0D0                                                       ADBP2532
      RETURN                                                            ADBP2533
      END                                                               ADBP2534
C  **************************************************************       ADBP2535
C                       SUBROUTINE FCOUL                                ADBP2536
C  **************************************************************       ADBP2537
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                       ADBP2538
C                                                                       ADBP2539
C     CALCULATION OF COULOMB FUNCTIONS FOR REAL ETA, RLAMB.GT.-1        ADBP2540
C  AND X LARGER THAN, OR OF THE ORDER OF XTP0 (THE TURNING POINT        ADBP2541
C  FOR RLAMB=0). STEED'S CONTINUED FRACTION METHOD IS COMBINED          ADBP2542
C  WITH RECURSION RELATIONS AND AN ASYMPTOTIC EXPANSION. THE            ADBP2543
C  OUTPUT VALUE ERR=1.0D0 INDICATES THAT THE ADOPTED EVALUATION         ADBP2544
C  ALGORITHM IS NOT APPLICABLE (X IS TOO SMALL).                        ADBP2545
C                                                                       ADBP2546
C  INPUT ARGUMENTS:                                                     ADBP2547
C     ETA ...... SOMMERFELD'S PARAMETER.                                ADBP2548
C     RLAMB .... ANGULAR MOMENTUM.                                      ADBP2549
C     X ........ VARIABLE (=WAVE NUMBER TIMES RADIAL DISTANCE).         ADBP2550
C                                                                       ADBP2551
C  OUTPUT ARGUMENTS:                                                    ADBP2552
C     F, FP .... REGULAR FUNCTION AND ITS DERIVATIVE.                   ADBP2553
C     G, GP .... IRREGULAR FUNCTION AND ITS DERIVATIVE.                 ADBP2554
C     ERR ...... RELATIVE NUMERICAL UNCERTAINTY. A VALUE OF THE         ADBP2555
C                ORDER OF 10**(-N) MEANS THAT THE CALCULATED            ADBP2556
C                FUNCTIONS ARE ACCURATE TO N DECIMAL FIGURES.           ADBP2557
C                THE MAXIMUM ACCURACY ATTAINABLE WITH DOUBLE            ADBP2558
C                PRECISION ARITHMETIC IS ABOUT 1.0D-15.                 ADBP2559
C                                                                       ADBP2560
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINE SUM2F0 AND                 ADBP2561
C                                 FUNCTIONS DELTAC AND CLGAM.           ADBP2562
C                                                                       ADBP2563
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2564
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,        ADBP2565
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)                               ADBP2566
      COMMON/OFCOUL/DELTA                                               ADBP2567
C                                                                       ADBP2568
      IF(RLAMB.LT.-0.999D0) THEN                                        ADBP2569
        WRITE(6,'(1X,''*** ERROR IN RCOUL: RLAMB.LT.-0.999'')')         ADBP2570
        STOP                                                            ADBP2571
      ENDIF                                                             ADBP2572
      IF(X.LT.EPS) GO TO 10                                             ADBP2573
C                                                                       ADBP2574
C  ****  NUMERICAL CONSTANTS.                                           ADBP2575
C                                                                       ADBP2576
      CI=DCMPLX(0.0D0,1.0D0)                                            ADBP2577
      CI2=2.0D0*CI                                                      ADBP2578
      CIETA=CI*ETA                                                      ADBP2579
      X2=X*X                                                            ADBP2580
      ETA2=ETA*ETA                                                      ADBP2581
C                                                                       ADBP2582
C  ****  TURNING POINT (XTP). (44)                                      ADBP2583
C                                                                       ADBP2584
      IF(RLAMB.GE.0.0D0) THEN                                           ADBP2585
        XTP=ETA+DSQRT(ETA2+RLAMB*(RLAMB+1.0D0))                         ADBP2586
      ELSE                                                              ADBP2587
        XTP=EPS                                                         ADBP2588
      ENDIF                                                             ADBP2589
      ERRS=10.0D0                                                       ADBP2590
      IF(X.LT.XTP) GO TO 1                                              ADBP2591
C                                                                       ADBP2592
C  ************  ASYMPTOTIC EXPANSION. (71-75)                          ADBP2593
C                                                                       ADBP2594
C  ****  COULOMB PHASE-SHIFT.                                           ADBP2595
      DELTA=DELTAC(ETA,RLAMB)                                           ADBP2596
C                                                                       ADBP2597
      CPA=CIETA-RLAMB                                                   ADBP2598
      CPB=CIETA+RLAMB+1.0D0                                             ADBP2599
      CPZ=CI2*X                                                         ADBP2600
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)                                ADBP2601
      CQA=CPA+1.0D0                                                     ADBP2602
      CQB=CPB+1.0D0                                                     ADBP2603
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)                               ADBP2604
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)                                 ADBP2605
C  ****  FUNCTIONS.                                                     ADBP2606
      THETA=X-ETA*DLOG(2.0D0*X)-RLAMB*PIH+DELTA                         ADBP2607
      IF(THETA.GT.1.0D4) THETA=DMOD(THETA,TPI)                          ADBP2608
      CEITH=CDEXP(CI*THETA)                                             ADBP2609
      CGIF=C2F0*CEITH                                                   ADBP2610
      G=CGIF                                                            ADBP2611
      F=-CI*CGIF                                                        ADBP2612
C  ****  DERIVATIVES.                                                   ADBP2613
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH                         ADBP2614
      GP=CGIFP                                                          ADBP2615
      FP=-CI*CGIFP                                                      ADBP2616
C  ****  GLOBAL UNCERTAINTY. THE WRONSKIAN MAY DIFFER FROM 1 DUE        ADBP2617
C        TO TRUNCATION AND ROUNDOFF ERRORS.                             ADBP2618
      ERR=DMAX1(ERR1,ERR2,DABS(G*FP-F*GP-1.0D0))                        ADBP2619
      IF(ERR.LE.EPS) RETURN                                             ADBP2620
      ERRS=ERR                                                          ADBP2621
C                                                                       ADBP2622
C  ************  STEED'S CONTINUED FRACTION METHOD.                     ADBP2623
C                                                                       ADBP2624
    1 CONTINUE                                                          ADBP2625
      CIETA2=CIETA+CIETA                                                ADBP2626
      ETAX=ETA*X                                                        ADBP2627
C                                                                       ADBP2628
C  ****  CONTINUED FRACTION FOR F. (60-70)                              ADBP2629
C                                                                       ADBP2630
      INULL=0                                                           ADBP2631
      RLAMBN=RLAMB+1.0D0                                                ADBP2632
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN                      ADBP2633
      B0=(RLAMBN/X)+(ETA/RLAMBN)                                        ADBP2634
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              ADBP2635
      FA3=B0                                                            ADBP2636
      FA2=B0*B1+A1                                                      ADBP2637
      FB3=1.0D0                                                         ADBP2638
      FB2=B1                                                            ADBP2639
      RF=FA3                                                            ADBP2640
C                                                                       ADBP2641
      DO 2 N=2,NTERM                                                    ADBP2642
      RFO=RF                                                            ADBP2643
      DAF=DABS(RF)                                                      ADBP2644
      RLAMBN=RLAMB+N                                                    ADBP2645
      AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2                         ADBP2646
      BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              ADBP2647
      FA1=FA2*BN+FA3*AN                                                 ADBP2648
      FB1=FB2*BN+FB3*AN                                                 ADBP2649
      TST=DABS(FB1)                                                     ADBP2650
C                                                                       ADBP2651
      IF(TST.LT.1.0D-25) THEN                                           ADBP2652
        IF(INULL.GT.0) STOP                                             ADBP2653
        INULL=1                                                         ADBP2654
        FA3=FA2                                                         ADBP2655
        FA2=FA1                                                         ADBP2656
        FB3=FB2                                                         ADBP2657
        FB2=FB1                                                         ADBP2658
        RF=RFO                                                          ADBP2659
      ELSE                                                              ADBP2660
        FA3=FA2/TST                                                     ADBP2661
        FA2=FA1/TST                                                     ADBP2662
        FB3=FB2/TST                                                     ADBP2663
        FB2=FB1/TST                                                     ADBP2664
        RF=FA2/FB2                                                      ADBP2665
        IF(DABS(RF-RFO).LT.EPS*DAF) GO TO 3                             ADBP2666
      ENDIF                                                             ADBP2667
    2 CONTINUE                                                          ADBP2668
    3 CONTINUE                                                          ADBP2669
      IF(DAF.GT.1.0D-25) THEN                                           ADBP2670
        ERRF=DABS(RF-RFO)/DAF                                           ADBP2671
      ELSE                                                              ADBP2672
        ERRF=EPS                                                        ADBP2673
      ENDIF                                                             ADBP2674
      IF(ERRF.GT.ERRS) THEN                                             ADBP2675
        ERR=ERRS                                                        ADBP2676
        RETURN                                                          ADBP2677
      ENDIF                                                             ADBP2678
C                                                                       ADBP2679
C  ****  DOWNWARD RECURSION FOR F AND FP. ONLY IF RLAMB.GT.1 AND        ADBP2680
C        X.LT.XTP. (48,49)                                              ADBP2681
C                                                                       ADBP2682
      RLAMB0=RLAMB                                                      ADBP2683
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN                              ADBP2684
        ISHIFT=0                                                        ADBP2685
      ELSE                                                              ADBP2686
        FT=1.0D0                                                        ADBP2687
        FTP=RF                                                          ADBP2688
        IS0=RLAMB0+1.0D-6                                               ADBP2689
        TST=X*(X-2.0D0*ETA)                                             ADBP2690
        DO 4 I=1,IS0                                                    ADBP2691
        ETARL0=ETA/RLAMB0                                               ADBP2692
        RL=DSQRT(1.0D0+ETARL0**2)                                       ADBP2693
        SL=(RLAMB0/X)+ETARL0                                            ADBP2694
        RLAMB0=RLAMB0-1.0D0                                             ADBP2695
        FTO=FT                                                          ADBP2696
        FT=(SL*FT+FTP)/RL                                               ADBP2697
        FTP=SL*FT-RL*FTO                                                ADBP2698
        IF(FT.GT.1.0D10) THEN                                           ADBP2699
          FTP=FTP/FT                                                    ADBP2700
          FT=1.0D0                                                      ADBP2701
        ENDIF                                                           ADBP2702
        RL1T=RLAMB0*(RLAMB0+1.0D0)                                      ADBP2703
        IF(TST.GT.RL1T) THEN                                            ADBP2704
          ISHIFT=I                                                      ADBP2705
          GO TO 5                                                       ADBP2706
        ENDIF                                                           ADBP2707
    4   CONTINUE                                                        ADBP2708
        ISHIFT=IS0                                                      ADBP2709
    5   CONTINUE                                                        ADBP2710
        XTPC=ETA+DSQRT(ETA2+RL1T)                                       ADBP2711
        RFM=FTP/FT                                                      ADBP2712
      ENDIF                                                             ADBP2713
C                                                                       ADBP2714
C  ****  CONTINUED FRACTION FOR P+CI*Q WITH RLAMB0. (76-79)             ADBP2715
C                                                                       ADBP2716
      INULL=0                                                           ADBP2717
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)                              ADBP2718
      CB0=X-ETA                                                         ADBP2719
      CBN=2.0D0*(X-ETA+CI)                                              ADBP2720
      CFA3=CB0                                                          ADBP2721
      CFA2=CB0*CBN+CAN                                                  ADBP2722
      CFB3=1.0D0                                                        ADBP2723
      CFB2=CBN                                                          ADBP2724
      CPIQ=CFA3                                                         ADBP2725
C                                                                       ADBP2726
      DO 6 N=2,NTERM                                                    ADBP2727
      CPIQO=CPIQ                                                        ADBP2728
      DAPIQ=CDABS(CPIQ)                                                 ADBP2729
      CAN=CAN+CIETA2+(N+N-2)                                            ADBP2730
      CBN=CBN+CI2                                                       ADBP2731
      CFA1=CFA2*CBN+CFA3*CAN                                            ADBP2732
      CFB1=CFB2*CBN+CFB3*CAN                                            ADBP2733
      TST=CDABS(CFB1)                                                   ADBP2734
C                                                                       ADBP2735
      IF(TST.LT.1.0D-25) THEN                                           ADBP2736
        IF(INULL.GT.0) STOP                                             ADBP2737
        INULL=1                                                         ADBP2738
        CFA3=CFA2                                                       ADBP2739
        CFA2=CFA1                                                       ADBP2740
        CFB3=CFB2                                                       ADBP2741
        CFB2=CFB1                                                       ADBP2742
        CPIQ=CPIQO                                                      ADBP2743
      ELSE                                                              ADBP2744
        CFA3=CFA2/TST                                                   ADBP2745
        CFA2=CFA1/TST                                                   ADBP2746
        CFB3=CFB2/TST                                                   ADBP2747
        CFB2=CFB1/TST                                                   ADBP2748
        CPIQ=CFA2/CFB2                                                  ADBP2749
        IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 7                      ADBP2750
      ENDIF                                                             ADBP2751
    6 CONTINUE                                                          ADBP2752
    7 CONTINUE                                                          ADBP2753
      IF(DAPIQ.GT.1.0D-25) THEN                                         ADBP2754
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ                                  ADBP2755
      ELSE                                                              ADBP2756
        ERRPIQ=EPS                                                      ADBP2757
      ENDIF                                                             ADBP2758
      IF(ERRPIQ.GT.ERRS) THEN                                           ADBP2759
        ERR=ERRS                                                        ADBP2760
        RETURN                                                          ADBP2761
      ENDIF                                                             ADBP2762
      CPIQ=CI*CPIQ/X                                                    ADBP2763
C                                                                       ADBP2764
      RP=CPIQ                                                           ADBP2765
      RQ=-CI*CPIQ                                                       ADBP2766
      IF(RQ.LE.1.0D-25) GO TO 10                                        ADBP2767
      ERR=DMAX1(ERRF,ERRPIQ)                                            ADBP2768
C                                                                       ADBP2769
C  ****  INVERTING STEED'S TRANSFORMATION. (57,58)                      ADBP2770
C                                                                       ADBP2771
      IF(ISHIFT.LT.1) THEN                                              ADBP2772
        RFP=RF-RP                                                       ADBP2773
        F=DSQRT(RQ/(RFP**2+RQ**2))                                      ADBP2774
        IF(FB2.LT.0.0D0) F=-F                                           ADBP2775
        FP=RF*F                                                         ADBP2776
        G=RFP*F/RQ                                                      ADBP2777
        GP=(RP*RFP-RQ**2)*F/RQ                                          ADBP2778
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 10                            ADBP2779
      ELSE                                                              ADBP2780
        RFP=RFM-RP                                                      ADBP2781
        FM=DSQRT(RQ/(RFP**2+RQ**2))                                     ADBP2782
        G=RFP*FM/RQ                                                     ADBP2783
        GP=(RP*RFP-RQ**2)*FM/RQ                                         ADBP2784
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 10                          ADBP2785
C  ****  UPWARD RECURSION FOR G AND GP (IF ISHIFT.GT.0). (50,51)        ADBP2786
        DO 8 I=1,ISHIFT                                                 ADBP2787
        RLAMB0=RLAMB0+1.0D0                                             ADBP2788
        ETARL0=ETA/RLAMB0                                               ADBP2789
        RL=DSQRT(1.0D0+ETARL0**2)                                       ADBP2790
        SL=(RLAMB0/X)+ETARL0                                            ADBP2791
        GO=G                                                            ADBP2792
        G=(SL*GO-GP)/RL                                                 ADBP2793
        GP=RL*GO-SL*G                                                   ADBP2794
        IF(G.GT.1.0D35) GO TO 10                                        ADBP2795
    8   CONTINUE                                                        ADBP2796
    9   W=RF*G-GP                                                       ADBP2797
        F=1.0D0/W                                                       ADBP2798
        FP=RF/W                                                         ADBP2799
      ENDIF                                                             ADBP2800
C  ****  THE WRONSKIAN MAY DIFFER FROM 1 DUE TO ROUNDOFF ERRORS.        ADBP2801
      ERR=DMAX1(ERR,DABS(FP*G-F*GP-1.0D0))                              ADBP2802
      RETURN                                                            ADBP2803
C                                                                       ADBP2804
   10 F=0.0D0                                                           ADBP2805
      FP=0.0D0                                                          ADBP2806
      G=1.0D35                                                          ADBP2807
      GP=-1.0D35                                                        ADBP2808
      ERR=1.0D0                                                         ADBP2809
      RETURN                                                            ADBP2810
      END                                                               ADBP2811
C  **************************************************************       ADBP2812
C                       SUBROUTINE SUM2F0                               ADBP2813
C  **************************************************************       ADBP2814
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)                                ADBP2815
C                                                                       ADBP2816
C     SUMMATION OF THE 2F0(CA,CB;CS) HYPERGEOMETRIC ASYMPTOTIC          ADBP2817
C  SERIES. THE POSITIVE AND NEGATIVE CONTRIBUTIONS TO THE REAL          ADBP2818
C  AND IMAGINARY PARTS ARE ADDED SEPARATELY TO OBTAIN AN ESTIMATE       ADBP2819
C  OF ROUNDING ERRORS.                                                  ADBP2820
C                                                                       ADBP2821
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2822
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)                    ADBP2823
      RRP=1.0D0                                                         ADBP2824
      RRN=0.0D0                                                         ADBP2825
      RIP=0.0D0                                                         ADBP2826
      RIN=0.0D0                                                         ADBP2827
      CDF=1.0D0                                                         ADBP2828
      ERR2=0.0D0                                                        ADBP2829
      ERR3=1.0D0                                                        ADBP2830
      DO 1 I=1,NTERM                                                    ADBP2831
      J=I-1                                                             ADBP2832
      CDF=CDF*(CA+J)*(CB+J)/(I*CZ)                                      ADBP2833
      ERR1=ERR2                                                         ADBP2834
      ERR2=ERR3                                                         ADBP2835
      ERR3=CDABS(CDF)                                                   ADBP2836
      IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 2                         ADBP2837
      AR=CDF                                                            ADBP2838
      IF(AR.GT.0.0D0) THEN                                              ADBP2839
        RRP=RRP+AR                                                      ADBP2840
      ELSE                                                              ADBP2841
        RRN=RRN+AR                                                      ADBP2842
      ENDIF                                                             ADBP2843
      AI=DCMPLX(0.0D0,-1.0D0)*CDF                                       ADBP2844
      IF(AI.GT.0.0D0) THEN                                              ADBP2845
        RIP=RIP+AI                                                      ADBP2846
      ELSE                                                              ADBP2847
        RIN=RIN+AI                                                      ADBP2848
      ENDIF                                                             ADBP2849
      CF=DCMPLX(RRP+RRN,RIP+RIN)                                        ADBP2850
      AF=CDABS(CF)                                                      ADBP2851
      IF(AF.GT.1.0D25) THEN                                             ADBP2852
        CF=0.0D0                                                        ADBP2853
        ERR=1.0D0                                                       ADBP2854
        RETURN                                                          ADBP2855
      ENDIF                                                             ADBP2856
      IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN                        ADBP2857
         ERR=EPS                                                        ADBP2858
         RETURN                                                         ADBP2859
      ENDIF                                                             ADBP2860
    1 CONTINUE                                                          ADBP2861
C  ****  ROUNDOFF ERROR.                                                ADBP2862
    2 CONTINUE                                                          ADBP2863
      TR=DABS(RRP+RRN)                                                  ADBP2864
      IF(TR.GT.1.0D-25) THEN                                            ADBP2865
        ERRR=(RRP-RRN)*ACCUR/TR                                         ADBP2866
      ELSE                                                              ADBP2867
        ERRR=1.0D0                                                      ADBP2868
      ENDIF                                                             ADBP2869
      TI=DABS(RIP+RIN)                                                  ADBP2870
      IF(TI.GT.1.0D-25) THEN                                            ADBP2871
        ERRI=(RIP-RIN)*ACCUR/TI                                         ADBP2872
      ELSE                                                              ADBP2873
        ERRI=1.0D0                                                      ADBP2874
      ENDIF                                                             ADBP2875
C  ****  ... AND TRUNCATION ERROR.                                      ADBP2876
      IF(AR.GT.1.0D-25) THEN                                            ADBP2877
      ERR=DMAX1(ERRR,ERRI)+ERR2/AF                                      ADBP2878
      ELSE                                                              ADBP2879
      ERR=DMAX1(ERRR,ERRI)                                              ADBP2880
      ENDIF                                                             ADBP2881
      RETURN                                                            ADBP2882
      END                                                               ADBP2883
C  **************************************************************       ADBP2884
C                         FUNCTION DELTAC                               ADBP2885
C  **************************************************************       ADBP2886
      FUNCTION DELTAC(ETA,RLAMB)                                        ADBP2887
C                                                                       ADBP2888
C     CALCULATION OF COULOMB PHASE SHIFT (MODULUS 2*PI). (47)           ADBP2889
C                                                                       ADBP2890
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2891
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                     ADBP2892
      CI=DCMPLX(0.0D0,1.0D0)                                            ADBP2893
C  ****  COULOMB PHASE-SHIFT.                                           ADBP2894
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)                              ADBP2895
      IF(DELTAC.GE.0.0D0) THEN                                          ADBP2896
        DELTAC=DMOD(DELTAC,TPI)                                         ADBP2897
      ELSE                                                              ADBP2898
        DELTAC=-DMOD(-DELTAC,TPI)                                       ADBP2899
      ENDIF                                                             ADBP2900
      RETURN                                                            ADBP2901
      END                                                               ADBP2902
C  **************************************************************       ADBP2903
C                       FUNCTION CLGAM                                  ADBP2904
C  **************************************************************       ADBP2905
      FUNCTION CLGAM(CZ)                                                ADBP2906
C                                                                       ADBP2907
C     THIS FUNCTION GIVES LOG(GAMMA(CZ)) FOR COMPLEX ARGUMENTS.         ADBP2908
C                                                                       ADBP2909
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-        ADBP2910
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 255-257.           ADBP2911
C                                                                       ADBP2912
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2913
      PARAMETER (PI=3.1415926535897932D0)                               ADBP2914
      CZA=CZ                                                            ADBP2915
      ICONJ=0                                                           ADBP2916
      AR=CZA                                                            ADBP2917
      CLGAM=36.84136149D0                                               ADBP2918
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  ADBP2919
C                                                                       ADBP2920
      AI=CZA*DCMPLX(0.0D0,-1.0D0)                                       ADBP2921
      IF(AI.GT.0.0D0) THEN                                              ADBP2922
        ICONJ=0                                                         ADBP2923
      ELSE                                                              ADBP2924
        ICONJ=1                                                         ADBP2925
        CZA=DCONJG(CZA)                                                 ADBP2926
      ENDIF                                                             ADBP2927
C                                                                       ADBP2928
      CZFAC=1.0D0                                                       ADBP2929
      CZFL=0.0D0                                                        ADBP2930
    1 CZFAC=CZFAC/CZA                                                   ADBP2931
      IF(CDABS(CZFAC).GT.1.0D8) THEN                                    ADBP2932
        CZFL=CZFL+CDLOG(CZFAC)                                          ADBP2933
        CZFAC=1.0D0                                                     ADBP2934
      ENDIF                                                             ADBP2935
      CZA=CZA+1.0D0                                                     ADBP2936
      AR=CZA                                                            ADBP2937
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  ADBP2938
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2                  ADBP2939
      GO TO 1                                                           ADBP2940
C  ****  STIRLING'S EXPANSION OF CDLOG(GAMMA(CZA)).                     ADBP2941
    2 CZI2=1.0D0/(CZA*CZA)                                              ADBP2942
      CZS=(43867.0D0/244188.0D0)*CZI2                                   ADBP2943
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2                                ADBP2944
      CZS=(CZS+1.0D0/156.0D0)*CZI2                                      ADBP2945
      CZS=(CZS-691.0D0/360360.0D0)*CZI2                                 ADBP2946
      CZS=(CZS+1.0D0/1188.0D0)*CZI2                                     ADBP2947
      CZS=(CZS-1.0D0/1680.0D0)*CZI2                                     ADBP2948
      CZS=(CZS+1.0D0/1260.0D0)*CZI2                                     ADBP2949
      CZS=(CZS-1.0D0/360.0D0)*CZI2                                      ADBP2950
      CZS=(CZS+1.0D0/12.0D0)/CZA                                        ADBP2951
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS        ADBP2952
     1     +CZFL+CDLOG(CZFAC)                                           ADBP2953
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)                                ADBP2954
      RETURN                                                            ADBP2955
      END                                                               ADBP2956
C  **************************************************************       ADBP2957
C                         FUNCION BESJN                                 ADBP2958
C  **************************************************************       ADBP2959
      FUNCTION BESJN(JY,N,X)                                            ADBP2960
C                                                                       ADBP2961
C      THIS FUNCTION COMPUTES THE SPHERICAL BESSEL FUNCTIONS OF         ADBP2962
C   THE FIRST KIND AND SPHERICAL BESSEL FUNCTIONS OF THE SECOND         ADBP2963
C   KIND (ALSO KNOWN AS SPHERICAL NEUMANN FUNCTIONS) FOR REAL           ADBP2964
C   POSITIVE ARGUMENTS.                                                 ADBP2965
C                                                                       ADBP2966
C      INPUT:                                                           ADBP2967
C         JY ...... KIND: 1(BESSEL) OR 2(NEUMANN).                      ADBP2968
C         N ....... ORDER (INTEGER).                                    ADBP2969
C         X ....... ARGUMENT (REAL AND POSITIVE).                       ADBP2970
C                                                                       ADBP2971
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-        ADBP2972
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 435-478.           ADBP2973
C                                                                       ADBP2974
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP2975
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                     ADBP2976
      IF(X.LT.0) THEN                                                   ADBP2977
        WRITE(6,1000)                                                   ADBP2978
 1000   FORMAT(1X,'*** NEGATIVE ARGUMENT IN FUNCTION BESJN.')           ADBP2979
        STOP                                                            ADBP2980
      ENDIF                                                             ADBP2981
C  ****  ORDER AND PHASE CORRECTION FOR NEUMANN FUNCTIONS.              ADBP2982
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.15.                            ADBP2983
      IF(JY.EQ.2) THEN                                                  ADBP2984
        NL=-N-1                                                         ADBP2985
        IPH=2*MOD(IABS(N),2)-1                                          ADBP2986
      ELSE                                                              ADBP2987
        NL=N                                                            ADBP2988
        IPH=1                                                           ADBP2989
      ENDIF                                                             ADBP2990
C  ****  SELECTION OF CALCULATION MODE.                                 ADBP2991
      IF(NL.LT.0) GO TO 10                                              ADBP2992
      IF(X.GT.1.0D0*NL) GO TO 7                                         ADBP2993
      XI=X*X                                                            ADBP2994
      IF(XI.GT.NL+NL+3.0D0) GO TO 4                                     ADBP2995
C  ****  POWER SERIES FOR SMALL ARGUMENTS AND POSITIVE ORDERS.          ADBP2996
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.2.                             ADBP2997
      F1=1.0D0                                                          ADBP2998
      IP=1                                                              ADBP2999
      IF(NL.NE.0) THEN                                                  ADBP3000
        DO 1 I=1,NL                                                     ADBP3001
        IP=IP+2                                                         ADBP3002
    1   F1=F1*X/IP                                                      ADBP3003
      ENDIF                                                             ADBP3004
      XI=0.5D0*XI                                                       ADBP3005
      BESJN=1.0D0                                                       ADBP3006
      PS=1.0D0                                                          ADBP3007
      DO 2 I=1,500                                                      ADBP3008
      IP=IP+2                                                           ADBP3009
      PS=-PS*XI/(I*IP)                                                  ADBP3010
      BESJN=BESJN+PS                                                    ADBP3011
      IF(DABS(PS).LT.1.0D-18*DABS(BESJN)) GO TO 3                       ADBP3012
    2 CONTINUE                                                          ADBP3013
    3 BESJN=IPH*F1*BESJN                                                ADBP3014
      RETURN                                                            ADBP3015
C  ****  MILLER'S METHOD FOR POSITIVE ORDERS AND INTERMEDIATE           ADBP3016
C        ARGUMENTS. ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                 ADBP3017
    4 XI=1.0D0/X                                                        ADBP3018
      F2=0.0D0                                                          ADBP3019
      F3=1.0D-35                                                        ADBP3020
      IP=2*(NL+31)+3                                                    ADBP3021
      DO 5 I=1,31                                                       ADBP3022
      F1=F2                                                             ADBP3023
      F2=F3                                                             ADBP3024
      IP=IP-2                                                           ADBP3025
      F3=IP*XI*F2-F1                                                    ADBP3026
      IF(DABS(F3).GT.1.0D30) THEN                                       ADBP3027
        F2=F2/F3                                                        ADBP3028
        F3=1.0D0                                                        ADBP3029
      ENDIF                                                             ADBP3030
    5 CONTINUE                                                          ADBP3031
      BESJN=1.0D0                                                       ADBP3032
      F2=F2/F3                                                          ADBP3033
      F3=1.0D0                                                          ADBP3034
      DO 6 I=1,NL                                                       ADBP3035
      F1=F2                                                             ADBP3036
      F2=F3                                                             ADBP3037
      IP=IP-2                                                           ADBP3038
      F3=IP*XI*F2-F1                                                    ADBP3039
      IF(DABS(F3).GT.1.0D30) THEN                                       ADBP3040
        BESJN=BESJN/F3                                                  ADBP3041
        F2=F2/F3                                                        ADBP3042
        F3=1.0D0                                                        ADBP3043
      ENDIF                                                             ADBP3044
    6 CONTINUE                                                          ADBP3045
      BESJN=IPH*XI*DSIN(X)*BESJN/F3                                     ADBP3046
      RETURN                                                            ADBP3047
C  ****  RECURRENCE RELATION FOR ARGUMENTS GREATER THAN ORDER.          ADBP3048
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                            ADBP3049
    7 XI=1.0D0/X                                                        ADBP3050
      F3=XI*DSIN(X)                                                     ADBP3051
      IF(NL.EQ.0) GO TO 9                                               ADBP3052
      F2=F3                                                             ADBP3053
      F3=XI*(F2-DCOS(X))                                                ADBP3054
      IF(NL.EQ.1) GO TO 9                                               ADBP3055
      IP=1                                                              ADBP3056
      DO 8 I=2,NL                                                       ADBP3057
      F1=F2                                                             ADBP3058
      F2=F3                                                             ADBP3059
      IP=IP+2                                                           ADBP3060
    8 F3=IP*XI*F2-F1                                                    ADBP3061
    9 BESJN=IPH*F3                                                      ADBP3062
      RETURN                                                            ADBP3063
C  ****  RECURRENCE RELATION FOR NEGATIVE ORDERS.                       ADBP3064
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                            ADBP3065
   10 NL=IABS(NL)                                                       ADBP3066
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN              ADBP3067
        BESJN=-1.0D35                                                   ADBP3068
        RETURN                                                          ADBP3069
      ENDIF                                                             ADBP3070
      XI=1.0D0/X                                                        ADBP3071
      F3=XI*DSIN(X)                                                     ADBP3072
      F2=XI*(F3-DCOS(X))                                                ADBP3073
      IP=3                                                              ADBP3074
      DO 11 I=1,NL                                                      ADBP3075
      F1=F2                                                             ADBP3076
      F2=F3                                                             ADBP3077
      IP=IP-2                                                           ADBP3078
      F3=IP*XI*F2-F1                                                    ADBP3079
      IF(DABS(F3).GT.1.0D35) THEN                                       ADBP3080
        BESJN=-1.0D35                                                   ADBP3081
        RETURN                                                          ADBP3082
      ENDIF                                                             ADBP3083
   11 CONTINUE                                                          ADBP3084
      BESJN=IPH*F3                                                      ADBP3085
      RETURN                                                            ADBP3086
      END                                                               ADBP3087
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP3088
CCCCCCCCCC          CUBIC SPLINE INTERPOLATION         CCCCCCCCCC       ADBP3089
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP3090
C  **************************************************************       ADBP3091
C                       SUBROUTINE SPLINE                               ADBP3092
C  **************************************************************       ADBP3093
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)                            ADBP3094
C                                                                       ADBP3095
C      CUBIC SPLINE INTERPOLATION BETWEEN TABULATED DATA.               ADBP3096
C   INPUT:                                                              ADBP3097
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           ADBP3098
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        ADBP3099
C     Y(I) (I=1, ...,N) ........ CORRESPONDING FUNCTION VALUES.         ADBP3100
C     S1,SN ..... SECOND DERIVATIVES AT X(1) AND X(N).                  ADBP3101
C            (THE NATURAL SPLINE CORRESPONDS TO TAKING S1=SN=0).        ADBP3102
C     N ........................ NUMBER OF GRID POINTS.                 ADBP3103
C      THE INTERPOLATING POLYNOMIAL IN THE I-TH INTERVAL, FROM          ADBP3104
C   X(I) TO X(I+1), IS                                                  ADBP3105
C            PI(X) = A(I)+X*(B(I)+X*(C(I)+X*D(I)))                      ADBP3106
C   OUTPUT:                                                             ADBP3107
C     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.                   ADBP3108
C                                                                       ADBP3109
C      REF.: M.J. MARON, 'NUMERICAL ANALYSIS: A PRACTICAL               ADBP3110
C            APPROACH', MACMILLAN PUBL. CO., NEW YORK 1982.             ADBP3111
C                                                                       ADBP3112
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP3113
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)                           ADBP3114
      IF(N.LT.4) THEN                                                   ADBP3115
      WRITE(6,10) N                                                     ADBP3116
   10 FORMAT(5X,'SPLINE INTERPOLATION CANNOT BE PERFORMED WITH',        ADBP3117
     1I4,' POINTS. STOP.')                                              ADBP3118
      STOP                                                              ADBP3119
      ENDIF                                                             ADBP3120
      N1=N-1                                                            ADBP3121
      N2=N-2                                                            ADBP3122
C  ****  AUXILIARY ARRAYS H(=A) AND DELTA(=D).                          ADBP3123
      DO 1 I=1,N1                                                       ADBP3124

      IF(X(I+1)-X(I).LT.1.0D-10) THEN                                   ADBP3125
      WRITE(6,11)                                                       ADBP3126
   11 FORMAT(5X,'SPLINE X VALUES NOT IN INCREASING ORDER. STOP.')       ADBP3127
      STOP                                                              ADBP3128
      ENDIF                                                             ADBP3129
      A(I)=X(I+1)-X(I)                                                  ADBP3130
    1 D(I)=(Y(I+1)-Y(I))/A(I)                                           ADBP3131
C  ****  SYMMETRIC COEFFICIENT MATRIX (AUGMENTED).                      ADBP3132
      DO 2 I=1,N2                                                       ADBP3133
      B(I)=2.0D0*(A(I)+A(I+1))                                          ADBP3134
      K=N1-I+1                                                          ADBP3135
    2 D(K)=6.0D0*(D(K)-D(K-1))                                          ADBP3136
      D(2)=D(2)-A(1)*S1                                                 ADBP3137
      D(N1)=D(N1)-A(N1)*SN                                              ADBP3138
C  ****  GAUSS SOLUTION OF THE TRIDIAGONAL SYSTEM.                      ADBP3139
      DO 3 I=2,N2                                                       ADBP3140
      R=A(I)/B(I-1)                                                     ADBP3141
      B(I)=B(I)-R*A(I)                                                  ADBP3142
    3 D(I+1)=D(I+1)-R*D(I)                                              ADBP3143
C  ****  THE SIGMA COEFFICIENTS ARE STORED IN ARRAY D.                  ADBP3144
      D(N1)=D(N1)/B(N2)                                                 ADBP3145
      DO 4 I=2,N2                                                       ADBP3146
      K=N1-I+1                                                          ADBP3147
    4 D(K)=(D(K)-A(K)*D(K+1))/B(K-1)                                    ADBP3148
      D(N)=SN                                                           ADBP3149
C  ****  SPLINE COEFFICIENTS.                                           ADBP3150
      SI1=S1                                                            ADBP3151
      DO 5 I=1,N1                                                       ADBP3152
      SI=SI1                                                            ADBP3153
      SI1=D(I+1)                                                        ADBP3154
      H=A(I)                                                            ADBP3155
      HI=1.0D0/H                                                        ADBP3156
      A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)                        ADBP3157
     1    +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))                                 ADBP3158
     2    +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))                               ADBP3159
      B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)                        ADBP3160
     1    +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)                          ADBP3161
      C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))                              ADBP3162
    5 D(I)=(HI/6.0D0)*(SI1-SI)                                          ADBP3163
      RETURN                                                            ADBP3164
      END                                                               ADBP3165
C  **************************************************************       ADBP3166
C                       SUBROUTINE FINDI                                ADBP3167
C  **************************************************************       ADBP3168
      SUBROUTINE FINDI(X,XC,N,I)                                        ADBP3169
C                                                                       ADBP3170
C      FINDS THE INTERVAL (X(I),X(I+1)) CONTAINING THE VALUE XC.        ADBP3171
C   INPUT:                                                              ADBP3172
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           ADBP3173
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        ADBP3174
C     XC ....................... POINT TO BE LOCATED.                   ADBP3175
C     N ........................ NUMBER OF GRID POINTS.                 ADBP3176
C   OUTPUT:                                                             ADBP3177
C     I ........................ INTERVAL INDEX.                        ADBP3178
C                                                                       ADBP3179
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP3180
      DIMENSION X(N)                                                    ADBP3181
      IF(XC.GT.X(N)) THEN                                               ADBP3182
      I=N-1                                                             ADBP3183
      RETURN                                                            ADBP3184
      ENDIF                                                             ADBP3185
      IF(XC.LT.X(1)) THEN                                               ADBP3186
      I=1                                                               ADBP3187
      RETURN                                                            ADBP3188
      ENDIF                                                             ADBP3189
      I=1                                                               ADBP3190
      I1=N                                                              ADBP3191
    1 IT=(I+I1)/2                                                       ADBP3192
      IF(XC.GT.X(IT)) I=IT                                              ADBP3193
      IF(XC.LE.X(IT)) I1=IT                                             ADBP3194
      IF(I1-I.GT.1) GO TO 1                                             ADBP3195
      RETURN                                                            ADBP3196
      END                                                               ADBP3197
C  **************************************************************       ADBP3198
C                       SUBROUTINE INTEG                                ADBP3199
C  **************************************************************       ADBP3200
      SUBROUTINE INTEG(X,A,B,C,D,XL,XU,SUM,N)                           ADBP3201
C                                                                       ADBP3202
C      INTEGRAL OF A CUBIC SPLINE FUNCTION.                             ADBP3203
C   INPUT:                                                              ADBP3204
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           ADBP3205
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        ADBP3206
C     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.                   ADBP3207
C     N ........................ NUMBER OF GRID POINTS.                 ADBP3208
C     XL ....................... LOWER LIMIT IN THE INTEGRAL.           ADBP3209
C     XU ....................... UPPER LIMIT IN THE INTEGRAL.           ADBP3210
C   OUTPUT:                                                             ADBP3211
C     SUM ...................... VALUE OF THE INTEGRAL.                 ADBP3212
C                                                                       ADBP3213
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP3214
      DIMENSION X(N),A(N),B(N),C(N),D(N)                                ADBP3215
C  ****   SET INTEGRATION LIMITS IN INCREASING ORDER.                   ADBP3216
      IF(XU.GT.XL) THEN                                                 ADBP3217
        XLL=XL                                                          ADBP3218
        XUU=XU                                                          ADBP3219
        SIGN=1.0D0                                                      ADBP3220
      ELSE                                                              ADBP3221
        XLL=XU                                                          ADBP3222
        XUU=XL                                                          ADBP3223
        SIGN=-1.0D0                                                     ADBP3224
      ENDIF                                                             ADBP3225
C  ****  CHECK INTEGRAL LIMITS.                                         ADBP3226
      IWR=0                                                             ADBP3227
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) IWR=1                              ADBP3228
C  ****  FIND INVOLVED INTERVALS.                                       ADBP3229
      SUM=0.0D0                                                         ADBP3230
      CALL FINDI(X,XLL,N,IL)                                            ADBP3231
      CALL FINDI(X,XUU,N,IU)                                            ADBP3232
C  ****  ONLY A SINGLE INTERVAL INVOLVED.                               ADBP3233
      IF(IL.EQ.IU) THEN                                                 ADBP3234
      X1=XLL                                                            ADBP3235
      X2=XUU                                                            ADBP3236
      SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))           ADBP3237
     1   -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))           ADBP3238
      GO TO 2                                                           ADBP3239
      ENDIF                                                             ADBP3240
C  ****  CONTRIBUTIONS FROM DIFFERENT INTERVALS.                        ADBP3241
      X1=XLL                                                            ADBP3242
      X2=X(IL+1)                                                        ADBP3243
      SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))           ADBP3244
     1   -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))           ADBP3245
      IL=IL+1                                                           ADBP3246
      DO 1 I=IL,IU                                                      ADBP3247
      X1=X(I)                                                           ADBP3248
      X2=X(I+1)                                                         ADBP3249
      IF(I.EQ.IU) X2=XUU                                                ADBP3250
      SUMP=X2*(A(I)+X2*((B(I)/2)+X2*((C(I)/3)+X2*D(I)/4)))              ADBP3251
     1    -X1*(A(I)+X1*((B(I)/2)+X1*((C(I)/3)+X1*D(I)/4)))              ADBP3252
    1 SUM=SUM+SUMP                                                      ADBP3253
    2 SUM=SIGN*SUM                                                      ADBP3254
C  ****  INTEGRAL LIMITS OUT OF RANGE.                                  ADBP3255
      IF(IWR.EQ.1) WRITE(6,10)                                          ADBP3256
   10 FORMAT(/'*** WARNING: INTEGRAL LIMITS OUT OF RANGE. ***')         ADBP3257
      RETURN                                                            ADBP3258
      END                                                               ADBP3259
C  **************************************************************       ADBP3260
C                       SUBROUTINE INTEG2                               ADBP3261
C  **************************************************************       ADBP3262
      SUBROUTINE INTEG2(X,A,B,C,D,XL,XU,SUM,N)                          ADBP3263
C                                                                       ADBP3264
C      INTEGRAL OF A SQUARED CUBIC SPLINE FUNCTION.                     ADBP3265
C   INPUT:                                                              ADBP3266
C     X(I) (I=1, ...,N) ........ GRID POINTS.                           ADBP3267
C                    (THE X VALUES MUST BE IN INCREASING ORDER).        ADBP3268
C     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.                   ADBP3269
C     N ........................ NUMBER OF GRID POINTS.                 ADBP3270
C     XL ....................... LOWER LIMIT IN THE INTEGRAL.           ADBP3271
C     XU ....................... UPPER LIMIT IN THE INTEGRAL.           ADBP3272
C   OUTPUT:                                                             ADBP3273
C     SUM ...................... VALUE OF THE INTEGRAL.                 ADBP3274
C                                                                       ADBP3275
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP3276
      DIMENSION X(N),A(N),B(N),C(N),D(N)                                ADBP3277
C  ****   SET INTEGRATION LIMITS IN INCREASING ORDER.                   ADBP3278
      IF(XU.GT.XL) THEN                                                 ADBP3279
        XLL=XL                                                          ADBP3280
        XUU=XU                                                          ADBP3281
        SIGN=1.0D0                                                      ADBP3282
      ELSE                                                              ADBP3283
        XLL=XU                                                          ADBP3284
        XUU=XL                                                          ADBP3285
        SIGN=-1.0D0                                                     ADBP3286
      ENDIF                                                             ADBP3287
C  ****  CHECK INTEGRAL LIMITS.                                         ADBP3288
      IWR=0                                                             ADBP3289
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) IWR=1                              ADBP3290
C  ****  FIND INVOLVED INTERVALS.                                       ADBP3291
      SUM=0.0D0                                                         ADBP3292
      CALL FINDI(X,XLL,N,IL)                                            ADBP3293
      CALL FINDI(X,XUU,N,IU)                                            ADBP3294
C  ****  ONLY A SINGLE INTERVAL INVOLVED.                               ADBP3295
      IF(IL.EQ.IU) THEN                                                 ADBP3296
      X1=XLL                                                            ADBP3297
      X2=XUU                                                            ADBP3298
      SUM=X2*(A(IL)*(A(IL)+X2*B(IL))                                    ADBP3299
     1   +X2*X2*(((2*A(IL)*C(IL)+B(IL)**2)/3.0D0)                       ADBP3300
     2   +X2*(((B(IL)*C(IL)+A(IL)*D(IL))/2.0D0)                         ADBP3301
     3   +X2*(((2*B(IL)*D(IL)+C(IL)**2)/5.0D0)                          ADBP3302
     4   +X2*D(IL)*((C(IL)/3.0D0)+X2*D(IL)/7.0D0)))))                   ADBP3303
     5   -X1*(A(IL)*(A(IL)+X1*B(IL))                                    ADBP3304
     6   +X1*X1*(((2*A(IL)*C(IL)+B(IL)**2)/3.0D0)                       ADBP3305
     7   +X1*(((B(IL)*C(IL)+A(IL)*D(IL))/2.0D0)                         ADBP3306
     8   +X1*(((2*B(IL)*D(IL)+C(IL)**2)/5.0D0)                          ADBP3307
     9   +X1*D(IL)*((C(IL)/3.0D0)+X1*D(IL)/7.0D0)))))                   ADBP3308
      GO TO 2                                                           ADBP3309
      ENDIF                                                             ADBP3310
C  ****  CONTRIBUTIONS FROM DIFFERENT INTERVALS.                        ADBP3311
      X1=XLL                                                            ADBP3312
      X2=X(IL+1)                                                        ADBP3313
      SUM=X2*(A(IL)*(A(IL)+X2*B(IL))                                    ADBP3314
     1   +X2*X2*(((2*A(IL)*C(IL)+B(IL)**2)/3.0D0)                       ADBP3315
     2   +X2*(((B(IL)*C(IL)+A(IL)*D(IL))/2.0D0)                         ADBP3316
     3   +X2*(((2*B(IL)*D(IL)+C(IL)**2)/5.0D0)                          ADBP3317
     4   +X2*D(IL)*((C(IL)/3.0D0)+X2*D(IL)/7.0D0)))))                   ADBP3318
     5   -X1*(A(IL)*(A(IL)+X1*B(IL))                                    ADBP3319
     6   +X1*X1*(((2*A(IL)*C(IL)+B(IL)**2)/3.0D0)                       ADBP3320
     7   +X1*(((B(IL)*C(IL)+A(IL)*D(IL))/2.0D0)                         ADBP3321
     8   +X1*(((2*B(IL)*D(IL)+C(IL)**2)/5.0D0)                          ADBP3322
     9   +X1*D(IL)*((C(IL)/3.0D0)+X1*D(IL)/7.0D0)))))                   ADBP3323
      IL=IL+1                                                           ADBP3324
      DO 1 I=IL,IU                                                      ADBP3325
      X1=X(I)                                                           ADBP3326
      X2=X(I+1)                                                         ADBP3327
      IF(I.EQ.IU) X2=XUU                                                ADBP3328
      SUMP=X2*(A(I)*(A(I)+X2*B(I))                                      ADBP3329
     1    +X2*X2*(((2*A(I)*C(I)+B(I)**2)/3.0D0)                         ADBP3330
     2    +X2*(((B(I)*C(I)+A(I)*D(I))/2.0D0)                            ADBP3331
     3    +X2*(((2*B(I)*D(I)+C(I)**2)/5.0D0)                            ADBP3332
     4    +X2*D(I)*((C(I)/3.0D0)+X2*D(I)/7.0D0)))))                     ADBP3333
     5    -X1*(A(I)*(A(I)+X1*B(I))                                      ADBP3334
     6    +X1*X1*(((2*A(I)*C(I)+B(I)**2)/3.0D0)                         ADBP3335
     7    +X1*(((B(I)*C(I)+A(I)*D(I))/2.0D0)                            ADBP3336
     8    +X1*(((2*B(I)*D(I)+C(I)**2)/5.0D0)                            ADBP3337
     9    +X1*D(I)*((C(I)/3.0D0)+X1*D(I)/7.0D0)))))                     ADBP3338
      IF(SUMP.LT.0.0D0) SUMP=0.0D0                                      ADBP3339
    1 SUM=SUM+SUMP                                                      ADBP3340
    2 SUM=SIGN*SUM                                                      ADBP3341
C  ****  INTEGRAL LIMITS OUT OF RANGE.                                  ADBP3342
      IF(IWR.EQ.1) WRITE(6,10)                                          ADBP3343
   10 FORMAT(/'*** WARNING: INTEGRAL LIMITS OUT OF RANGE. ***')         ADBP3344
      RETURN                                                            ADBP3345
      END                                                               ADBP3346




c
c*******************************************************************
c
	SUBROUTINE AROPM(RP,zPOT,A,B,C,D,npts,rmax)
	IMPLICIT NONE

	INTEGER NPOT,I,NPTS
	DOUBLE PRECISION Z0,ZINF,VLIM,X,Y,RP,zPOT,A,B,C,D,dx
      integer iread
      real*8 rmax, zasymp

	PARAMETER(NPOT=500)

	DIMENSION RP(NPOT),zPOT(NPOT),A(NPOT),B(NPOT),C(NPOT),D(NPOT)

	OPEN (UNIT=2,FILE='charge.dat',STATUS='OLD')
c.... Asumming the input charge goes to 1 asymptotically
      zasymp = 1.d0

      npts = 0
	DO I=1, NPOT 
	 READ(2,*,end=5000)X,Y
c  	 RP(I)=(i-1)*0.05
c 	 zPOT(I)=-1.0d0
       RP(i) = X
       zPOT(i) = -Y
       npts = npts + 1
	enddo
5000  CLOSE(2)
      rp(npts) = rmax
      zpot(npts) = -zasymp
	CALL SPLINE(rp,zPOT,A,B,C,D,0.D0,0.D0,npts)

      return
	END
c
c*******************************************************************
c
	SUBROUTINE AROPMold(RP,POT,A,B,C,D)
	IMPLICIT NONE

	INTEGER NPOT,I
	DOUBLE PRECISION Z0,ZINF,VLIM,X,Y,RP,POT,TOP,A,B,C,D

	PARAMETER(NPOT=400)

	DIMENSION RP(NPOT),POT(NPOT),A(NPOT),B(NPOT),C(NPOT),D(NPOT)

	OPEN (UNIT=2,FILE='potential.inp',STATUS='OLD')
	Z0=18.D0 !Hidrogeno
	ZINF=1.D0
	VLIM=1.D-6
	DO I=1, NPOT-1
	READ(2,*)X,Y
 	RP(I+1)=X !EMPIEZA CON I+1 PORQUE RP(1)=0 (SALVAT)
	POT(I+1)=X*Y/2.D0-Z0
! 	POT(I)=-Z0
	TOP=DABS(POT(I+1)+VLIM)
	IF (TOP .LT. VLIM)POT(I+1)=-ZINF
	END DO
	RP(1)=0.D0
	POT(1)=-Z0
	CLOSE(2)
	CALL SPLINE(RP,POT,A,B,C,D,0.D0,0.D0,NPOT)
	END
c
c*******************************************************************
c
	REAL*8 FUNCTION POTE(XC,RP,POT,A,B,C,D,npts)
	IMPLICIT REAL*8(a-h,o-z) 
	IMPLICIT INTEGER (I-N)
	PARAMETER (NPOT=500)
	REAL*8 a(NPOT),b(NPOT),c(NPOT),d(NPOT),rp(NPOT),pot(NPOT)
      CALL FINDI(RP,XC,npts,IN)
      POTE=A(IN)+B(IN)*XC+C(IN)*XC*XC+D(IN)*XC*XC*XC
	END
c
c*******************************************************************
c


