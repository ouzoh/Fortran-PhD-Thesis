PROGRAM FITCOSINE
IMPLICIT NONE
INTEGER,PARAMETER::Nmax=100 !max data points to store, no of coeffs
INTEGER::I,J,K,wc,Imax,N,Tmax
REAL::Ycal(1,Nmax)! ,thetha(Nmax+1),
CHARACTER:: dummy
INTEGER::Mmax
REAL, ALLOCATABLE, DIMENSION(:) ::deg
REAL, ALLOCATABLE, DIMENSION(:) :: E_J
REAL::sumXTX,suumXTX,SUMy,SIGMA, D_SUM
REAL, PARAMETER :: ONE       = -1.D0
REAL, PARAMETER :: PI        = ACOS(ONE)       ! some constants
REAL, PARAMETER :: Degree180 = 180.0
REAL, PARAMETER :: RtoD    = Degree180/PI
REAL, PARAMETER :: DtoR    = PI/Degree180
character(len=10) :: filename

!#####################################################################!
!Author: Ogaga Uzoh                                                   !
!Uni: University College London                                       !
!Date: 3-3-2013                                                       !
!Part of my published work: Analysis of the conformational profiles of!
!fenamates shows route towards novel, higher accuracy, force-fields   !
!for pharmaceuticals (paper) & Modelling molecular flexibility for    !
!crystal structure prediction (thesis)                                !
!Files needed in conjunction with code:                               !
!                        1. summary.r         -> ab-initio data       !
!                        2. 2 nag libraries                           !
!                             i. g02ddfe      -> linear regression    !
!#####################################################################!




OPEN (40,FILE='summary.r')
!count the no data {Imax}
CALL SYSTEM ('wc -l < summary.r > temp_Imax.txt')
OPEN (41,FILE='temp_Imax.txt')
READ (41,*) Tmax
Imax=Tmax-2
CALL SYSTEM ('rm -f temp_Imax.txt')


ALLOCATE(deg(Imax))
ALLOCATE(E_J(Imax))

!Remove dominant.r file -> summary of the output
CALL SYSTEM('rm -f dominant.r')

!Read degree and ab initio energies
READ (40,*) Mmax
READ (40,*) dummy

DO I = 1, Imax
READ (40,*) deg(I), E_J(I)
END DO


!Header for final output!
 OPEN (217,FILE='temp.txt')
 DO N=1,10 !Imax
     WRITE(217, "(1X,A3,I1,1X)",ADVANCE='NO') "cos", N-1
     !WRITE(217,'(A3)',ADVANCE='NO') trim (filename)
     !WRITE(217,'(1X,3A,I2,A2,1X)',ADVANCE='NO') 'cos',N,'X'
 END DO
 DO N=11,Imax
     WRITE(217, "(1X,A3,I2,1X)",ADVANCE='NO') "cos", N-1
 END DO

!3008 FORMAT(1(1X,3A,I10.2,1A,1X))

 WRITE (217,'(35A)',ADVANCE='NO') "   SErF1   MinEC2      BLnEgy3     RE@GM4    SErI5"
 CLOSE (217)


 CALL SYSTEM ('cat temp.txt > dominant.r')


 CALL SYSTEM ('rm -f  summary_log_input.d')
 DO N=1,Imax
         !For each cosine series
         !perform regression
         OPEN (212,FILE='g02ddfe.d')
         WRITE(212,*)   'INPUT data for multiple regression calculation'
         WRITE(212,*)   Imax
         WRITE(212,3009) Imax,N,"'U'"
         DO J=1,Imax
             WRITE(212,3010) &
                         & (DBLE(COS(DtoR*(K-1)*deg(J))), K=1,N,+1),    &
                         & E_J(J) 
         END DO
         
         CALL SYSTEM ('cat g02ddfe.d >> summary_log_input.d')

         !2.modify the regression code
         OPEN(214,FILE='g02ddfe.f')

         WRITE(214,"(a)" ) "*     G02DDF Example Program Text"
         WRITE(214,"(a)" ) "*     Mark 14 Release. NAG Copyright 1989."
         WRITE(214,"(a)" ) "*     .. Parameters .."
         WRITE(214,"(a)" ) "      INTEGER          MMAX, NMAX, LDQ, PMAX"
         WRITE(214,' (A29,I3,A6,I3,A1)') "      PARAMETER        (MMAX=",N+1,",NMAX=",Imax,")"
         WRITE(214,"(a)" ) "      PARAMETER        (LDQ=NMAX)"
         WRITE(214,"(a)" ) "      INTEGER          NIN, NOUT"
         WRITE(214,"(a)" ) "      PARAMETER        (NIN=5,NOUT=6)"
         WRITE(214,"(a)" ) "*     .. Local Scalars .."
         WRITE(214,"(a)" ) "      DOUBLE PRECISION RSS, TOL,MATS,SumF,SumT"
         WRITE(214,"(a)" ) "      INTEGER          I, IDF, IFAIL, IP, IRANK, J, M, N, K"
         WRITE(214,"(a)" ) "      INTEGER          Param(MMAX,1)"
         WRITE(214,"(a)" ) "      LOGICAL          SVD"
         WRITE(214,"(a)" ) "      CHARACTER        WEIGHT"
         WRITE(214,"(a)" ) "*     .. Local Arrays .."
         WRITE(214,"(a)" ) "      DOUBLE PRECISION B(MMAX), COV(MMAX*(MMAX+1)/2), P(MMAX*(MMAX+2)),"
         WRITE(214,"(a)" ) "     +                 Q(LDQ,MMAX+1), SE(MMAX), WK(MMAX*MMAX+5*MMAX),"
         WRITE(214,"(a)" ) "     +                 WT(NMAX), X(NMAX,MMAX),EMjODiEL(LDQ,MMAX),"
         WRITE(214,"(a)" ) "     +                 QA(LDQ,MMAX+1),eJODI(LDQ,MMAX),E_J(LDQ,MMAX),"
         WRITE(214,"(a)" ) "     +                 deg(LDQ,MMAX), BP(MMAX+2) "
         WRITE(214,"(a)" ) "*     .. External Subroutines .."
         WRITE(214,"(a)" ) "      EXTERNAL         G02DDF, G02DEF"
         WRITE(214,"(a)" ) "*     .. Executable Statements .."
         WRITE(214,"(a)" ) "*      WRITE (NOUT,*) 'G02DDF Example Program Results'"
         WRITE(214,"(a)" ) "*     Skip heading in data file"
         WRITE(214,"(a)" ) "      READ (NIN,*)"
         WRITE(214,"(a)" ) "      WRITE(NOUT,*) 'SUMMARY FOR', IP, 'COSINE '"
         WRITE(214,"(a)" ) "      READ (NIN,*) PMAX"
         WRITE(214,"(a)" ) "      READ (NIN,*) N, M, WEIGHT"
         WRITE(214,"(a)" ) "      IF (N.LE.NMAX .AND. M.LT.MMAX) THEN"
         WRITE(214,"(a)" ) "         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN"
         WRITE(214,"(a)" ) "            DO 20 I = 1, N"
         WRITE(214,"(a)" ) "               READ (NIN,*) (X(I,J),J=1,M), Q(I,1), WT(I)"
         WRITE(214,"(a)" ) "   20       CONTINUE"
         WRITE(214,"(a)" ) "         ELSE"
         WRITE(214,"(a)" ) "            DO 40 I = 1, N"
         WRITE(214,"(a)" ) "               READ (NIN,*) (X(I,J),J=1,M), Q(I,1)"
         WRITE(214,"(a)" ) "   40       CONTINUE"

         WRITE(214,"(a)" ) "            DO I = 1, N"
         WRITE(214,"(a)" ) "               DO J=1,M"
         WRITE(214,"(a)" ) "                      E_J(I,1)=Q(I,1)"
         WRITE(214,"(a)" ) "               END DO"
         WRITE(214,"(a)" ) "            END DO"
         WRITE(214,"(a)" ) "         END IF"
         WRITE(214,"(a)" ) "*        Set tolerance"
         WRITE(214,"(a)" ) "         TOL = 0.000001D0"
         WRITE(214,"(a)" ) "         IP = 0"
         WRITE(214,"(a)" ) "         DO 60 I = 1, M"
         WRITE(214,"(a)" ) "            IFAIL = 1"
         WRITE(214,"(a)" ) "*"
         WRITE(214,"(a)" ) "*           Fit model using G02DEF"
         WRITE(214,"(a)" ) "            CALL G02DEF(WEIGHT,N,IP,Q,LDQ,P,WT,X(1,I),RSS,TOL,IFAIL)"
         WRITE(214,"(a)" ) "*"
         WRITE(214,"(a)" ) "            IF (IFAIL.EQ.0) THEN"
         WRITE(214,"(a)" ) "               IP = IP + 1"
         WRITE(214,"(a)" ) "            ELSE IF (IFAIL.EQ.3) THEN"
         WRITE(214,"(a)" ) "               WRITE (NOUT,*) ' * New variable not added *'"
         WRITE(214,"(a)" ) "            ELSE"
         WRITE(214,"(a)" ) "               WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "               WRITE (NOUT,99996) ' ** G02DEF returned with IFAIL = ',"
         WRITE(214,"(a)" ) "     +           IFAIL"
         WRITE(214,"(a)" ) "               GO TO 100"
         WRITE(214,"(a)" ) "            END IF"
         WRITE(214,"(a)" ) "   60    CONTINUE"
         WRITE(214,"(a)" ) "         RSS = 0.0D0"
         WRITE(214,"(a)" ) "         IFAIL = 1"
         WRITE(214,"(a)" ) "*"
         WRITE(214,"(a)" ) "         CALL G02DDF(N,IP,Q,LDQ,RSS,IDF,B,SE,COV,SVD,IRANK,P,TOL,WK,"
         WRITE(214,"(a)" ) "     +               IFAIL)"
         WRITE(214,"(a)" ) "*"
         WRITE(214,"(a)" ) "         IF (IFAIL.EQ.0) THEN"
         WRITE(214,"(a)" ) "            WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "            IF (SVD) THEN"
         WRITE(214,"(a)" ) "               WRITE (NOUT,*) 'Model not of full rank'"
         WRITE(214,"(a)" ) "               WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "            END IF"
         WRITE(214,"(a)" ) "*     !Addition  : Calculate the intra energy of model"
         WRITE(214,"(a)" ) "*     !Author    : Ogaga G. Uzoh"
         WRITE(214,"(a)" ) "*     !University: UCL"
         WRITE(214,"(a)" ) "*     !Date      : 29-11-2013"
         WRITE(214,"(a)" ) "*"

         !The model's relative and absolute intramolecular energy
         WRITE(214,"(a)" ) "            WRITE (NOUT,*)'Intra model,   Energy+constant' "
         WRITE(214,"(a)" ) "         DO J=1,N"
         WRITE(214,"(a)" ) "            mats=0.0"         
         !1. Calculating the model's absolute intramolecular energy
         !1a.Calculating the Coeff*Regressors (Beta*X)
         WRITE(214,"(a)" ) "            DO I =1,IP"
         WRITE(214,"(a)" ) "                mats=mats+X(J,I)*B(I)"
         WRITE(214,"(a)" ) "            END DO"
         !1b.Calculating the Coeff*Regressor + contributions not fitted(QA)
         WRITE(214,"(a)" ) "            QA(J,1)=0"
         WRITE(214,"(a)" ) "            EMjODiEL(J,1)=QA(J,1)+mats"
         WRITE(214,"(a)" ) "         END DO"
         !Relative intramolecular energy = eJodi
         WRITE(214,"(a)" ) "         DO J=1,N"
         WRITE(214,"(a)" ) "            eJODI(J,1)=EMjODiEL(J,1)-MINVAL(EMjODiEL(1:N,1))"
         !Absolute intramolecular energy
         WRITE(214,"(a)" ) "            WRITE(*,'(F12.2,A1,F12.2)') eJODI(J,1),','"!,"
         !Testing 
!         WRITE(214,"(a)" ) "     +EMjODiEL(J,1)-0"!B(1)"
         WRITE(214,"(a)" ) "         END DO"
         
         WRITE(214,"(a)" ) "         mats=0.0"
         WRITE(214,"(a)" ) "         DO J=1,N"
         WRITE(214,"(a)" ) "            mats=mats+(E_J(J,1)-eJODI(J,1))**2"
         WRITE(214,"(a)" ) "         END DO"
         WRITE(214,"(a)" ) "            WRITE (NOUT,99999) 'Residual sum of squares (fit)  = ', RSS"
         
         !OUTPUTS
         !#standard error fit          (1)
         WRITE(214,"(a)" ) "            WRITE (NOUT,99999) 'Standard error (fit)           = ',"
         WRITE(214,"(a)" ) "     +(RSS/(IDF))**0.5" !#1
         WRITE(214,"(a)" ) "            WRITE (NOUT,99998) 'Degrees of freedom             = ', IDF"
         
         !#minimum energy conformation (2)
         WRITE(214,"(a)" ) "         WRITE(*,*)  'Minimum energy conformation'"
         WRITE(214,"(a)" ) "         WRITE(*,*) deg(MINLOC(eJODI(1:N,1)),1)" !#2
         
         !#baseline energy             (3)
         WRITE(214,"(a)" ) "         WRITE(*,*)  'Baseline energy'"
         WRITE(214,"(a)" ) "         WRITE(*,*) ((MINVAL(EMjODiEL(1:N,1)))-B(1))" !#3

         !#model's relative energy @ global minimum (4)
         WRITE(214,"(a)" ) "         WRITE(*,*)  'Model relative energy at global minimum'"
         WRITE(214,"(a)" ) "         WRITE(*,*) eJODI((MINLOC(E_J(1:N,1))),1)" !#4
         
         !#standard error intra        (5)
         WRITE(214,"(a)" ) "            WRITE (NOUT,99999) 'Residual sum of squares (intra)= ', mats"
         WRITE(214,"(a)" ) "            WRITE (NOUT,99999) 'Standard error (intra)         = ',"
         WRITE(214,"(a)" ) "     +(mats/(IDF))**0.5" !#5        


         WRITE(214,"(a)" ) "            WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "            WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "     +        'Variable,  Parameter estimate,  Standard error'"
         
         !#rescale parameters          (0)
!!!         WRITE(214,"(a)" ) "            OPEN (215,FILE='Rpara.txt')"
         WRITE(214,"(a)" ) "            OPEN (216,FILE='Table.txt')"
         WRITE(214,"(a)" ) "            DO 80 J = 1, IP"
         WRITE(214,"(a)" ) "               WRITE (NOUT,99997) J-1,',', B(J),',', SE(J)"
!!!         WRITE(214,"(a)" ) "               WRITE (215,'(1X,F20.4)',ADVANCE='NO') B(J)"
         WRITE(214,"(a)" ) "   80       CONTINUE"
         WRITE(214,"(a)" ) "            ELSEIF (J.EQ.1) THEN"
         WRITE(214,"(a)" ) "               WRITE (NOUT,99995) 'c',',',B(J),',', SE(J)"
!!!         WRITE(214,"(a)" ) "               WRITE (215,'(1X,F20.4)',ADVANCE='NO') B(J)"
         WRITE(214,"(a)" ) "            ELSE"
         WRITE(214,"(a)" ) "            WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "            WRITE (NOUT,99996) ' ** G02DDF returned with IFAIL = ',"
         WRITE(214,"(a)" ) "     +        IFAIL"
         WRITE(214,"(a)" ) "         END IF"
         WRITE(214,"(a)" ) "      END IF"
         WRITE(214,"(a)" ) "  100 CONTINUE"

         !WRITE FORMATED PARAMETER IN A SINGLE LINE
         WRITE(214,"(a)" ) "             WRITE (NOUT,*) '#SUMMARY#'"
         
         WRITE(214,"(a)" ) "             DO 81 I=1,IP"
         WRITE(214,"(a)" ) "                WRITE (216,'(1X,F10.4)',ADVANCE='NO') B(I)"
         WRITE(214,"(a)" ) "   81        CONTINUE"
         
         WRITE(214,"(a)" ) "             DO 82 I=IP+1,PMAX"
         WRITE(214,"(a)" ) "             WRITE (216,'(2X,1A,2X)',ADVANCE='NO') '-'"
         WRITE(214,"(a)" ) "   82        CONTINUE"

         WRITE(214,"(a)" ) "               WRITE (216,'(1X,5F10.2)')"
         !#standard error fit          (1)
         WRITE(214,"(a)" ) "     +(RSS/(IDF))**0.5,"
         !#minimum energy conformation (2)
         WRITE(214,"(a)" ) "     +deg(MINLOC(eJODI(1:N,1)),1)," !#2
         !#baseline energy             (3)
         WRITE(214,"(a)" ) "     +((MINVAL(EMjODiEL(1:N,1)))-B(1))," !#3
         !#model's relative energy @ global minimum (4)
         WRITE(214,"(a)" ) "     +eJODI((MINLOC(E_J(1:N,1))),1)," !#4
         !#standard error intra        (5)
         WRITE(214,"(a)" ) "     +(mats/(IDF))**0.5" !#5
         
         WRITE(214,"(a)" ) "              CALL SYSTEM ('cat Table.txt >> dominant.r')"     

         WRITE(214,"(a)" ) "                 WRITE (NOUT,*) ''"        
         
         WRITE(214,"(a)" ) "*"
         WRITE(214,"(a)" ) "99999 FORMAT (1X,A,2F16.10)"
         WRITE(214,"(a)" ) "99998 FORMAT (1X,A,I4)"
         WRITE(214,"(a)" ) "99997 FORMAT (1X,I6,A1,F20.4,A1,F20.4)"
         WRITE(214,"(a)" ) "99996 FORMAT (1X,A,I5)"
         WRITE(214,"(a)" ) "99995 FORMAT (1X,A6,A1,F20.4,A1,F20.4)"
         WRITE(214,"(a)" ) "      END"
            

         
         !
         CALL SYSTEM ('gfortran -m64 g02ddfe.f /usr/local/nag/fll6a22dfl/lib/libnag_nag.a &
                      & -o g02ddfe.exe')
         CALL SYSTEM ('./g02ddfe.exe < g02ddfe.d ')!> g02ddfe.r')
         REWIND(212)
         REWIND(214)

 END DO
 !clean up
CALL SYSTEM('rm -f g02ddfe.exe  summary_log_input.d  Table.txt')
CALL SYSTEM('rm -f g02ddfe.d g02ddfe.f temp.txt')
CALL SYSTEM('echo "---------------------------------------"')
 !end message
CALL SYSTEM('echo "Fitting routine terminated successfully"')
CALL SYSTEM('echo "---------------------------------------"')
CALL SYSTEM('echo ""')
 !Formats
3000 FORMAT(F10.8,G18.8,G18.8,G18.8,I6,A5,A5,I4,I4)
3001 FORMAT(F10.8,G18.8,G18.8,G18.8,I6,A5,A5,I4,I4)
3002 FORMAT(A8,1X,A10,6X,A10,15X,A10,4X,A6,A5,1X,A4,A4,A4)
3003 FORMAT(A9,1X,A13,6X,A13,3X,A9,3X,A11,5X,A11)
3004 FORMAT(G16.10,G16.10,G16.10,G16.10,2x,G14.8,G14.8)
3005 FORMAT(A2,I2,A1,A2,I2,F5.1,4(F8.2,F8.2,A2,F5.2,A1))
3006 FORMAT(A9,A5,12(A8))
3007 FORMAT(100F14.7)
3009 FORMAT(2(I5),A5) !Integers
3010 FORMAT(100(F14.7))
3011 FORMAT(100(A5,A5,A4))

!!!!!!!!!!!!!IMPROVEMENT: DATA EXTRACTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM FITCOSINE
