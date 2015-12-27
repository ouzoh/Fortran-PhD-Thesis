PROGRAM read_mol2
IMPLICIT NONE
INTEGER                                  :: I,J,IOSTAT,L,K,N,IK
INTEGER                                  :: Imax,M,Gmax,Ncon,Mmax
INTEGER                                  :: RI,RJ,RL,RK
INTEGER                                  :: sumL,NcombL !no of regression/comb
CHARACTER(LEN=20)                        :: dummyC
CHARACTER(LEN=30)                        :: mol2
CHARACTER(LEN=10)                        :: degr
REAL                                     :: E_H,dE_H,dE_J
REAL                                     :: R_DIFF
REAL                                     :: SUM_GAV
INTEGER                                  :: int_I, int_K
REAL                                     :: DX,DY,DZ,E_POT !TO CHECK POTENTIAL ENRGY VALUE
REAL                                     :: E_POTD,E_POTR,E_POTD2,E_POTR2
REAL,ALLOCATABLE,DIMENSION(:)            :: paraA,paraB,paraC
REAL,ALLOCATABLE,DIMENSION(:)            :: deg,E_J
REAL,ALLOCATABLE,DIMENSION(:,:)          :: X,Y,Z
REAL,ALLOCATABLE,DIMENSION(:,:)          :: A,B,C
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: R
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: E_GAV
REAL,ALLOCATABLE,DIMENSION(:)            :: E_GAV_ND,E_GAVR_ND,E_GAVD_ND,E_GAVR2
REAL,ALLOCATABLE,DIMENSION(:)            :: E_GAV2,E_GAVD2
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: E_GAVR,E_GAVD
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: E_GAVR3
INTEGER                                  :: at_at
REAL,ALLOCATABLE,DIMENSION(:,:)          :: neigh,E_GAVR_at_Type
INTEGER,ALLOCATABLE                      :: I_st(:), K_st(:),c_int(:)
CHARACTER(LEN=5),ALLOCATABLE,DIMENSION(:):: atom_ty
CHARACTER(LEN=2),ALLOCATABLE,DIMENSION(:):: atom_ele
CHARACTER(LEN=2)                         :: atom_ele1,atom_ele2
CHARACTER(LEN=2),ALLOCATABLE,DIMENSION(:):: atom_G_ele1,atom_G_ele2
CHARACTER(LEN=2),ALLOCATABLE,DIMENSION(:):: atom_D_ele1,atom_D_ele2
LOGICAL                                  :: expr1,expr2,expr3,expr4,expr12,expr34,expr42
LOGICAL                                  :: exprh1,exprh2,exprh3,exprh4
LOGICAL                                  :: ans1=.TRUE., ans2=.FALSE.

!Energy cutoff
!If the difference between interactions over all conformations varies then
!the interaction is classed as dominant
REAL                                     :: E_cutoff=1.00

!#####################################################################!
!Author: Ogaga Uzoh                                                   !
!Uni: University College London                                       !
!Date: 3-3-2013                                                       !
!Part of my published work: Analysis of the conformational profiles of!
!fenamates shows route towards novel, higher accuracy, force-fields   !
!for pharmaceuticals (paper) & Modelling molecular flexibility for    !
!crystal structure prediction (thesis)                                !
!Files needed in conjunction with code:                               !
!                        1. summary.r         -> ab-initio            !
!                        2. mol2 files        -> conformations        !
!                        3. gavezzotti_para.d -> Gavezotti's parameter!
!                        4. 2 nag libraries                           !
!                             i. g02ddfe      -> linear regression    !
!                             ii.h03adfe      -> nearest neigbour     !
!Code also contains subrountine combi.f90 for combitorial mathematics !
!#####################################################################!

!open summary.r, read and store data
 CALL SYSTEM ('cp -r inputs/* .')
 CALL SYSTEM ('wc -l < summary.r > temp_n')

 OPEN (200,FILE='temp_n')
 !number of conformations
 READ (200,*) Imax
 Imax=Imax-2
 CLOSE(200)

 CALL SYSTEM ('wc -l < gavezzotti_para.d > temp_n')
 OPEN (200,FILE='temp_n')
 !number of parameters
 READ (200,*) Gmax
 Gmax=Gmax-2
 CLOSE(200)
 
!Allocations
 ALLOCATE(deg(Imax))
 ALLOCATE(E_J(Imax))

 OPEN(UNIT=201,FILE='summary.r')
 !M is the number of atoms
 READ(201,*) M
 READ(201,*) dummyC

 ALLOCATE(atom_ele(M))             ! atom element
 ALLOCATE(atom_G_ele1(Gmax))       ! atom's Gavezzotti's element
 ALLOCATE(atom_G_ele2(Gmax))       !  ,,
 ALLOCATE(atom_D_ele1(Gmax))       !  ,,
 ALLOCATE(atom_D_ele2(Gmax))       !  ,,
 ALLOCATE(paraA(Gmax))             ! parameter A
 ALLOCATE(paraB(Gmax))             !           B
 ALLOCATE(paraC(Gmax))             !           C
 ALLOCATE(X(Imax,M))               ! cartesian
 ALLOCATE(Y(Imax,M))               !   not
 ALLOCATE(Z(Imax,M))               !   stored
 ALLOCATE(A(M-1,M))                ! para A  !stored
 ALLOCATE(B(M-1,M))                ! para B  !  once
 ALLOCATE(C(M-1,M))                ! para C  !  A,B,C is the same for one atom_ty
 ALLOCATE(atom_ty(M))              ! atom_type
 ALLOCATE(R(Imax,M-1,M))           ! interatomic distance
 ALLOCATE(neigh(M-1,M))            !nearest neighbours
 ALLOCATE(E_GAV(Imax,M-1,M))       !Gav. Energy
 ALLOCATE(E_GAVR(Imax,M-1,M))      !Gav. Repulson Energy
 ALLOCATE(E_GAVD(Imax,M-1,M))      !Gav. Dispersion Energy



! ALLOCATE(at_at(Imax*Imax,(M-1)*(M),M*M))
 ALLOCATE(E_GAVR_at_Type(Imax,Imax))


 ALLOCATE(E_GAV_ND(Imax))
 ALLOCATE(E_GAVR_ND(Imax))
 ALLOCATE(E_GAVD_ND(Imax))

 ALLOCATE(E_GAV2(Imax))
 ALLOCATE(E_GAVR3(Imax,Imax+1,Imax))
 ALLOCATE(E_GAVR2(Imax))
 ALLOCATE(E_GAVD2(Imax))

 ALLOCATE(I_st(Mmax))
 ALLOCATE(K_st(Mmax))

 ALLOCATE(c_int(Gmax))

 !Remove dominant.r file
 CALL SYSTEM('rm -f dominant.r')
 DO J=1,Imax
     !read summary
     !contains deg
     !reads intramolecular energies in H and KJ/mol & its derivatives
     READ(201,*,END=100) deg(J),E_H,E_J(J),dE_H,dE_J
     WRITE(degr,'(F10.0)') deg(J)

     !open each mol2 file
     mol2=adjustl("zmat")//trim(adjustl(degr))//adjustl("mol2")
     WRITE (*,*) mol2
     OPEN(202,FILE=mol2)

     !loop over each mol
     DO I=1,7
         READ(202,*) dummyC
     END DO
     !Allocate the atom's element and type and the cartesian coord
     !  for each conformation present
     WRITE(*,*) "Conformation",J
     DO I=1,M
         READ(202,*) dummyC, atom_ele(I),X(J,I),Y(J,I),Z(J,I) &
         &,atom_ty(I),dummyC,dummyC,dummyC
         !WRITE(*,*) "Conformation",J
      !   WRITE(*,*) atom_ele(I),X(J,I),Y(J,I),Z(J,I),atom_ty(I)
     END DO
     CLOSE(202)
     100 CONTINUE
 END DO
 CLOSE(201)

 !Find nearest neighbour
 !Open only one mol2 file
 WRITE(degr,'(F10.0)') deg(1)
 mol2=adjustl("zmat")//trim(adjustl(degr))//adjustl("mol2")
 OPEN(202,FILE=mol2)
 !loop over the mol2 file
 DO I=1,8+M
     READ(202,*) dummyC
 END DO

 !extract connectivity table
 OPEN (203,FILE='temp_con')
 OPEN (2003,FILE='conn_atm')
 DO I=1,M*(M-1)
     READ (202,*,END=101) RL,RJ,RK,dummyC
     !format must be adjusted for atom numbers greater than 1,000
     WRITE (203,'(I4,I4)') RJ,RK
     !Find polar Hydrogen
     !N-H,O-H incorporated
     !Don't include H-S
     WRITE (2003,*) atom_ele(RJ),atom_ele(RK)
     exprh1=(atom_ele(RJ).EQ.'H').AND.(atom_ele(RK).EQ.'O')
     exprh2=(atom_ele(RJ).EQ.'O').AND.(atom_ele(RK).EQ.'H')
     exprh3=(atom_ele(RJ).EQ.'H').AND.(atom_ele(RK).EQ.'N')
     exprh4=(atom_ele(RJ).EQ.'N').AND.(atom_ele(RK).EQ.'H')
     IF (exprh1.EQV.ans1) THEN
         !polar H present
         WRITE(*,*) RJ,'is polar'
         atom_ele(RJ)='HB'
     ELSEIF (exprh2.EQV.ans1) THEN
         !polar H present
         WRITE(*,*) RK,'is polar'
         atom_ele(RK)='HB'
     ELSEIF (exprh3.EQV.ans1) THEN
         !polar H present
         WRITE(*,*) RJ,'is polar'
         atom_ele(RJ)='HB'
     ELSEIF (exprh4.EQV.ans1) THEN
         !polar H present
         WRITE(*,*) RK,'is polar'
         atom_ele(RK)='HB'
     ENDIF
 101 CONTINUE
 END DO
 CLOSE(202)

 !find the number of non-zero element in the distance matrix D
 !This equals the number of connectivity in the connectivity table
 CALL SYSTEM ('wc -l < temp_con > temp_Cnum')
 OPEN (204,FILE='temp_Cnum')
 READ (204,*) Ncon
 CLOSE (204)

 !Number of atoms=M
 !number of connectivity=Ncon
 !Start atom and end atom determined by loop I,J
 !I .NE. J

 !Reorder the conn file
 REWIND(203)
 OPEN (205,FILE='temp_conR')
 DO I=1,Ncon
     READ(203,'(I4,I4)') RJ,RK
     IF (RK.GT.RJ) THEN
         WRITE (205,'(A3,I4,I4)') '1.0',RJ,RK
     ELSE
         WRITE (205,'(A3,I4,I4)') '1.0',RK,RJ
     END IF
 END DO
 CLOSE(203)


 !Create file to be used by h03adf subrountine
 !uses temp_conR
 !M,I,J,Ncon,F
 OPEN (208,FILE='conn_table')
 DO I=1,M-1
     DO K=I+1,M
         OPEN (206,FILE='temp_top')
         WRITE(206,*) 'This is the input required to find neigh'
         WRITE(206,*) M,I,K,Ncon,'F'
         CALL SYSTEM ('cat temp_top temp_conR > h03adfe.d')
         REWIND (206)
         !'Length of shortest path using nag library h03adfe.d'
         CALL SYSTEM ('./h03adfe.exe < h03adfe.d > h03adfe.r')
         !grep results
         CALL SYSTEM ("grep 'Length of shortest path' h03adfe.r > conn.r")
         !awk nearest neighbour value
         CALL SYSTEM ("awk '{print $6}' conn.r > conn1.r")
         OPEN (207,FILE='conn1.r')
         READ (207,*) neigh(I,K)
        ! WRITE(208,'(1X,I4,1X,I4,3X,F8.0)') I, K, neigh(I,K)
         REWIND (207)
     END DO
 END DO


 !Find all intramolecular distances and store
 DO J=1,Imax
     WRITE(208,*) "Conformation",J
     DO I=1,M-1
         DO K=I+1,M
             DX=X(J,I)-X(J,K)
             DY=Y(J,I)-Y(J,K)
             DZ=Z(J,I)-Z(J,K)
             R(J,I,K)=SQRT((DX)**2+(DY)**2+(DZ)**2)
             !WRITE(*,*) I,K,R(J,I,K),neigh(I,K)
             WRITE(208,*) atom_ele(I),atom_ele(K),'/',R(J,I,K),'/',neigh(I,K)
         END DO
     END DO
 END DO


 !read all Gavezzotti's parameter
 OPEN(203,FILE='gavezzotti_para.d')
 READ(203,*) dummyC
 READ(203,*) dummyC
 DO J=1,Gmax
     READ(203,*) atom_G_ele1(J), atom_G_ele2(J), paraA(J),paraB(J),paraC(J)
     !WRITE(*,*) atom_G_ele1(J), atom_G_ele2(J), paraA(J),paraB(J),paraC(J)
 END DO

 !use the read Gavezzotti's parameter to determine all intramolecular parameters
 DO I=1,M-1
     DO K=I+1,M
        ! WRITE(*,*) 'at_I=',atom_ele(I),'at_J=',atom_ele(K)
         DO J=1,Gmax
            expr1=(atom_ele(I).EQ.atom_G_ele1(J)) .AND. (atom_ele(K).EQ.atom_G_ele2(J))
            expr2=(atom_ele(K).EQ.atom_G_ele1(J)) .AND. (atom_ele(I).EQ.atom_G_ele2(J))
            expr12=(expr1.OR.expr2)
            expr3=(atom_ele(I).NE.atom_G_ele1(J)) .AND. (atom_ele(K).NE.atom_G_ele2(J))
            expr4=(atom_ele(K).NE.atom_G_ele1(J)) .AND. (atom_ele(I).NE.atom_G_ele2(J))
            expr34=(expr3.OR.expr4)
            IF (expr12.EQV.ans1) THEN
                A(I,K)=paraA(J)
                B(I,K)=paraB(J)
                C(I,K)=paraC(J)
         !       WRITE(*,*) atom_ele(I),atom_ele(K), paraA(J),paraB(J),paraC(J)
                EXIT
            ELSEIF  (expr34.EQV.ans2) THEN
                CYCLE !WRITE(*,*) 'not yet found'
            END IF
        END DO
     END DO
 END DO

 !We want to determine the below commented line!
 !#############################################!
 !evaluate exp-6 atom atom potential
 !#############################################!

 DO J=1,Imax
     E_GAVR(J,I,K)=0
     E_GAVD(J,I,K)=0
     E_GAV(J,I,K)=0
     DO I=1,M-1 !SUM OVER ATOM 
         DO K=I+1,M !SUM OVER ATOM J
             !Condition 1: 1-4 interactions and above

             IF (neigh(I,K).GE.3) THEN
                 !Calculate the repulsion contribution
                 E_GAVR(J,I,K)=(A(I,K)*EXP(-B(I,K)*(R(J,I,K))))
                 
                 !Calculate the dispersion contribution
                 E_GAVD(J,I,K)=-(C(I,K)/(R(J,I,K)**6))

                 !Calculate the repulsion & dispersion dispersion
                 E_GAV(J,I,K)=E_GAVR(J,I,K)+E_GAVD(J,I,K)
                 !WRITE(*,*) J,I,K,E_GAV(J,I,K),neigh(I,K)
             ELSE
                 CYCLE
             END IF
         END DO
     END DO
 END DO

 !Formats
3000 FORMAT(F10.8,G18.8,G18.8,G18.8,I6,A5,A5,I4,I4) 
3001 FORMAT(F10.8,G18.8,G18.8,G18.8,I6,A5,A5,I4,I4)
3002 FORMAT(A8,1X,A10,6X,A10,15X,A10,4X,A6,A5,1X,A4,A4,A4)
3003 FORMAT(A9,1X,A13,6X,A13,3X,A9,3X,A11,5X,A11)
3004 FORMAT(G16.10,G16.10,G16.10,G16.10,2x,G14.8,G14.8)
3005 FORMAT(A2,I2,A1,A2,I2,F5.1,4(F8.2,F8.2,A2,F5.2,A1))
3006 FORMAT(A9,A5,12(A8))
3007 FORMAT(100F14.7)
3008 FORMAT(50(A2,A1,A2,A1,3X))
3009 FORMAT(2(I5),A5) !Integers
3010 FORMAT(A5,100(F14.7))
3011 FORMAT(100(A5,A5,A4))




 PRINT*,''
 PRINT*,'####################################################################'
 PRINT*,'############DOMINANT atomType---atomType interactions###############'
 PRINT*,'####################################################################'
 PRINT*,''


OPEN(209,FILE='dominant.txt')
 WRITE(*,3006) 'at I-at K','N','Rmn','Rmx','Rdiff','ERmn','ERmx','ERdiff','EDmn','EDmx','EDdiff','ERDmn','ERDmx','ERDdiff' 
 DO I=1,M-1
     DO K=I+1,M
         IF (neigh(I,K).GE.3) THEN
             !1. Determine the dominant interactions!
             !The energy cutoff it used             !
             IF ((MAXVAL(E_GAV(1:Imax,I,K)))-(MINVAL(E_GAV(1:Imax,I,K))).GE.E_cutoff) THEN
                 WRITE(*,3005) atom_ele(I),I,'-',atom_ele(K),K,neigh(i,k)+1    &
                 !Rmin,Rmax,Rmax-Rmin
                 &,MINVAL(R(1:Imax,I,K)),MAXVAL(R(1:Imax,I,K)),                &
                 &'(',MAXVAL(R(1:Imax,I,K))-MINVAL(R(1:Imax,I,K)),')'          &
                 !ERepMin,ERepMax,ERepMax-ERepMin
                 &,MINVAL(E_GAVR(1:Imax,I,K)),MAXVAL(E_GAVR(1:Imax,I,K)),      &
                 &'(',MAXVAL(E_GAVR(1:Imax,I,K))-MINVAL(E_GAVR(1:Imax,I,K)),')'&
                 !EDisMin,EDisMax,EDisMax-EDisMin
                 &,MINVAL(E_GAVD(1:Imax,I,K)),MAXVAL(E_GAVD(1:Imax,I,K)),      &
                 &'(',MAXVAL(E_GAVD(1:Imax,I,K))-MINVAL(E_GAVD(1:Imax,I,K)),')'&
                 !ERDMin,ERDMax,ERDMax-ERDMin
                 &,MINVAL(E_GAV(1:Imax,I,K)),MAXVAL(E_GAV(1:Imax,I,K)),        &
                 &'(',MAXVAL(E_GAV(1:Imax,I,K))-MINVAL(E_GAV(1:Imax,I,K)),')'
                 
                 !2a. atomtype---atomtype interactions of the dominant
                 !interactions
                 IF (atom_ele(K).GT.atom_ele(I)) THEN
                     WRITE(209,*) atom_ele(I),' ',atom_ele(K)
                 ELSE
                     WRITE(209,*) atom_ele(K),' ',atom_ele(I)
                 END IF
             ELSE
                 CYCLE
             END IF
         ELSE
             CYCLE
         END IF
     END DO
 END DO

 !2b. Unique atomtype---atomtype interactions of the dominant
 CALL SYSTEM ('sort -u dominant.txt > dominant_at-at.txt')
 CALL SYSTEM ('wc -l < dominant_at-at.txt > temp_n')

 
 


 OPEN (210,FILE='dominant_at-at.txt')
 !number of atType-atType interactions
 OPEN (200,FILE='temp_n')
 READ (200,*) Mmax
 CLOSE(200)
 
 !store dominant interaction
 DO N=1,Mmax
     READ(210,*) atom_D_ele1(N),atom_D_ele2(N)
 END DO 


 !Header for final output!
  OPEN (217,FILE='temp.txt')
  DO N=1,Mmax
      WRITE(217,3008,ADVANCE='NO') trim(atom_D_ele1(N)),'-',trim(atom_D_ele2(N))
  END DO

  WRITE (217,'(35A)',ADVANCE='NO') "   cos      const       SErF1   MinEC2      BLnEgy3     RE@GM4    SErI5"
  CLOSE (217)


  CALL SYSTEM ('cat temp.txt > dominant.r')
 ! CALL SYSTEM ('rm temp.txt')

 !Determine all inputs required for regression analysis
 !Condition 1 applies
 
 DO J=1,Imax
     DO I=1,M-1 !SUM OVER ATOM 
         DO K=I+1,M !SUM OVER ATOM J
             expr2=.FALSE.
             IF (neigh(I,K).GE.3) THEN
                 !REWIND (210)
                 DO N=1,Mmax
                     
                     !#######################################!
                     !check if the interaction is significant!
                     !i.e included in major_int file         !
                     !#######################################!

                     !#######################################!
                     !true if significant                    !
                     !#######################################!
                     !READ(210,*) atom_D_ele1(N), atom_D_ele2(N)
                     expr12=(atom_D_ele1(N).EQ.atom_ele(I)) .AND. (atom_D_ele2(N).EQ.atom_ele(K))
                     expr34=(atom_D_ele2(N).EQ.atom_ele(I)) .AND. (atom_D_ele1(N).EQ.atom_ele(K))
                     expr42=expr12.OR.expr34
                     IF (expr42.EQV.ans1) THEN
                         !I,K contribute significantly
                         expr1=.TRUE.
                         at_at=N
                         !!!!!!!!!!!!!!!!WRITE(*,*) 'tHIS IS N VALUE',N
                         expr2=expr2.OR.expr1
                     ELSE
                         !########################!
                         !false if not significant!
                         !########################!
                         !WRITE(*,*) atom_ele(I),atom_ele(K),'F'
                         expr3=.FALSE.
                         expr2=expr2.OR.expr3
                     ENDIF
                 END DO
             ELSE
                 CYCLE
             ENDIF

             !##############################!
             !if (false) i.e. not significant
             !##############################! 
             IF ((neigh(I,K).GE.3).AND.(expr2.EQV.ans2)) THEN
                 !WRITE(*,*) atom_ele(I),atom_ele(K),'F'
                 !if nearest neigbour is 4 and I,K interactions 
                 !does not contribute significantly then do:                 

                 !1. calculate the repulsion for interactions not being
                 !fitted
                 E_POTR=0
                 E_POTR=(A(I,K)*EXP(-B(I,K)*(R(J,I,K))))
                 E_GAVR_ND(J)=E_GAVR_ND(J)+E_POTR
 
                 !2. calculate the dispersion for interactions not being
                 !fitted
                 E_POTD=0
                 E_POTD=-(C(I,K)/(R(J,I,K)**6))
                 E_GAVD_ND(J)=E_GAVD_ND(J)+E_POTD
 
                 !3. calculate the sum of repulsion-dispersion not being
                 !fitted
                 E_GAV_ND(J)=E_GAV_ND(J)+E_POTR+E_POTD
                 
                 !#########################################!
                 !if (true) i.e. interaction is significant 
                 !#########################################!
             ELSE IF ((neigh(I,K).GE.3).AND.(expr2.EQV.ans1)) THEN
                 !WRITE(*,*) atom_ele(I),atom_ele(K),'F'
                 !if nearest neigbour is 4 and interaction I,K
                 !contributes significantly, then do the following:

                 !1. calculate the repulsion for interactions to be fitted
                 E_POTR2=(A(I,K)*EXP(-B(I,K)*(R(J,I,K))))
                 !Slipt into atomtype-atomtype contribution!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 E_GAVR2(J)=E_GAVR2(J)+E_POTR2
                 !!!!!!
                 E_GAVR_at_Type(J,at_at)=E_GAVR_at_Type(J,at_at)+E_POTR2
   
                 !2. calculate the dispersion for interactions not being
                 !fitted
                 E_POTD2=-(C(I,K)/(R(J,I,K)**6))
                 E_GAVD2(J)=E_GAVD2(J)+E_POTD2

                 !store repulsion to be fitted in easy to print format
                 !E_GAVR3(J,I,K)=E_POTR2
             ELSE
                 CYCLE
             END IF
         END DO
     END DO
 END DO

 WRITE(*,3003) 'REPUL_ND','DISP_ND','REPUL_D','DISP_D'
 DO J=1,Imax
     WRITE(*,3004) E_GAVR_ND(J),E_GAVD_ND(J),E_GAVR2(J),E_GAVD2(J)
 END DO
 
 PRINT*,''
 PRINT*,'####################################################################'
 PRINT*,'Repulsion energy contribution group atomType---atomType interactions'
 PRINT*,'####################################################################'
 PRINT*,''
 PRINT 3011,(atom_D_ele1(N),'----',atom_D_ele2(N), N=1,Mmax)&
            &,'','EAb-','EGRes'

 DO J=1,Imax
     PRINT 3007,(E_GAVR_at_Type(J,at_at),at_at=1,Mmax),&
     &((E_J(J))-(E_GAVR_ND(J)+E_GAVD_ND(J)+E_GAVD2(J)))
     !&E_GAVR_ND(J)+E_GAVD_ND(J)+E_GAVD2(J),E_J(J)
 END DO

!PRINT*, ((E_J(1:Imax))-(E_GAVR_ND(1:Imax)+E_GAVD_ND(1:Imax)+E_GAVD2(1:Imax)))&
!&/(E_GAVR_at_Type(1:Imax,1:Mmax))

 CALL SYSTEM ('rm -f  summary_log_input.d')
 !looping from 1 to maximum atT--atT interaction
 DO N=1,Mmax
     !combination code
     OPEN(211,FILE='combi.f90')
     WRITE(211,*) 'program combinations'
     WRITE(211,*) '  implicit none'
     WRITE(211,*) '  integer, parameter :: m_max =', N
     WRITE(211,*) '  integer, parameter :: n_max =', Mmax
     WRITE(211,*) '    integer, dimension (m_max) :: comb'
     WRITE(211,*) "    character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', m_max - 1) // ')'"
     WRITE(211,*) '   '
     WRITE(211,*) '    call gen (1)'
     WRITE(211,*) '   '
     WRITE(211,*) '  contains'
     WRITE(211,*) '   '
     WRITE(211,*) '    recursive subroutine gen (m)'
     WRITE(211,*) '   '
     WRITE(211,*) '      implicit none'
     WRITE(211,*) '      integer, intent (in) :: m'
     WRITE(211,*) '      integer :: n'
     WRITE(211,*) '   '
     WRITE(211,*) '      if (m > m_max) then'
     WRITE(211,*) '        write (*, fmt) comb'
     WRITE(211,*) '      else'
     WRITE(211,*) '        do n = 1, n_max'
     WRITE(211,*) '          if ((m == 1) .or. (n > comb (m - 1))) then'
     WRITE(211,*) '            comb (m) = n'
     WRITE(211,*) '            call gen (m + 1)'
     WRITE(211,*) '          end if'
     WRITE(211,*) '        end do'
     WRITE(211,*) '      end if'
     WRITE(211,*) '   '
     WRITE(211,*) '    end subroutine gen'
     WRITE(211,*) '   '
     WRITE(211,*) '  end program combinations'
     
     !compile combination code

     PRINT*,''
     !PRINT 3008,(atom_D_ele1(L),'-----',atom_D_ele2(L),L=1,Mmax)
     CALL SYSTEM('gfortran -m64 combi.f90 -o comb.exe')
     REWIND(211)
     CALL SYSTEM('./comb.exe > comb.out')
     !CALL SYSTEM('./comb.exe')
     PRINT*,''
     
     !number of regressions to perform
     !Reads the number of line(s) output by combination code
     CALL SYSTEM ('wc -l < comb.out > temp_n')
     OPEN (200,FILE='temp_n')
     READ (200,*) NcombL
     CLOSE(200)

     !OPEN (212,FILE='g02ddfe.d')
      

     !number of atomT-atomT interaction in each row (combination)=N
     !regression
     DO I=1,NcombL
         !For each row
         !modify regression code
         !perform regression
         
         !store the interactions that will undergo regression
         OPEN(213,FILE='comb.out')
         READ(213,*) (c_int(K),K=1,N)
         
         !WRITE(*,*) '###########################'
         !WRITE(*,*) 'REGRESSION ANALYSIS: RESULT'
         !WRITE(*,*) '---------------------------'
         WRITE(*,*) '##############################################################'
         !WRITE(*,*) (c_int(K),K=1,N)
         WRITE(*,*) 'The following atom-type...atom-type interaction(s) were&
         & rescale'
         WRITE(*, 3008) (atom_D_ele1(c_int(K)),'-',atom_D_ele2(c_int(K)),',',K=1,N)

        ! DO J=1,Imax
        !PRINT 3007,(E_GAVR_at_Type(J,at_at),at_at=1,Mmax),&
        !&((E_J(J))-(E_GAVR_ND(J)+E_GAVD_ND(J)+E_GAVD2(J)))
        !END DO


         !1. create input
         !!!!!!!!!!!!!!!!!!!!!!!
         !FORMAT !!!!!!!!!!!!!!!
         !1 REGRESSORS Y!
         !G02DDF Example Program Data
         !15 10 'U'
         !1  102.7592468     7.5416231    18.0316296    21.2146263    11.4237890    19.1057587    10.4240065     6.2759800     4.4999952   149.1976013
         !1  103.2297897     7.5259795    17.8553143    19.3802319    11.0623226    19.1911964     9.5916662     5.0133109     4.6530294   145.0759888
         OPEN (212,FILE='g02ddfe.d')
         WRITE(212,*)   'INPUT data for multiple regression calculation'
         WRITE(212,*)   (c_int(K),K=1,N)
         WRITE(212,*)   Mmax
         WRITE(212,3009) Imax,N+1,"'U'"
         DO L=1,Imax
                          !1.00 for the intercept
                          !followed by the regressor (repulsion contribution in
                          !selected subset)
             WRITE(212,3010) '1.00',(E_GAVR_at_Type(L,c_int(K)),K=1,N),             & 
                          !Ab initio minus All other contribution not fitted 
                          &((E_J(L))-(E_GAVR_ND(L)+E_GAVD_ND(L)+E_GAVD2(L)))-       &
                          &((SUM(E_GAVR_at_Type(L,1:Mmax)))-                        &
                          &(SUM(E_GAVR_at_Type(L,c_int(1:N))))),                    & 
                          !Contributions not fitted
                          &((E_GAVR_ND(L)+E_GAVD_ND(L)+E_GAVD2(L))+                 &
                          &((SUM(E_GAVR_at_Type(L,1:Mmax)))-                        &
                          &(SUM(E_GAVR_at_Type(L,c_int(1:N))))) ),                  &
                          !Torsion
                          &deg(L),                                                  &
                          !Ab intio energy
                          &E_J(L) 
            ! WRITE(*,3010) '1.00',(E_GAVR_at_Type(L,c_int(K)),K=1,N),             &
            !              &((E_J(L))-(E_GAVR_ND(L)+E_GAVD_ND(L)+E_GAVD2(L)))-       &
            !              &((SUM(E_GAVR_at_Type(L,1:Mmax)))-                        &
            !              &(SUM(E_GAVR_at_Type(L,c_int(1:N)))))
         
         END DO
         ! WRITE(*,*) "--------------------"
         CALL SYSTEM ('cat g02ddfe.d >> summary_log_input.d')

         !2.modify the regression code
         OPEN(214,FILE='g02ddfe.f')

         WRITE(214,"(a)" ) "*     G02DDF Example Program Text"
         WRITE(214,"(a)" ) "*     Mark 14 Release. NAG Copyright 1989."
         WRITE(214,"(a)" ) "*     .. Parameters .."
         WRITE(214,"(a)" ) "      INTEGER          MMAX, NMAX, LDQ, PMAX"
         WRITE(214,' (A29,I3,A6,I3,A1)') "      PARAMETER        (MMAX=",N+2,",NMAX=",Imax,")"
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
         WRITE(214,"(a)" ) "     +                 deg(LDQ,MMAX), "
         WRITE(214,"(A27,I3,A1)" ) "     +                  BP(",MMAX,")"
         WRITE(214,"(a)" ) "*     .. External Subroutines .."
         WRITE(214,"(a)" ) "      EXTERNAL         G02DDF, G02DEF"
         WRITE(214,"(a)" ) "*     .. Executable Statements .."
         WRITE(214,"(a)" ) "*      WRITE (NOUT,*) 'G02DDF Example Program Results'"
         WRITE(214,"(a)" ) "*     Skip heading in data file"
         WRITE(214,"(a)" ) "      READ (NIN,*)"
         WRITE(214,"(A35,I3,A3)" ) "      READ (NIN,*) (Param(K,1),K=1,",N,")"
         WRITE(214,"(A36,I3,A3)" ) "      WRITE(NOUT,*) (Param(K,1),K=1,",N,")"
         WRITE(214,"(a)" ) "      WRITE(NOUT,*) 'SUMMARY FOR ABOVE'"
         WRITE(214,"(a)" ) "      READ (NIN,*) PMAX"
         WRITE(214,"(a)" ) "      READ (NIN,*) N, M, WEIGHT"
         WRITE(214,"(a)" ) "      IF (N.LE.NMAX .AND. M.LT.MMAX) THEN"
         WRITE(214,"(a)" ) "         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN"
         WRITE(214,"(a)" ) "            DO 20 I = 1, N"
         WRITE(214,"(a)" ) "               READ (NIN,*) (X(I,J),J=1,M), Q(I,1), WT(I)"
         WRITE(214,"(a)" ) "   20       CONTINUE"
         WRITE(214,"(a)" ) "         ELSE"
         WRITE(214,"(a)" ) "            DO 40 I = 1, N"
         WRITE(214,"(a)" ) "               READ (NIN,*) (X(I,J),J=1,M), Q(I,1),QA(I,1),deg(I,1)"
         WRITE(214,"(a)" ) "     +,E_J(I,1) "
         WRITE(214,"(a)" ) "   40       CONTINUE"
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
         WRITE(214,"(a)" ) "            EMjODiEL(J,1)=QA(J,1)+mats"
         WRITE(214,"(a)" ) "         END DO"
         !Relative intramolecular energy = eJodi
         WRITE(214,"(a)" ) "         DO J=1,N"
         WRITE(214,"(a)" ) "            eJODI(J,1)=EMjODiEL(J,1)-MINVAL(EMjODiEL(1:N,1))"
         !Absolute intramolecular energy
         WRITE(214,"(a)" ) "            WRITE(*,'(F12.2,A1,F12.2)') eJODI(J,1),',',"
         !Testing 
         WRITE(214,"(a)" ) "     +EMjODiEL(J,1)-B(1)"
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
         WRITE(214,"(a)" ) "            OPEN (215,FILE='Rpara.txt')"
         WRITE(214,"(a)" ) "            OPEN (216,FILE='Table.txt')"
         WRITE(214,"(a)" ) "            DO 80 J = 1, IP"
         WRITE(214,"(a)" ) "               WRITE (NOUT,99997) J-1,',', B(J),',', SE(J)"
         WRITE(214,"(a)" ) "               WRITE (215,'(1X,F20.4)',ADVANCE='NO') B(J)"
!         WRITE(214,"(a)" ) "                 DO K=1,N-1"!I=1,PMAX"
!         WRITE(214,"(a)" ) "                    DO I=1,PMAX"!K=1,N-1"
!         WRITE(214,"(a)" ) "                        IF (I.EQ.Param(K,1)) THEN"
!         WRITE(214,"(a)" ) "                            BP(I)=B(J)"
!         WRITE(214,"(a)" ) "                        ELSE"
!         WRITE(214,"(a)" ) "                            BP(I)=0.0D0"
!         WRITE(214,"(a)" ) "                        END IF"
!         WRITE(214,"(a)" ) "                    END DO"
!         WRITE(214,"(a)" ) "                 END DO"
         WRITE(214,"(a)" ) "   80       CONTINUE"
         WRITE(214,"(a)" ) "            ELSEIF (J.EQ.1) THEN"
         WRITE(214,"(a)" ) "               WRITE (NOUT,99995) 'c',',',B(J),',', SE(J)"
         WRITE(214,"(a)" ) "               WRITE (215,'(1X,F20.4)',ADVANCE='NO') B(J)"
         WRITE(214,"(a)" ) "            ELSE"
         WRITE(214,"(a)" ) "            WRITE (NOUT,*)"
         WRITE(214,"(a)" ) "            WRITE (NOUT,99996) ' ** G02DDF returned with IFAIL = ',"
         WRITE(214,"(a)" ) "     +        IFAIL"
         WRITE(214,"(a)" ) "         END IF"
         WRITE(214,"(a)" ) "      END IF"
         WRITE(214,"(a)" ) "  100 CONTINUE"


!         WRITE(214,"(a)" ) "               WRITE (216,'(1X,5F20.4)')"
!         !#standard error fit          (1)
!         WRITE(214,"(a)" ) "     +(RSS/(IDF))**0.5,"
!         !#minimum energy conformation (2)
!         WRITE(214,"(a)" ) "     +deg(MINLOC(eJODI(1:N,1)),1)," !#2
!         !#baseline energy             (3)
!         WRITE(214,"(a)" ) "     +((MINVAL(EMjODiEL(1:N,1)))-B(1))," !#3
!         !#model's relative energy @ global minimum (4)
!         WRITE(214,"(a)" ) "     +eJODI((MINLOC(E_J(1:N,1))),1)," !#4
!         !#standard error intra        (5)
!         WRITE(214,"(a)" ) "     +(mats/(IDF))**0.5" !#5
        


         !WRITE FORMATED PARAMETER IN A SINGLE LINE
         WRITE(214,"(a)" ) "             WRITE (NOUT,*) '#SUMMARY#'"
         WRITE(214,"(a)" ) "             DO 82 K=1,IP-1"
         WRITE(214,"(a)" ) "                    BP(Param(K,1))=B(K+1)"
!         WRITE(214,"(a)" ) "              WRITE (NOUT,'(1X,F10.4)',ADVANCE='NO') BP(K)"
         WRITE(214,"(a)" ) "   82        CONTINUE"
         
         WRITE(214,"(a)" ) "             DO 81 I=1,PMAX"
         WRITE(214,"(a)" ) "                    IF (ANY(I.EQ.Param).EQV..TRUE.) THEN"
         WRITE(214,"(a)" ) "                        WRITE (216,'(1X,F10.4)',ADVANCE='NO') BP(I)"
         WRITE(214,"(a)" ) "                    ELSE"
         WRITE(214,"(a)" ) "                        WRITE (216,'(1X,A4)',ADVANCE='NO') '-'"
         WRITE(214,"(a)" ) "                    END IF"     
         WRITE(214,"(a)" ) "   81        CONTINUE"
         
         WRITE(214,"(a)" ) "             WRITE (216,'(1X,F10.2)',ADVANCE='NO') B(1)"


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
      REWIND(213)
 END DO
  !clean up
   CALL SYSTEM('rm -f comb.exe combi.f90 comb.out conn1.r conn.r conn_atm  conn_table')
   CALL SYSTEM('rm -f dominant_at-at.txt dominant.txt') 
   CALL SYSTEM('rm -f g02ddfe.exe  summary_log_input.d  Table.txt')
   CALL SYSTEM('rm -f g02ddfe.d g02ddfe.f temp.txt')
   CALL SYSTEM('rm -f temp_* Rpara.txt h03adfe.d h03adfe.r')
   CALL SYSTEM ('rm -f zmat*mol2 gavezzotti_para.d h03adfe.exe summary.r') 
   CALL SYSTEM('echo "---------------------------------------"')
   !end message
   CALL SYSTEM('echo "Fitting routine terminated successfully"')
   CALL SYSTEM('echo "---------------------------------------"')
   CALL SYSTEM('echo ""')

END PROGRAM read_mol2
