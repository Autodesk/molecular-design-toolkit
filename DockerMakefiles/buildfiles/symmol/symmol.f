      program SYMMOL
      implicit double precision (a-h,o-z)
c     logical esiste
      PARAMETER (NMA=1000)
      CHARACTER*2 ATOM(NMA)
      CHARACTER*80 RIGA
c     CHARACTER*60 inpfile
      common/matcel/PC(7),PCR(7),O(3,3),OI(3,3),G(3,3),GI(3,3),CS(12)
      COMMON/AT0/NA,NMO
      COMMON/AT1/X(3,NMA),AMAS(NMA),RCOV(NMA),MSP(NMA)
      COMMON/AT2/SX(3,NMA),SIG(NMA),DXM(NMA),DCM,DCME,indwgh,indtol
      COMMON/AT3/MOL(NMA),MLG(NMA)
      CHARACTER*6 NAME,SPEC
      COMMON/AT4/NAME(NMA)
      integer out
      COMMON/output/out
c________________________________________________________________________
c if you want address all the output on the standard output, suppress the two
c sequent instruction and activate the third
c
      out=9
      OPEN(9,file='symmol.out')
c     out=6
c________________________________________________________________________
      if(out.eq.9)then
        write(*,*)'                  ============'
        write(*,*)'                     SYMMOL'
        write(*,*)' A PROGRAM FOR THE SYMMETRIZATION OF GROUPS OF ATOMS'
        write(*,*)'       By Tullio Pilati and Alessandra Forni'
        write(*,*)'               Version November 4th 2002'
        write(*,*)' ==================================================='
        write(*,*)
      endif
      write(out,*)'                  ============'
      write(out,*)'                     SYMMOL'
      write(out,*)' A PROGRAM FOR THE SYMMETRIZATION OF GROUPS OF ATOMS'
      write(out,*)'       By Tullio Pilati and Alessandra Forni'
      write(out,*)'               Version November 4th 2002'
      write(out,*)' ==================================================='
      write(out,*)
c 900 write(*,*)'input file name?'
c     read(*,'(a)')inpfile
c     inquire(file=inpfile,exist=esiste)
c     if(.not.esiste)then
c        write(*,'(a60,a)')inpfile,' does not exist'
c        stop
c     else
c        open(1,file=inpfile)
c     endif
 1000 read(*,'(A80)')RIGA
      if(RIGA(1:1).eq.'#') go to 1000
      READ(RIGA,*)(PC(i),i=1,6)
c********************************************
 1010 read(*,'(A80)')RIGA
      if(RIGA(1:1).eq.'#') go to 1010
      READ(RIGA,*)indwgh,indtol,DCM,DCME
      if(DCME.lt.DCM)DCME=DCM
      write(out,*)
      if(indwgh.le.1)then
      indwgh=1
      write(out,*)' INDWGH=1 ===> WEIGHTS AS ATOMIC MASS'
      endif
      if(indwgh.eq.2)write(out,*)' INDWGH=2 ===> UNITARY WEIGHTS'
      if(indwgh.eq.3)write(out,*)' INDWGH=3 ===> WEIGHTS BASED ON S.U.'
      if(indtol.le.1)then
      indtol=1
      write(out,*)' INDTOL=1 ===> TOLERANCE=CONSTANT'
      endif
      if(indtol.eq.2)write(out,*)
     *               ' INDTOL=2 ===> TOLERANCE BASED ON DISTANCES'
      if(indtol.eq.3)
     *         write(out,*)' INDTOL=3 ===> TOLERANCE BASED ON S.U.'
      write(out,1)DCM,DCME
    1 format('  CONSTANTS OF TOLERANCE=',2f7.3)
      write(out,*)
c********************************************
      i=1
      MOmax=0
 1100 read(*,'(A80)',END=1200)RIGA
      if(i.gt.NMA)then
        write(out,*)'ERROR: TOO MANY ATOMS IN INPUT'
        write(out,*)'Enlarge PARAMETER NMA'
        write(out,*)'Actual value: PARAMETER (NMA=1000)'
        stop
      endif
c     write(out,*)riga
      if(RIGA(1:1).eq.'#')go to 1100
      READ(RIGA,'(A6,I2,6f9.5)')NAME(i),MOL(i),(X(k,i),k=1,3),
     *(SX(k,i),k=1,3)
      if(MOL(i).gt.MOmax)Momax=MOL(i)
      ATOM(i)='  '
      SPEC=NAME(i)
      call COMPATTA(SPEC,6,K)
      K=1
      if(SPEC(2:2).ge.'A'.and.SPEC(2:2).le.'Z')K=2
      if(SPEC(2:2).ge.'a'.and.SPEC(2:2).le.'z')K=2
      call MAIUSCOL(SPEC(1:1))
      IF(K .GT. 1) then
      call MAIUSCOL(SPEC(2:2))
      ATOM(i)=SPEC(1:2)
      else
      ATOM(i)(2:2)=SPEC(1:1)
      endif
      i=i+1
      go to 1100
 1200 NA=i-1
      if(indwgh.eq.3.or.indtol.eq.3) then
        sxmax=0.000001d0
        do i=1,NA
            do k=1,3
              if(sx(k,i).gt.sxmax) sxmax=sx(k,i)
            enddo
        enddo
        if(sxmax.lt.0.00001d0) then
          write(out,*)
          write(out,*) ' ERROR: S.U. ON COORDINATES ARE MISSING.'
          if(indtol.eq.3)
     *       write(out,*) '       INTRODUCE THE S.U. OR CHANGE INDTOL'
          if(indwgh.eq.3)
     *       write(out,*) '       INTRODUCE THE S.U. OR CHANGE INDWGH'
          stop
        endif
      endif
C
      call mass(NA,ATOM)
c     call costanti
      call cella
      write(out,*)
      write(out,*)' ATOM GROUP  INPUT COORDINATES AND THEIR S.U.'
      write(out,*)
      do I=1,NA
      write(out,'(1x,A6,i2,1x,6f9.5,2x,i2)')NAME(I),MOL(I)
     *,(X(k,I),k=1,3),(SX(k,I),k=1,3)
      enddo
C
      do 1300 MO=1,MOmax
      call work(MO)
 1300 continue
      end
      subroutine arccos(cost,ang)
C arco coseno in gradi
      implicit double precision (a-h,o-z)
      RAD=57.29577951308232D0
      if(cost.gt.1.d0)cost=1.d0
      if(cost.lt.-1.d0)cost=-1.d0
      ang=RAD*DACOS(cost)
      return
      end
      subroutine asymunit(MK,IASU,N,NMS)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      DIMENSION MK(NMA,NMG),IASU(NMA)
      do i =1,N
        IASU(I)=2
      enddo
      DO 2800 I=1,N
        if(IASU(I).ne.2)go to 2800
      DO 2790 J=1,NMS
      K=IABS(MK(I,J))
      if(K.eq.I)go to 2790
      if(K.eq.0)then
        IASU(I)=0
        go to 2800
      endif
        IASU(K)=1
 2790 CONTINUE
 2800 CONTINUE
      return
      end
      subroutine ax_order(A,i,m,msign,invers)
C
C m = group order for the matrix SIM(i)
C msign= 1 asse di rotazione propria
C msign=-1 asse di rotazione impropria
C     
      PARAMETER (maxorder=8)
      implicit double precision (a-h,o-z)
      dimension A(3,3,*),B(3,3),C(3,3,2)
      integer out
      COMMON/output/out
      call azzera(C,0.d0,18)
      C(1,1,1)=1.d0
      C(2,2,1)=1.d0
      C(3,3,1)=1.d0
      C(1,1,2)=-1.d0
      C(2,2,2)=-1.d0
      C(3,3,2)=-1.d0
      invers=2
      msign=NINT(det(A,i))
        call prodmm(A,C,B,i,1,1)
      do m=1,2*maxorder
        if(ium(C,B,1.d-2,2,1).eq.1)invers=1
        if(ium(C,B,1.d-2,1,1).eq.1)go to 1000
        call prodmm(A,B,B,i,1,1)
      enddo
      write(out,*)'INPUT PARAMETER DCM probably too HIGH. Reduce it!'
      stop
1000  if(m.lt.6.or.msign.eq.1)return
      m1=(m/4)*4
      if(m1.ne.m)m=m/invers
c     write(out,*)'ax_order,i,m,msign,invers'
c     write(out,'(10i5)')i,m,msign,invers
c     write(out,'(3f10.6)')((A(ii,kk,i),kk=1,3),ii=1,3)
      return
      end
      subroutine azzera(A,COST,N)
      implicit double precision (a-h,o-z)
C azzeraMENTO
      DIMENSION A(1)
      DO 1 I=1,N
    1 A(I)=COST
      RETURN
      END
      subroutine bond(XO,IASU,N)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      integer out
      COMMON/output/out
      COMMON/AT1/X(3,NMA),AMAS(NMA),RCOV(NMA),MSP(NMA)
      CHARACTER*6 NAME
      COMMON/AT4/NAME(NMA)
      COMMON/AT5/LEGAMI(20,NMA),NLEG(NMA)
      DIMENSION XO(3,NMA),IASU(NMA)
      DIMENSION T(3)
      do 1200 i=1,N
      r1=RCOV(i)
      N1=0
      do 1100 k=1,N
      if(IASU(k).eq.0)go to 1100
      if(k.eq.i)go to 1100
      r2=(r1+RCOV(k))*1.2
      call combli(XO,XO,T,1.d0,-1.d0,i,k,1)
      call prods(t,t,d1,1,1,2)
c     write(out,'(a20,2i3,2f10.5)')'i,k,d1,r2',i,k,d1,r2
      if(d1.gt.r2)go to 1100
      N1=N1+1
      if(N1.gt.20)then
        go to 1200
      endif
      LEGAMI(N1,I)=k
 1100 continue
      NLEG(i)=N1
 1200 continue
c     write(out,*)'LEGAMI'
c     do ii=1,N
c     N1=NLEG(II)
c     write(out,'(20i3)')(LEGAMI(kk,ii),kk=1,N1)
c     enddo
      return
      end
      double precision function brent(ax,bx,cx,func,tol,xmin)
      implicit double precision (a-h,o-z)
      INTEGER ITMAX
c     REAL brent,ax,bx,cx,tol,xmin,func,CGOLD,ZEPS
      EXTERNAL func
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
c     REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=func(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=func(u)
        write(6,*)'brent,fu,fx',fu,fx
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
      subroutine cella
      implicit double precision (a-h,o-z)
      integer out
      COMMON/output/out
      common/matcel/PC(7),PCR(7),O(3,3),OI(3,3),G(3,3),GI(3,3),CS(12)
      RAD=57.29577951308232D0
      ARAD=1.D0/RAD
      PC(7)=1.D0
      PCR(7)=1.D0
      COM=1.D0
      DO 1030 I=1,3
      K=I+3
      CS(I)=COS(PC(K)*ARAD)
      CS(K)=SIN(PC(K)*ARAD)
      IF(PC(K).GT.0.D0)GO TO 1030
      PC(K)=90.D0
      CS(I)=0.D0
      CS(K)=1.D0
 1030 COM=COM-CS(I)*CS(I)
      COM=SQRT(COM+2.D0*CS(1)*CS(2)*CS(3))
      O(1,1)=PC(1)*CS(6)
      O(1,2)=0.D0
      O(1,3)=PC(3)*(CS(2)-CS(1)*CS(3))/CS(6)
      O(2,1)=PC(1)*CS(3)
      O(2,2)=PC(2)
      O(2,3)=PC(3)*CS(1)
      O(3,1)=0.D0
      O(3,2)=0.D0
      O(3,3)=PC(3)*COM/CS(6)
      IF(ABS(O(2,1)).LT.1.D-10)O(2,1)=0.D0
      IF(ABS(O(2,3)).LT.1.D-10)O(2,3)=0.D0
      IF(ABS(O(1,3)).LT.1.D-10)O(1,3)=0.D0
      PC(7)=det(O,1)
      call traspo(O,G,1,1)
      call prodmm(G,O,G,1,1,1)
      call inv33(O,OI,1,1)
      PCR(7)=det(OI,1)
      call inv33(G,GI,1,1)
      PCR(1)=SQRT(GI(1,1))
      PCR(2)=SQRT(GI(2,2))
      PCR(3)=SQRT(GI(3,3))
      CS(7)=GI(2,3)/(PCR(2)*PCR(3))
      CS(8)=GI(1,3)/(PCR(1)*PCR(3))
      CS(9)=GI(1,2)/(PCR(1)*PCR(2))
      PCR(4)=RAD*DACOS(CS(7))
      PCR(5)=RAD*DACOS(CS(8))
      PCR(6)=RAD*DACOS(CS(9))
      CS(10)=SIN(ARAD*PCR(4))
      CS(11)=SIN(ARAD*PCR(5))
      CS(12)=SIN(ARAD*PCR(6))
      write(out,*)
      write(out,*)' CELL'
      write(out,'(3f10.5,3f10.3,f11.5,/)')PC
c     write(out,*)' RECIPROCAL CELL'
c     write(out,'(3f10.5,3f10.3,f11.5,/)')PCR
c     write(out,*)' METRIC MATRIX'
c     write(out,'(3g16.6)')((G(i,k),k=1,3),i=1,3)
c     write(out,*)
c     write(out,*)' ORTHOGONALIZATION MATRIX'
c     write(out,'(3g16.6)')((O(i,k),k=1,3),i=1,3)
c     write(out,*)
      return
      END
      subroutine combli(A,B,C,D,E,I,J,K)
COMBINAZIONE LINEARE D*A+B*E=C
      implicit double precision (a-h,o-z)
      DIMENSION A(3,1),B(3,1),C(3,1)
      C(1,K)=A(1,I)*D+B(1,J)*E
      C(2,K)=A(2,I)*D+B(2,J)*E
      C(3,K)=A(3,I)*D+B(3,J)*E
      RETURN
      END
      subroutine compatta(WORD,L,K)
      implicit double precision (a-h,o-z)
      CHARACTER*(*) WORD
      K=0
      DO 1 I=1,L
      IF(WORD(I:I).LE.' '.OR.WORD(I:I).GT.'~')GO TO 1
      K=K+1
      WORD(K:K)=WORD(I:I)
    1 CONTINUE
      N=K+1
      DO 2 I=N,L
    2 WORD(I:I)=' '
      RETURN
      END
      subroutine compl_gr(MK,N,*)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      DIMENSION CO(3,3),MK(NMA,NMG)
complete the group
 2520 NN=NMS
      DO 2600 I=1,NN
      DO 2600 J=1,NN
      call prodmm(SIM,SIM,CO,I,J,1)
      DO 2560 JJ=1,NN
      L=ium(CO,SIM,1.d-2,1,JJ)
      IF(L.EQ.1)GO TO 2590
 2560 CONTINUE
      GO TO 2610
C group multiplication table
 2590 MTG(I,J)=JJ
 2600 CONTINUE
c Gruppo completato
      RETURN 
C HA TROVATO UNA NUOVA MATRICE
 2610 NMS=NMS+1
      IF(NMS.LE.NMG) GO TO 2630
      write(6,2)
    2 format(' ERROR: TOO MANY MATRICES FOUND')
      write(6,'(3(3F10.5,/),/)')(((SIM(I,J,K),J=1,3),I=1,3),K=1,NMS)
c ritorno per errore
      RETURN 1
 2630 call trasfm(CO,SIM,1,NMS)
      do 2640 k=1,N
      k1=MK(k,J)
      MK(k,NMS)=MK(k1,I)
 2640 continue
      go to 2520
      end
      double precision function crms(t)
      implicit double precision (a-h,o-z)
C
C ppu=vettore XO ripetuto NMS volte
C ppo=vettore XS trasformato da MK(NMA,NMS)
C t(3)=vettore variazione angolare in radianti
C
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      PARAMETER (NMV=NMA*NMG)
      common/rmsmin/RV(3,3),ppo(3,NMV),ppu(3,NMV),npu
      dimension t(3),po(3),pu(3),rpu(3)
      dimension RX(3,3),RY(3,3),RZ(3,3),RR(3,3)
c     write(6,'(6F9.5)')t
      call azzera(RX,0.d0,9)
      call azzera(RY,0.d0,9)
      call azzera(RZ,0.d0,9)
      call azzera(RR,0.d0,9)
      RX(1,1)=1.d0
      RX(2,2)=dcos(t(1))
      RX(2,3)=-dsin(t(1))
      RX(3,2)=-RX(2,3)
      RX(3,3)=dcos(t(1))
      RY(2,2)=1.d0
      RY(1,1)=dcos(t(2))
      RY(1,3)=-dsin(t(2))
      RY(3,1)=-RY(1,3)
      RY(3,3)=dcos(t(2))
      RZ(1,1)=dcos(t(3))
      RZ(1,2)=-dsin(t(3))
      RZ(2,1)=-RZ(1,2)
      RZ(2,2)=dcos(t(3))
      RZ(3,3)=1.d0
c     write(6,'(3g12.5)')RX
c     write(6,*)
c     write(6,'(3g12.5)')RY
c     write(6,*)
c     write(6,'(3g12.5)')RZ
c     write(6,*)
      call prodmm(RY,RX,RR,1,1,1)
      call prodmm(RZ,RR,RV,1,1,1)
C ortonormalizzazione di precisione
      call prodv(RV,RV,RV,1,2,3)
      call prodv(RV,RV,RV,3,1,2)
      call prodv(RV,RV,RV,2,3,1)
      call norm(RV,1)
      call norm(RV,2)
      call norm(RV,3)
c     write(6,'(3g12.5)')RV
C fine ortonormalizzazione di precisione
      func=0.d0
      do i=1,npu
       po(1)=ppo(1,i)
       po(2)=ppo(2,i)
       po(3)=ppo(3,i)
       pu(1)=ppu(1,i)
       pu(2)=ppu(2,i)
       pu(3)=ppu(3,i)
       call prodmv(RV,pu,pu,1,1,1)
c     write(6,'(6F9.5)')po,pu
       func=func+(po(1)-pu(1))**2+(po(2)-pu(2))**2+(po(3)-pu(3))**2
      enddo
      crms=func
c     write(6,*)func
      return
      end
      FUNCTION det(X,N)
      implicit double precision (a-h,o-z)
C DETERMINANTE DELLA MATRICE X
      DIMENSION X(3,3,1)
      DET=  +X(1,1,N)*X(2,2,N)*X(3,3,N)-X(1,1,N)*X(2,3,N)*X(3,2,N)
      DET=DET+X(1,2,N)*X(2,3,N)*X(3,1,N)-X(1,2,N)*X(2,1,N)*X(3,3,N)
      DET=DET+X(1,3,N)*X(2,1,N)*X(3,2,N)-X(1,3,N)*X(2,2,N)*X(3,1,N)
      RETURN
      END
      subroutine eigen(A,VEC,EIG,W,GAM,BET,BSQ,P,Q,IPO,IORD,IVP,NN)
      implicit double precision (a-h,o-z)
C
C     --------------
C     QCPE VERSION
C     DECEMBER 1971
C     --------------
C
C     MATRIX DIAGNOLIZATION ROUTINE FOR REAL SYMMETRIC CASE
C     HOUSEHOLDER METHOD
C     RHO=UPPER LIMIT FOR OFF-DIAGONAL ELEMENT
C     NN SIZE OF MATRIX
C     A=MATRIX (ONLY LOWER TRIANGLE IS USED+THIS IS DESTROYED
C     EIG=RETURNED eigenVALUES IN ALGEBRAIC DESCENDING ORDER
C     VEC=RETURNED eigenVECTORS IN COLUMNS
C
C
      DIMENSION A(NN,NN),VEC(NN,NN),EIG(NN),W(NN),GAM(NN),BET(NN),
     *BSQ(NN),P(NN),Q(NN),IPO(NN),IORD(NN),IVP(NN)
C
      DATA RHOSQ/1.0D-12/
      ADUE=.50D0
      ZERO=0.D0
      UNO=1.0D0
      DUE=2.0D0
C
      N=NN
      IF(N)10,550,10
   10 N1=N-1
      N2=N-2
      GAM(1)=A(1,1)
      IF(N2) 180,170,20
   20 DO 160 NR=1,N2
      B=A(NR+1,NR)
      S=ZERO
      DO 30 I=NR,N2
   30 S=S+A(I+2,NR)**2
C     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
      A(NR+1,NR)=ZERO
      IF (S) 150,150,40
   40 S=S+B*B
      SGN=UNO
      IF (B) 50,60,60
   50 SGN=-UNO
   60 SQRTS=SQRT(S)
      D=SGN/(SQRTS+SQRTS)
      TEMP=SQRT(ADUE+B*D)
      W(NR)=TEMP
      A(NR+1,NR)=TEMP
      D=D/TEMP
      B=-SGN*SQRTS
C     D IS FACTOR OF PROPORTIONALITY. NOW COMPUTE AND SAVE W VECTOR.
C     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED.
      DO 70 I=NR,N2
      TEMP=D*A(I+2,NR)
      W(I+1)=TEMP
   70 A(I+2,NR)=TEMP
C     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.
C     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)
      WTAW=ZERO
      DO 120 I=NR,N1
      I1=I+1
      SUM=ZERO
      DO 80 J=NR,I
   80 SUM=SUM+A(I1 ,J+1)*W(J)
      IF(N1-I1) 110,90,90
   90 DO 100 J=I1,N1
  100 SUM=SUM+A(J+1,I1 )*W(J)
  110 P(I)=SUM
      WWWI=W(I)
  120 WTAW=WTAW+SUM*WWWI
C     P VECTOR AND SCALAR K  NOW STORED. NEXT COMPUTE Q VECTOR
      DO 130 I=NR,N1
  130 Q(I)=P(I)-WTAW*W(I)
C     NOW FORM PAP MATRIX, REQUIRED PART
      DO 140 J=NR,N1
      QJ=Q(J)
      WJ=W(J)
      DO 140 I=J,N1
  140 A(I+1,J+1)=A(I+1,J+1)-DUE*(W(I)*QJ+WJ*Q(I))
  150 BET(NR)=B
      BSQ(NR)=B*B
  160 GAM(NR+1)=A(NR+1,NR+1)
  170 B=A(N,N1)
      BET(N1)=B
      BSQ(N1)=B*B
      GAM(N)=A(N,N)
  180 BSQ(N)=ZERO
C     ADJOIN AN IDENTIFY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS.
      DO 200 I=1,N
      DO 190 J=1,N
  190 VEC(I,J)=ZERO
  200 VEC(I,I)=UNO
      M=N
      SUM=ZERO
      NPAS=1
      GO TO 330
  210 SUM=SUM+SHIFT
      COSA=UNO
      G=GAM(1)-SHIFT
      PP=G
      PPBS=PP*PP+BSQ(1)
      PPBR=SQRT(PPBS)
      DO 300 J=1,M
      COSAP=COSA
      IF(PPBS)230,220,230
  220 SINA=ZERO
      SINA2=ZERO
      COSA=UNO
      GO TO 270
  230 SINA=BET(J)/PPBR
      SINA2=BSQ(J)/PPBS
      COSA=PP/PPBR
C     POSTMULTIPLY IDENTITY BY P-TRANSPOSE MATRIX
      NT=J+NPAS
      IF(NT-N)250,240,240
  240 NT=N
  250 DO 260 I=1,NT
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)
  260 VEC(I,J)=TEMP
  270 DIA=GAM(J+1)-SHIFT
      U=SINA2*(G+DIA)
      GAM(J)=G+U
      G=DIA-U
      PP=DIA*COSA-SINA*COSAP*BET(J)
      IF(J-M)290,280,290
  280 BET(J)=SINA*PP
      BSQ(J)=SINA2*PP*PP
      GO TO 310
  290 PPBS=PP*PP+BSQ(J+1)
      PPBR=SQRT(PPBS)
      BET(J)=SINA*PPBR
  300 BSQ(J)=SINA2*PPBS
  310 GAM(M+1)=G
C     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
      NPAS=NPAS+1
      IF(BSQ(M)-RHOSQ)320,320,350
  320 EIG(M+1)=GAM(M+1)+SUM
  330 BET(M)=ZERO
      BSQ(M)=ZERO
      M=M-1
      IF(M)340,380,340
  340 IF(BSQ(M)-RHOSQ)320,320,350
C     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE
C     AS ESTIMATE OF eigenVALUE TO USE FOR SHIFT
  350 A2=GAM(M+1)
      R2=ADUE*A2
      R1=ADUE*GAM(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP=SQRT(DIF*DIF+BSQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=ABS(A2-R1)-ABS(A2-R2)
      IF(DIF)370,360,360
  360 SHIFT=R2
      GO TO 210
  370 SHIFT=R1
      GO TO 210
  380 EIG(1)=GAM(1)+SUM
C     INITIALIZE AUXILIARY TABLES REQUIRED FOR REARRANGING THE VECTORS
      DO 390 J=1,N
      IPO(J)=J
      IVP(J)=J
  390 IORD(J)=J
C     USE A TRANSPOSITION SORT TO ORDER THE eigenVALUES
      M=N
      GO TO 430
  400 DO 420 J=1,M
      IF(EIG(J)-EIG(J+1))410,420,420
  410 TEMP=EIG(J)
      EIG(J)=EIG(J+1)
      EIG(J+1)=TEMP
      ITEMP=IORD(J)
      IORD(J)=IORD(J+1)
      IORD(J+1)=ITEMP
  420 CONTINUE
  430 M=M-1
      IF(M)400,440,400
  440 IF(N1)450,490,450
  450 DO 480 L=1,N1
      NV=IORD(L)
      NP=IPO(NV)
      IF(NP-L)460,480,460
  460 LV=IVP(L)
      IVP(NP)=LV
      IPO(LV)=NP
      DO 470 I=1,N
      TEMP=VEC(I,L)
      VEC(I,L)=VEC(I,NP)
  470 VEC(I,NP)=TEMP
  480 CONTINUE
C     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
  490 DO 540 NRR=1,N
      K=N1
  500 K=K-1
      IF(K)540,540,510
  510 SUM=ZERO
      DO 520 I=K,N1
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)
      SUM=SUM+SUM
      DO 530 I=K,N1
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)
      GO TO 500
  540 CONTINUE
  550 RETURN
      END
      subroutine eq_plane(t,u,v,a)
      implicit double precision (a-h,o-z)
      dimension t(3),u(3),v(3),AA(3,3),DD(3,3),a(4)
C   Equazione del piano passante per i punti t,u,v nella forma canonica
C   a(1).x+a(2).y+a(3).z+a(4)=0 
C   con a(1),a(2),a(3)=coseni direttori a(4)=-distanza piano-origine
C vedi International Tables for Crystallography II, p. 43, equaz. 2,3,4, 6
      DD(1,1)=t(1)
      DD(1,2)=t(2)
      DD(1,3)=t(3)
      DD(2,1)=u(1)
      DD(2,2)=u(2)
      DD(2,3)=u(3)
      DD(3,1)=v(1)
      DD(3,2)=v(2)
      DD(3,3)=v(3)
      do i=1,3
      call trasfm(DD,AA,1,1)
      call azzera(AA(1,i),1.d0,3)
      a(i)=det(AA,1)
      enddo
      call prods(a,a,dist,1,1,2)
      a(1)=a(1)/dist
      a(2)=a(2)/dist
      a(3)=a(3)/dist
      a(4)=-det(DD,1)/dist
      return
      end
      subroutine icosahed(XO,PESO,N,MN,MK,MD,II,MDEG,*)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      CHARACTER*6 NAME
      COMMON/AT4/NAME(NMA)
      common/matcel/PC(7),PCR(7),O(3,3),OI(3,3),G(3,3),GI(3,3),CS(12)
      COMMON/AT0/NA,NMO
      COMMON/AT2/SX(3,NMA),SIG(NMA),DXM(NMA),DCM,DCME,indwgh,indtol
      COMMON/ORIE/OR(3,3),OT(3,3),OTI(3,3),BARC(3),BARO(3),RIN(3)
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      DIMENSION XO(3,NMA),PESO(NMA),MD(NMA,2),MK(NMA,NMG),MN(NMA)
      DIMENSION DA(NMA),meq(NMA)
      dimension eqp(4),dp(5),vd(3,5),mp(5),io(5),A(3,3),B(3,3),V(3)
      logical ico
      RAD=57.29577951308232D0
      do 2000 I1=1,N-4
      if(MD(I1,1).ne.II)go to 2000
      mp(1)=I1
      io(1)=I1
      in1=I1+1
      do 1900 I2=in1,N-3
      if(MD(I1,1).ne.II)go to 1900
      mp(2)=I2
      in2=I2+1
      do 1800 I3=in2,N-2
      if(MD(I3,1).ne.II)go to 1800
      mp(3)=I3
      call eq_plane(XO(1,I1),XO(1,I2),XO(1,I3),eqp)
c     write(out,*)'equazione del piano'
c     write(out,'(4f10.6,3i5)')eqp,mp(1),mp(2),mp(3)
c     write(out,'(4f10.6,3i5)')DXM(I1),DXM(I2),DXM(I3),eqp(4)
C se il piano e' troppo vicino all'origine viene scartato
      if(eqp(4).lt.0.5)go to 1800
      in3=I3+1
      do 1700 I4=in3,N-1
      if(MD(I4,1).ne.II)go to 1700
      mp(4)=I4
      d2=eqp(1)*XO(1,I4)+eqp(2)*XO(2,I4)+eqp(3)*XO(3,I4)+eqp(4)
      if(DABS(d2).gt.DXM(I4))go to 1700
c     write(out,'(A20,i5,f10.5)')'IV atomo del piano',mp(4),d2
 1300 in4=I4+1
      do 1600 I5=in4,N
      if(MD(I5,1).ne.II)go to 1600
      mp(5)=I5
      d2=eqp(1)*XO(1,I5)+eqp(2)*XO(2,I5)+eqp(3)*XO(3,I5)+eqp(4)
      if(DABS(d2).gt.DXM(I5))go to 1600
c     write(out,'(A20,i5,f10.5)')' V atomo del piano',mp(5),d2
C in mp ci sono i possibili 5 atomi equivalenti rispetto all'asse 5
      dmin=100.
      dp(1)=dmin
      do i=2,5
        io(i)=0
        dp(i)=0.d0
        call combli(XO,XO,vd,1.d0,-1.d0,I1,mp(i),1)
        call prods(vd,vd,dp(i),1,1,2)
        if(dp(i).lt.dmin)then
          dmin=dp(i)
          jj=i
        endif
      enddo
      io(2)=mp(jj)
      dp(jj)=100.
      dmin=100.
      do i=2,5
        if(dp(i).lt.dmin.and.io(2).ne.mp(i))then
          dmin=dp(i)
          kk=i
        endif
      enddo 
      io(5)=mp(kk)
      do i=2,5
      if(mp(i).ne.io(2).and.mp(i).ne.io(5).and.io(3).eq.0)io(3)=mp(i)
      if(mp(i).ne.io(2).and.mp(i).ne.io(5).and.mp(i).ne.io(3))jj=mp(i)
      enddo
      call combli(XO,XO,vd,1.d0,-1.d0,io(2),io(3),1)
      call prods(vd,vd,d1,1,1,1)
      call combli(XO,XO,vd,1.d0,-1.d0,io(2),jj,1)
      call prods(vd,vd,d2,1,1,1)
      io(4)=jj
      if(d2.lt.d1)then
      io(4)=io(3)
      io(3)=jj
      endif
C ora i cinque atomi sono ordinati
C per controllo preliminare prima di passare alla verify, controllo 
che siano uguali ( a meno della tolleranza) le distanze 1-2
      dmed=0.
      do i=1,5
        k=i+1
        if(k.gt.5)k=k-5
        call combli(XO,XO,vd,1.d0,-1.d0,io(i),io(k),i)
        call prods(vd,vd,dp(i),i,i,2)
        dmed=dmed+dp(i)
      enddo
      dmed=dmed*0.2d0
c     write(out,'(a5,5i5)')'io',io
c     write(out,'(a5,5f10.5)')'dp',dp
c     write(out,'(a5,5f10.5)')'dmed',dmed
      call AZZERA(A,0.d0,9)
      do i=1,5
        if(DABS(dmed-dp(i)).gt.DXM(i)*.5)go to 1600
c asse C5 come somma dei vertici del pentagono
        A(1,3)=A(1,3)+XO(1,io(i))
        A(2,3)=A(2,3)+XO(2,io(i))
        A(3,3)=A(3,3)+XO(3,io(i))
      enddo
c     write(out,'(a5,3f10.6)')'A(3)',(A(i,3),i=1,3)
C mette il riferimento in modo che il primo asse 5 coincida con z e che
C il secondo sia nel piano xz
      A(1,1)=XO(1,io(1))+XO(1,io(2))-XO(1,io(3))-XO(1,io(4))+
     *XO(1,io(5))
      A(2,1)=XO(2,io(1))+XO(2,io(2))-XO(2,io(3))-XO(2,io(4))+
     *XO(2,io(5))
      A(3,1)=XO(3,io(1))+XO(3,io(2))-XO(3,io(3))-XO(3,io(4))+
     *XO(3,io(5))
C se la somma dei cinque atomi da' un vettore nullo, cioe' il pentagono
C e' sul cerchio massimo (caso possibile per un gruppo di 30 atomi sugli
C assi C2) passo ad altro pentagono
 1400 call norm(A,3)
      call norm(A,1)
      call prodv(A,A,A,3,1,2)
      call norm(A,2)
      call prodv(A,A,A,2,3,1)
      call traspo(A,A,1,1)
c     write(out,'(a5,3(/,3f10.6))')'A',((A(i,k),i=1,3),k=1,3)
      call r_frame(A,XO,N)
Costruzione dell'asse C5
      call azzera(B,0.d0,9)
      ROT=72.d0/RAD
      CA=COS(ROT)
      CB=SIN(ROT)
      B(1,1)=CA
      B(2,2)=CA
      B(1,2)=-CB
      B(2,1)=CB
      B(3,3)=1.d0
c     write(out,'(a5,3(/,3f10.6))')'B',((B(i,k),i=1,3),k=1,3)
c verifica e ottimizza il primo asse 5
      call verify(XO,B,MK,MN,MV,N)
c     write(out,*)'primo asse C5 MV,NMS',MV,NMS
c     write(out,'(30i3)')(MK(ll,NMS),ll=1,N)
      if(MV.eq.0)go to 1600
      call opt_axis(XO,PESO,V,MK,N,2)
c ottimizza l'sse x con la stessa tecnica precedente
      call azzera(A,0.d0,9)
      if(V(3).lt.0)then
      V(1)=-V(1)
      V(2)=-V(2)
      V(3)=-V(3)
      endif
      call trasfv(V,A,1,3)
      A(1,1)=XO(1,io(1))+XO(1,io(2))-XO(1,io(3))-XO(1,io(4))+
     *XO(1,io(5))
      A(2,1)=XO(2,io(1))+XO(2,io(2))-XO(2,io(3))-XO(2,io(4))+
     *XO(2,io(5))
      A(3,1)=XO(3,io(1))+XO(3,io(2))-XO(3,io(3))-XO(3,io(4))+
     *XO(3,io(5))
 1500 call norm(A,3)
      call norm(A,1)
      call prodv(A,A,A,3,1,2)
      call norm(A,2)
      call prodv(A,A,A,2,3,1)
      call traspo(A,A,1,1)
c     write(out,'(a5,3(/,3f10.6))')'A',((A(i,k),i=1,3),k=1,3)
      call r_frame(A,XO,N)
      go to 2010
 1600 continue
 1700 continue
 1800 continue
 1900 continue
 2000 continue
c uscita normale (non ha trovato nessun asse C5)
      if(NMS.eq.1)return
completa il gruppo C5. Indispensabile qui!!!
 2010 call compl_gr(MK,N,*2015)
C scelgo i due pentagoni piu' distanti dall'asse C5 per trovare
C  un asse C2. nel caso che gli atomi del sottoset non siano
C tutti indipendenti, postrebbero esserci piu' di cinque atomi
C sul piano parallelo al primo pentagono e quindi devo 
C testare tutti i possibili pentagoni
 2015 DM=0.d0
      do 2020 i=1,N
      if(MD(I1,1).ne.II)go to 2020
      DA(i)=XO(1,I)**2+XO(2,I)**2
      if(DM.gt.DA(i))go to 2020
      kk=i
      DM=DA(i)
 2020 continue
      io(1)=kk
      dp1=XO(3,kk)
c     write(out,*)'io(1),dp1',io(1),dp1
      do k=1,4
        k1=k+1
        io(k1)=MK(io(k),2)
        dp1=dp1+XO(3,io(k1))
      enddo
      dp1=-0.2d0*dp1
C in io(i) ci sono i cinque atomi del primo pentagono
C in meq(i) vanno gli atomi del piano parallelo
      me=0
      do  2030 i=1,N
      com=DABS(dp1-XO(3,i))
      if(com.gt.DXM(i))go to 2030
C eliminazione (eventuale) degli atomi del primo piano
      do k=1,5
        if(i.eq.io(k))go to 2030
      enddo
      me=me+1
      meq(me)=i
 2030 continue
c     write(out,*)'me,meq',me
c     write(out,'(30i4)')(meq(i),i=1,me)
C scelta degli atomi del secondo pentagono
 2040 do i=1,me
        if(meq(i).ne.0)go to 2050
      enddo
C non ci sono assi C2!!!
      MDEG=1
      return
C i e' il primo atomo del nuovo pentagono
 2050 mp(1)=meq(i)
      meq(i)=0
      k=1
C cerca gli altri atomi del secondo penatgono e annulla i relativi meq
      do k=1,4
        k1=k+1
        mp(k1)=MK(mp(k),2)
        do l=1,me
          if(mp(k1).eq.meq(l))meq(l)=0
        enddo
      enddo
C somma vettori del primo pentagono opportunamente ruotati
C costruzione dell'asse C2 perpendicolare a C5 e parallelo a x
 2210 call azzera(A,0.d0,9)
      call azzera(vd,0.d0,15)
c     write(out,*)'io,mp'
c     write(out,'(5i4)')io,mp
      do 2250 i=1,5
      k=io(i)
      do 2220 l=1,NMS
      if(MK(k,l).ne.io(1))go to 2220
      call prodmv(SIM,XO,vd,l,k,3)
      call combli(vd,vd,vd,1.d0,1.d0,3,1,1)
c     write(out,'(a20,i5,3f10.5)')'io(i),vd(1)',io(i),(vd(kk,1),kk=1,3)
      go to 2230
 2220 continue
 2230 k=mp(i)
      do 2240 l=1,NMS
      if(MK(k,l).ne.mp(1))go to 2240
      call prodmv(SIM,XO,vd,l,k,3)
      call combli(vd,vd,vd,1.d0,1.d0,3,2,2)
c     write(out,'(a20,i5,3f10.5)')'mp(i),vd(2)',mp(i),(vd(kk,2),kk=1,3)
      go to 2250
 2240 continue
 2250 continue
      vd(3,1)=0.d0
      vd(3,2)=0.d0
      call norm(vd,1)
      call norm(vd,2)
c     write(out,'(a20,3f10.5)')'vd(1)',(vd(kk,1),kk=1,3)
c     write(out,'(a20,3f10.5)')'vd(2)',(vd(kk,2),kk=1,3)
      cost=DCOS(36.1/RAD)
      do i=1,5
        call prodmv(SIM,vd,vd,2,2,2)
        call prods(vd,vd,cosa,1,2,1)
c     write(out,'(a14,8f9.5)')'vd(1),vd(2)',(vd(kk,1),kk=1,3),
c    *(vd(kk,2),kk=1,3),cosa,cost
        if(cosa.ge.cost)go to 2300
      enddo
 2300 call combli(vd,vd,A,1.d0,1.d0,1,2,1)
      call norm(A,1)
c     write(out,'(a20,3f10.5)') 'A(1)',(A(kk,1),kk=1,3)
      A(3,3)=1.d0
      call prodv(A,A,A,3,1,2)
      call norm(A,2)
      call prodv(A,A,A,2,3,1)
      call norm(A,1)
      call traspo(A,A,1,1)
      call r_frame(A,XO,N)
c     write(out,'(i4,3f9.5)')(lll,(XO(kk,lll),kk=1,3),lll=1,N)
      call azzera(B,0.d0,9)
      B(1,1)= 1.d0
      B(2,2)=-1.d0
      B(3,3)=-1.d0
      call verify(XO,B,MK,MN,MV,N)
c     write(out,*)'asse C2 MV,NMS',MV,NMS
      if(MV.eq.1)go to 2310
      call traspo(A,A,1,1)
      call r_frame(A,XO,N)
      call azzera(A,0.d0,9)
      A(1,2)=-1.d0
      A(2,1)= 1.d0
      A(3,3)= 1.d0
      call r_frame(A,XO,N)
      call verify(XO,B,MK,MN,MV,N)
c     write(out,*)'asse C2 MV,NMS',MV,NMS
      if(MV.eq.1)go to 2310
      call traspo(A,A,1,1)
      call r_frame(A,XO,N)
      go to 2040
c     ico=.true.
c     call compl_gr(MK,N,*5000)
c     write(out,*)'asse C2,NMS',NMS
Controlla il centro
 2310 call azzera(B,0.d0,9)
      B(1,1)=-1.d0
      B(2,2)=-1.d0
      B(3,3)=-1.d0
      call verify(XO,B,MK,MN,MV,N)
c     write(out,*)'CENTRO MV,NMS',MV,NMS
Costruzione dell'asse C5 a |x|=63.43494882 da z e la cui priezione e' a 18 da x
      call azzera(A,0.d0,9)
      ROT=-63.43494882d0/RAD
 3000 CA=COS(ROT)
      CB=SIN(ROT)
      call azzera(A,0.d0,9)
      A(1,1)=1.d0
      A(2,2)= CA
      A(3,3)= CA
      A(2,3)=-CB
      A(3,2)= CB
c rotazione 18 su z
      call azzera(B,0.d0,9)
      ROT18=-18.d0/RAD
 3100 CA=COS(ROT18)
      CB=SIN(ROT18)
      call azzera(B,0.d0,9)
      B(1,1)= CA
      B(2,2)= CA
      B(1,2)=-CB
      B(2,1)= CB
      B(3,3)=1.d0
      call prodmm(B,A,A,1,1,1)
c     write(out,*)' matrice di rotazione'
c     write(out,'(a12,2f10.6)')'ROT,ROT18',ROT,ROT18
c     write(out,'(3f10.6)')((A(i,k),k=1,3),i=1,3)
      call prodmm(SIM,A,B,2,1,1)
      call traspo(A,A,1,1)
      call prodmm(A,B,A,1,1,1)
      call verify(XO,A,MK,MN,MV,N)
c     write(out,*)'secondo asse C5 ,MV,NMS',MV,NMS
      if(MV.eq.1)go to 4000
      if(ROT18.lt.0.d0)then
      ROT18=-ROT18
      go to 3100
      endif
      if(ROT.gt.0.d0)return
      ROT=-ROT
      go to 3000
 4000 call compl_gr(MK,N,*5000)
 5000 return 1
      END
      subroutine intpe(RIGA,B,C,IA,NB,*,*)
      implicit double precision (a-h,o-z)
      CHARACTER*(*) RIGA
      CHARACTER*23 CAR,CARA
C INTERPRETA IL VETTORE  A COME POSIZIONE EQUIVALENTE O
C  PREPARA LO STESSO PER LA STAMPA A PARTIRE DA B E C
      DIMENSION B(3,3,1),C(3,1)
      CAR='1234567890 /,+-xyz[]XYZ'
      GO TO (10,500),NB
C INTERPRETAZIONE DI A
   10 MCM=1
      LR=LEN(RIGA)
      call compatta(RIGA,LR,KM)
      K=1
      COST=1.d0
      DO 15 I=1,3
      C(I,IA)=0.d0
      DO 15 J=1,3
   15 B(I,J,IA)=0.d0
      DO 200 I=1,KM
C IDENTIFICAZIONE DEL CARATTERE
      CARA=RIGA(I:I)
      DO 20 J=1,23
      IF(CAR(J:J).EQ.CARA)GO TO 30
   20 CONTINUE
   25 return 2
   30 IF(J-11)40,200,100
C CARATTERE NUMERICO
   40 IF(J.EQ.10)J=0
      GO TO(50,60),MCM
C NUMERATORE
   50 IF(C(K,IA).NE.0.d0)C(K,IA)=C(K,IA)*10
      C(K,IA)=C(K,IA)+COST*J
      COST=1.d0
      GO TO 200
C DENOMINATORE
   60 C(K,IA)=C(K,IA)/J
      MCM=1
      GO TO 200
CARATTERI ALFABETICI
  100 if(J.gt.20)J=J-5
      J=J-11
      GO TO(110,120,130,140,150,150,150,200,210),J
C BARRA
  110 MCM=2
      GO TO 200
CAMBIO RIGA MATRICE (VIRGOLA)
  120 K=K+1
      GO TO 200
C SEGNO +
  130 COST=1.d0
      GO TO 200
C SEGNO -
  140 COST=-1.d0
      GO TO 200
C  X,Y O Z
  150 J=J-4
      B(K,J,IA)=COST
      COST=1.d0
  200 CONTINUE
  210 RETURN 1
C SCRITTURA  DEL VETTORE POSIZIONE EQUIVALENTE ALFANUMERICA
  500 KM=LEN(RIGA)
      DO 510 I=1,KM
  510 RIGA(i:i)=car(11:11)
      RIGA(12:12)=','
      RIGA(23:23)=','
      DO 1000 I=1,3
      MCM=11*I
      DO 600 J=1,3
      L=4-J
      IF(NINT(B(I,L,IA)).EQ.0)GO TO 600
      MA=L+15
      RIGA(MCM:MCM)=CAR(MA:MA)
      K=14
      IF(B(I,L,IA).LT.(-0.1d0))K=15
      MCM=MCM-1
      RIGA(MCM:MCM)=CAR(K:K)
      MCM=MCM-1
  600 CONTINUE
      MS=11
      IF(ABS(C(I,IA)).LT.0.1d0)GO TO 1000
      IF(C(I,IA).LT.0.d0)MS=15
      DO 700 J=1,6
      CA=C(I,IA)*J
      K=NINT(CA)
      IF(ABS(FLOAT(K)-CA).LT.0.1d0)GO TO 800
  700 CONTINUE
  800 IF(J.EQ.1)GO TO 900
      RIGA(MCM:MCM)=CAR(J:J)
      MCM=MCM-1
      RIGA(MCM:MCM)=CAR(12:12)
      MCM=MCM-1
  900 J=IABS(MOD(K,10))
      RIGA(MCM:MCM)=CAR(J:J)
      K=K/10
      MCM=MCM-1
      IF(K.GT.0)GO TO 900
      RIGA(MCM:MCM)=CAR(MS:MS)
 1000 CONTINUE
      KM1=KM-1
      call compatta(RIGA(2:KM),KM1,K)
      if(RIGA(2:2).eq.'+')RIGA(2:2)=' '
      do 1010 i=1,2
      k=index(RIGA,',+')
      k1=k+1
      if(k.ne.0)RIGA(k1:k1)=' '
 1010 continue
      call compatta(RIGA(2:KM),KM1,K)
      RETURN 1
      END
      subroutine inv33(X,Y,K,L)
      implicit double precision (a-h,o-z)
      DIMENSION X(3,3,1),Y(3,3,1)
C INVERSIONE DI UNA MATRICE 3X3.
      C=1.0/det(X,K)
      Y(1,1,L)=(X(2,2,K)*X(3,3,K)-X(2,3,K)*X(3,2,K))*C
      Y(2,1,L)=(X(2,3,K)*X(3,1,K)-X(2,1,K)*X(3,3,K))*C
      Y(3,1,L)=(X(2,1,K)*X(3,2,K)-X(2,2,K)*X(3,1,K))*C
      Y(1,2,L)=(X(1,3,K)*X(3,2,K)-X(1,2,K)*X(3,3,K))*C
      Y(2,2,L)=(X(1,1,K)*X(3,3,K)-X(1,3,K)*X(3,1,K))*C
      Y(3,2,L)=(X(1,2,K)*X(3,1,K)-X(1,1,K)*X(3,2,K))*C
      Y(1,3,L)=(X(1,2,K)*X(2,3,K)-X(1,3,K)*X(2,2,K))*C
      Y(2,3,L)=(X(1,3,K)*X(2,1,K)-X(1,1,K)*X(2,3,K))*C
      Y(3,3,L)=(X(1,1,K)*X(2,2,K)-X(1,2,K)*X(2,1,K))*C
      RETURN
      END
      FUNCTION ium(A,B,C,M,N)
CONTROLLA  SE A=B A MENO DI C
C VERO ium=1      FALSO  ium=0
      implicit double precision (a-h,o-z)
      DIMENSION A(3,3,1),B(3,3,1)
      ium=0
      DO 1 I=1,3
      DO 1 J=1,3
      IF(ABS(A(I,J,M)-B(I,J,N)).GT.C)RETURN
    1 CONTINUE
      ium=1
      RETURN
      END
      subroutine maiuscol(c)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character c
      do 1 i=97,122
      if(c.ne.CHAR(i))go to 1
      c=CHAR(i-32)
      return
    1 continue
      return
      end
      subroutine mass(NAT,ATOM)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      integer out
      COMMON/output/out
      CHARACTER*2 ATOM(NMA)
      CHARACTER*2 SIMBO(103),ITT
      REAL*8 MASSA(103),RAGGI(103)
      COMMON/AT1/X(3,NMA),AMAS(NMA),RCOV(NMA),MSP(NMA)
      DATA (MASSA(i),i=1,2)/1.00794,4.002602/
      DATA (MASSA(i),i=3,10)/6.941,9.012182,10.811,12.011,14.00674,
     115.9994,18.9984032,20.1797/
      DATA (MASSA(i),i=11,18)/22.989768,24.3050,26.981539,28.0855,
     130.97362,32.066,35.4527,39.948/
      DATA (MASSA(i),i=19,36)/39.0983,40.078,44.955910,47.88,
     1 50.9415,51.9961,54.93085,55.847,58.93320,58.69,63.546,65.39,
     2 69.723,72.61,74.92159,78.96,79.904,83.80/
      DATA (MASSA(i),i=37,54)/85.4678,87.62,88.90585,91.224,
     1 92.90638,95.94,98,101.07,102.90550,106.42,107.8682,112.411,
     2 114.82,118.710,121.75,127.60,126.90447,131.29/
      DATA (MASSA(i),i=55,86)/132.90543,137.327,138.9055,140.115,
     1 140.90765,144.24,145,150.36,151.965,157.25,158.92534,162.50,
     2 164.93032,167.26,168.93421,173.04,174.967,178.49,180.9479,
     3 183.85,186.207,190.2,192.22,195.08,196.96654,200.59,204.3833,
     4 207.2,208.98037,209,210.0,222/
      DATA (MASSA(i),i=87,103)/223,226.025,227.028,232.0381,
     1 231.03588,238.0289,237.048,244.,243.,247.,247.,251.,252.,
     2 257.,258.,259.,260.0/
      DATA (RAGGI(i),i=1,2)/.32,.93/
      DATA (RAGGI(i),i=3,10)/1.23,0.90,.82,.770,.75,.73,.72,.71/
      DATA (RAGGI(i),i=11,18)/1.54,1.4,1.18,1.11,1.06,1.02,.99,.98/
      DATA (RAGGI(i),i=19,36)/2.03,1.74,1.44,1.32,1.22,1.18,1.17,1.17,
     1 1.16,1.15,1.17,1.25,1.26,1.22,1.20,1.16,1.14,1.12/
      DATA (RAGGI(i),i=37,54)/2.16,1.91,1.62,1.45,1.34,1.30,1.27,1.25,
     1 1.25,1.28,1.34,1.48,1.44,1.41,1.40,1.36,1.33,1.31/
C per i raggi degli elementi dal Ce al Lu ho interpolato
      DATA (RAGGI(i),i=55,86)/2.35,1.98,1.69,1.67,1.65,1.63,1.61,1.59,
     1 1.57,1.55,1.53,1.51,1.49,1.47,1.45,1.45,1.45,1.44,1.34,1.30,1.28,
     2 1.26,1.27,1.30,1.34,1.49,1.48,1.47,1.46,1.76,1.45,1.45/
C per i raggi degli elementi dal Ce a Lu ho interpolato
      DATA (RAGGI(i),i=87,103)/2.35,1.98,1.69,1.67,1.65,1.63,1.61,1.59,
     1 1.57,1.55,1.53,1.51,1.49,1.47,1.45,1.45,1.45/
      DATA (SIMBO(i),i=1,2)/' H','HE'/
      DATA (SIMBO(i),i=3,10)/'LI','BE',' B',' C',' N',' O',' F','NE'/
      DATA (SIMBO(i),i=11,18)/'NA','MG','AL','SI',' P',' S','CL','AR'/
      DATA (SIMBO(i),i=19,36)/' K','CA','SC','TI',' V','CR','MN','FE',
     1'CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR'/
      DATA (SIMBO(i),i=37,54)/'RB','SR',' Y','ZR','Nb','MO','TC','RU',
     1'RH','PD','AG','CD','IN','SN','SB','TE',' I','XE'/
      DATA (SIMBO(i),i=55,86)/'CS','BA','LA','CE','PR','ND','PM',
     1'SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',' W',
     2'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN'/
      DATA (SIMBO(i),i=87,103)/'FR','RA','AC','TH','PA',' U','NP',
     1'PU','AM','CM','BK','CF','ES','FM','MD','NO','LR'/
      DO J=1,NAT
        DO I=1,103
          IF(ATOM(J).EQ.SIMBO(I)) THEN
            AMAS(J)=MASSA(I)
            RCOV(J)=RAGGI(I)
c           write(out,'(a2,2f7.3)')ATOM(J),AMAS(J),RCOV(J)
            GOTO 2
          ENDIF
        ENDDO
        write(out,*)
        write(out,3) ATOM(J)
    3   format(' ERROR: MISSING DATA FOR THE ATOM-',A2,'-')
        STOP
    2   CONTINUE
      ENDDO
       DO J=1,NAT
       MSP(J)=0
       ENDDO
        M=1
        ITT=ATOM(1)
 10       DO 1 K=1,NAT
          IF(ATOM(K).EQ.ITT)MSP(K)=M
 1       CONTINUE
        M=M+1
        DO J=1,NAT
        IF(MSP(J).EQ.0)THEN
        ITT=ATOM(J)
        MSP(J)=M
        GOTO 10
        ENDIF
        ENDDO
      RETURN
        END
      subroutine momin(XO,P,RO,N)
      implicit double precision (a-h,o-z)
c da un set di coordinate XO che vengono moltiplicate per la matrice RO
C che puo' essere una matrice di ortogonalizzazione calcola i momenti di
C inerzia utilizzando i pesi P
      COMMON/ORIE/OR(3,3),OT(3,3),OTI(3,3),BARC(3),BARO(3),RIN(3)
      DIMENSION RO(3,3),AA(3,3),AMOM(6),CAM(3),ORI(3,3)
      DIMENSION C1(3),C2(3),C3(3),C4(3),C5(3),C6(3),IC1(3),IC2(3),IC3(3)
      DIMENSION XO(3,1),P(1)
      DO 100 I=1,3
      RIN(I)=0.d0
      AMOM(I)=0.d0
  100 AMOM(I+3)=0.d0
      SUM=0.d0
      DO 110 I=1,N
      call prodmv(RO,XO,XO,1,I,I)
      PESO=P(I)
      RIN(1)=RIN(1)+PESO*XO(1,I)
      RIN(2)=RIN(2)+PESO*XO(2,I)
      RIN(3)=RIN(3)+PESO*XO(3,I)
  110 SUM=SUM+PESO
      DO 115 I=1,3
      RIN(I)=RIN(I)/SUM
  115 BARO(I)=RIN(I)
      DO 120 I=1,N
      PESO=P(I)
      XO(1,I)=XO(1,I)-BARO(1)
      XO(2,I)=XO(2,I)-BARO(2)
      XO(3,I)=XO(3,I)-BARO(3)
      CAM(1)=XO(1,I)*PESO
      CAM(2)=XO(2,I)*PESO
      CAM(3)=XO(3,I)*PESO
      AMOM(1)=AMOM(1)+CAM(1)*XO(1,I)
      AMOM(2)=AMOM(2)-CAM(1)*XO(2,I)
      AMOM(3)=AMOM(3)-CAM(1)*XO(3,I)
      AMOM(4)=AMOM(4)+CAM(2)*XO(2,I)
      AMOM(5)=AMOM(5)-CAM(2)*XO(3,I)
      AMOM(6)=AMOM(6)+CAM(3)*XO(3,I)
  120 CONTINUE
      AA(1,1)=AMOM(4)+AMOM(6)
      AA(2,2)=AMOM(1)+AMOM(6)
      AA(3,3)=AMOM(1)+AMOM(4)
      AA(1,2)=AMOM(2)
      AA(2,1)=AMOM(2)
      AA(1,3)=AMOM(3)
      AA(3,1)=AMOM(3)
      AA(2,3)=AMOM(5)
      AA(3,2)=AMOM(5)
      call eigen(AA,OR,RIN,C1,C2,C3,C4,C5,C6,IC1,IC2,IC3,3)
C__________________
C modifica per rendere OR machine-independent
      am1=-1
      am2=-1
      do i=1,3
        if(DABS(OR(i,1)).gt.am1)then
          am1=DABS(OR(i,1))
          n1=i
        endif
        if(DABS(OR(i,2)).gt.am2)then
          am2=DABS(OR(i,2))
          n2=i
        endif
      enddo
      if(OR(n1,1).lt.0.0)then
        do i=1,3
          OR(i,1)=-OR(i,1)
        enddo
      endif
      if(OR(n2,1).lt.0.0)then
        do i=1,3
          OR(i,2)=-OR(i,2)
        enddo
      endif
      CALL prodv(OR,OR,OR,1,2,3)
C__________________
      call traspo(OR,OR,1,1)
      call prodmm(OR,RO,OT,1,1,1)
      call inv33(OT,OTI,1,1)
      call inv33(RO,ORI,1,1)
      call prodmv(ORI,BARO,BARC,1,1,1)
      do 130 i=1,N
  130 call prodmv(OR,XO,XO,1,I,I)
      RETURN
      END
      subroutine norm(A,I)
      implicit double precision (a-h,o-z)
C NORMALIZZAZIONE DEL  VETTORE A
      DIMENSION A(3,1)
      call prods(A,A,B,I,I,2)
      DO 1 J=1,3
    1 A(J,I)=A(J,I)/B
      RETURN
      END
      subroutine opt_axis(XO,PESO,V3,MK,N,IAX)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      DIMENSION XO(3,NMA),MK(NMA,NMG),V3(3),PESO(NMA),VS(3)
      DIMENSION IV(NMA),CO(3,NMA),W(NMA),nge(NMA)
C fa in modo che gli la somma delle distanze dall'elemento di simmetria IAX
C sia minimo
      do i=1,N
        iv(i)=0
      enddo
      ng=0
      cost=det(SIM,IAX)
 1000 do I=1,N
         if(iv(i).eq.0)go to 1100
      enddo
      go to 1300
 1100 ng=ng+1
C ng = numero del gruppo di atomi equivalente rispetto all'operazione IAX
      iv(i)=ng
      do j=1,3
         CO(j,ng)=XO(j,I)
      enddo
      W(ng)=PESO(I)
      nge(ng)=1
      k=I
 1200 l=MK(k,IAX)
      if(l.eq.I)go to 1000
      iv(l)=ng
C somma tutti gli atomi del gruppo con segno opportuno
      do j=1,3
         CO(j,ng)=CO(j,ng)+cost*XO(j,l)
      enddo
      W(ng)=W(ng)+PESO(l)
      nge(ng)=nge(ng)+1
      k=l
      go to 1200
C ha assegnato tutti gli atomi ad un gruppo ed il peso ai gruppi 
 1300 continue    
      do i=1,ng
      W(I)=W(i)/nge(I)
      enddo
      call azzera(V3,0.d0,3)
      dism=0.d0
cerca il gruppo piu' lontano dall'elemento di simmetria (im)
      do i=1,ng 
         call prodmv(SIM,CO,VS,IAX,i,1)
         call combli(CO,VS,VS,1.d0,cost,i,1,1)
         call prods(VS,VS,dis,1,1,1)
c        write(out,'(i3,6f10.5)')i,(CO(j,i),j=1,3),VS
         if(dis.gt.dism)then
            dism=dis
            im=i
         endif
      enddo
      call azzera(V3,0.d0,3)
         call prods(CO,CO,dim,im,im,2)
c      write(out,'(6him,dim,i3,3f10.5)')im,dim
c somma fra loro tutti i gruppi con un segno determinato dall'angolo
c fra il gruppo in esame ed il gruppo im
      do i=1,ng
         call prods(CO,CO,abcos,im,i,1)
         call prods(CO,CO,di,i,i,2)
         absabcos=abs(abcos)
         if(absabcos.gt.1.d-5)then
            cost=W(I)*abcos/absabcos
            call combli(V3,CO,V3,1.d0,cost,1,i,1)
         endif
      enddo
c il vettore risultante e' normalizzato e diventera' una riga della
c matrice di rotazione definitiva
      call norm(V3,1)
      if(V3(1).gt.0.d0)return
      V3(1)=-V3(1)
      V3(2)=-V3(2)
      V3(3)=-V3(3)
      return
      end
      subroutine out_bond(XO,IEQAT,IASU,N)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      integer out
      COMMON/output/out
      CHARACTER*6 NAME
      COMMON/AT4/NAME(NMA)
      COMMON/AT5/LEGAMI(20,NMA),NLEG(NMA)
      DIMENSION T(3),V(3),IASU(NMA),tableg(20,NMA),XO(3,NMA),IEQAT(NMA)
      CHARACTER*80 riga
      write(out,*)
c     write(out,'(a6,3f10.6)')(NAME(IEQAT(I)),(XO(k,i),k=1,3),i=1,N)
      RAD=57.29577951308232D0
      riga=' '
      nb=0
      do 1400 i=1,N
      if(IASU(i).ne.2)go to 1400
      N1=NLEG(i)
      do 1300 k=1,N1
      l=LEGAMI(k,i)
      call combli(XO,XO,T,1.d0,-1.d0,i,l,1)
      call prods(t,t,tableg(k,i),1,1,2)
      nc=nb*26+1
      nb=nb+1
      write(riga(nc:nc+25),'(a6,a1,a6,f7.4,6x)')NAME(IEQAT(I)),'-',
     *NAME(IEQAT(l)),tableg(k,i)
      if(nb.eq.3)then
        write(out,'(a)')riga
        nb=0
        riga=' '
      endif
 1300 continue
 1400 continue
      if(nb.ne.0)write(out,'(a)')riga
      write(out,*)
        riga=' '
        nb=0
      do 2000 i=1,N
      N1=NLEG(i)
      if(IASU(i).ne.2)go to 2000
      do 1900 k=1,N1-1
      k1=k+1
      i1=LEGAMI(k,i)
      call combli(XO,XO,T,1.d0,-1.d0,i1,i,1)
      do 1800 l=k1,N1
      i2=LEGAMI(l,i)
      call combli(XO,XO,V,1.d0,-1.d0,i2,i,1)
      call prods(T,V,cost,1,1,1)
      cost=cost/(tableg(k,i)*tableg(l,i))
      CALL arccos(cost,ang)
      nc=nb*35+1
      nb=nb+1
      write(riga(nc:nc+27),'(a6,2(a1,a6),f8.3,7x)')NAME(IEQAT(I1)),
     *'-',NAME(IEQAT(I)),'-',NAME(IEQAT(I2)),ang
      if(nb.eq.2)then
        write(out,'(a)')riga
        nb=0
        riga=' '
      endif
 1800 continue
 1900 continue
 2000 continue
      if(nb.ne.0)write(out,'(a)')riga
      write(out,*)
      return
      END
      subroutine pratba(B,A,C,I,J,K)
      implicit double precision (a-h,o-z)
C                    T
C TRASFORMAZIONE  C=A BA
      DIMENSION A(3,3,1),B(3,3,1),C(3,3,1),D(3,3)
      call PRODMM(B,A,C,I,J,K)
      call TRASPO(A,D,J,1)
      call PRODMM(D,C,C,1,K,K)
      RETURN
      END
      subroutine prodmm(A,B,C,I,J,N)
      implicit double precision (a-h,o-z)
      DIMENSION A(3,3,1),B(3,3,1),C(3,3,1),X(3,3)
      DO 1 K=1,3
      DO 1 L=1,3
      X(K,L)=0.0d0
      DO 1 M=1,3
      Y=A(K,M,I)
      Z=B(M,L,J)
    1 X(K,L)=X(K,L)+Y*Z
      DO 2 K=1,3
      DO 2 L=1,3
    2 C(K,L,N)=X(K,L)
      RETURN
      END
      subroutine prodms(A,B,C,M,N)
      implicit double precision (a-h,o-z)
      DIMENSION A(3,3,1), B(3,3,1)
      X=C
      DO 1 I=1,3
      DO 1 J=1,3
      Y=A(I,J,M)
    1 B(I,J,N)=X*Y
      RETURN
      END
      subroutine prodmv(A,B,C,I,J,K)
      implicit double precision (a-h,o-z)
      DIMENSION A(3,3,1),B(3,1),C(3,1),X(3)
      DO 1 L=1,3
      X(L)=0.0d0
      DO 1 M=1,3
      Y=A(L,M,I)
      Z=B(M,J)
    1 X(L)=X(L)+Y*Z
      DO 2 L=1,3
    2 C(L,K)=X(L)
      RETURN
      END
      subroutine prods(A,B,C,I,J,K)
      implicit double precision (a-h,o-z)
C PRODOTTO SCALARE A*B=C (K=1)
C  MODULO VETTORE A        (K=2)
      DIMENSION A(3,1),B(3,1)
      X=0.0d0
      DO 1 L=1,3
      Y=A(L,I)
      Z=B(L,J)
    1 X=X+Y*Z
      C=X
      IF(K.EQ.2)C=DSQRT(X)
      RETURN
      END
      subroutine prodv(A,B,C,I,J,K)
      implicit double precision (a-h,o-z)
C PRODOTTO VETTORIALE A*B=C
      DIMENSION A(3,1),B(3,1),C(3,1)
      X1=A(1,I)
      X2=A(2,I)
      X3=A(3,I)
      Y1=B(1,J)
      Y2=B(2,J)
      Y3=B(3,J)
      C(1,K)=X2*Y3-X3*Y2
      C(2,K)=-X1*Y3+X3*Y1
      C(3,K)=X1*Y2-X2*Y1
      RETURN
      END
      subroutine r_frame(CO,XO,N)
C ruota il riferimento e le coordinate 
      implicit double precision (a-h,o-z)
      COMMON/ORIE/OR(3,3),OT(3,3),OTI(3,3),BARC(3),BARO(3),RIN(3)
      DIMENSION CO(3,3),XO(3,*)
c     write(6,*)'rotazione del riferimento'
c     write(6,'(3f10.6)')((CO(i,k),k=1,3),i=1,3)
      DO 1010 I=1,N
 1010 call prodmv(CO,XO,XO,1,I,I)
      call prodmm(CO,OR,OR,1,1,1)
      call prodmm(CO,OT,OT,1,1,1)
      call inv33(OT,OTI,1,1)
c     write(6,*)'nuove coordinate'
c     write(6,'(i5,3f10.6)')(i,(XO(k,i),k=1,3),i=1,N)
      return 
      end
      subroutine rms_min(tm)
      implicit double precision (a-h,o-z)
C
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      PARAMETER (NMV=NMA*NMG)
      common/rmsmin/RV(3,3),ppo(3,NMV),ppu(3,NMV),npu
      dimension tm(3),t(3)
      passo=0.02
      call azzera(RV,0.d0,9)
      call azzera(tm,0.d0,3)
      rm=crms(tm)
    1 ind1=0
      do i1=1,3
       t(1)=tm(1)+(2-i1)*passo
       do i2=1,3
        t(2)=tm(2)+(2-i2)*passo
        do i3=1,3
         t(3)=tm(3)+(2-i3)*passo
         arms=crms(t)
C        write(6,*)t,arm,i1,i2,i3
         if(arms.lt.rm)then
          ind1=1
          rm=arms
c         write(6,*)'rm,passo',rm,passo
          call trasfv(t,tm,1,1)
         endif
        enddo
       enddo
      enddo
      if(ind1.eq.1)goto 1
      if(passo.lt.0.0002)return
      passo=passo*0.5
      go to 1
      end
      subroutine s_coor(SI,XO,XS,CO,D,DEL,DLM,CSM,MK,ID,ISU,nes,N)
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      COMMON/SIMME1/RMS(3,NMA),RMST
      DIMENSION CO(3,NMA),XO(3,NMA),XS(3,NMA),SI(3,3,NMG),D(NMA)
      DIMENSION DEL(4),DLM(4),V3(3),ID(4),MK(NMA,NMG),ISU(NMG)
      DIMENSION dc(NMA),cf(NMA)
calcola le coordinate simmetrizzate XS
c     write(out,'(i3,3f10.5,)')(kk,(XO(ll,kk),ll=1,3),kk=1,N)
      NAZ=3*N
      call azzera(CO,0.d0,NAZ)
      call azzera(XS,0.d0,NAZ)
      call azzera(RMS,0.d0,NAZ)
      call azzera(DEL,0.d0,4)
      call azzera(DLM,0.d0,4)
      call azzera(cf,0.d0,N)
C ISU=indice delle matrici del sottogruppo in considerazione
C nes=numero di elementi del sottogruppo
c     write(out,'(a3,12i3)')'ISU',(ISU(J),J=1,nes)
      DO 2800 I=1,N
      call prods(XO,XO,D(I),I,I,2)
      DO 2790 J=1,nes
       K=MK(I,ISU(j))
       call prodmv(SI,XO,V3,ISU(J),I,1)
       call combli(XS,V3,XS,1.d0,1.d0,K,1,K)
 2790 continue
c     write(out,'(i2,4f9.5)')I,(XO(kk,I),kk=1,3),D(I)
 2800 CONTINUE
      do i=1,N
       do J=1,3
        XS(J,I)=XS(J,I)/dfloat(nes)
       enddo
       call prods(XS,XS,dc(i),i,i,2)
      enddo
      do I=1,N
        DO J=1,nes
          K=MK(I,ISU(j))
          cf(I)=cf(I)+D(K)
        enddo
        if(dc(i).gt.1.d-01)cf(i)=cf(i)/(dc(i)*nes)
        do J=1,3
          XS(J,I)=XS(J,I)*cf(i)
        enddo
      enddo
c     write(out,'(i2,6f10.6)')1,(XS(kk,1),kk=1,3),dc(1),cf(1)
C il fattore di correzione cf=media delle distanze atomi equivalenti baricentro,
C corregge le coordinate simmetrizzate facendo si che queste equivalgano alla
C media per pura rotazione degli atomi equivalenti
calcola le rms e gli scostamenti medi (DEL) e massimi DLM
      do j=1,4
       ID(j)=0
      enddo
      DO 2820 i=1,N
       do j=1,3
        V3(j)=XO(j,i)-XS(j,i)
        com=dabs(V3(j))
        DEL(j)=DEL(J)+com
        if(com.gt.DLM(j))then
          DLM(j)=com
          ID(j)=I
        endif
       enddo
C massima deviazione dalla simmetria in angstrom
      call prods(V3,V3,com,1,1,2)
      DEL(4)=DEL(4)+com
      if(DLM(4).lt.com)then
       DLM(4)=com
       ID(4)=I
      endif
 2820 continue
calcola la continuous symmetry measure (CSM)
      CSM=0.d0
      do j=1,nes
       do i=1,N
        K=MK(i,ISU(j))
        call prodmv(SI,XO,V3,ISU(J),i,1)
        call combli(XS,V3,V3,1.d0,-1.d0,K,1,1)
        RMS(1,K)=RMS(1,K)+V3(1)*V3(1)
        RMS(2,K)=RMS(2,K)+V3(2)*V3(2)
        RMS(3,K)=RMS(3,K)+V3(3)*V3(3)
        call prods(V3,V3,com,1,1,1)
        CSM=CSM+com
       enddo
      enddo
      cost=1.d0/dfloat(nes)
      RMST=0.d0
      do I=1,N
       do j=1,3
       RMST=RMST+RMS(j,i)
       RMS(j,i)=sqrt(RMS(j,i)*cost)
       enddo
      enddo
      RMST=sqrt(RMST/dfloat(nes*N))
      CSM=CSM*1.d2/dfloat(nes*N)
      return
      end
      subroutine schoenfl(MDEG)
C From the group matrices gives the Schoenflies symbol
      implicit double precision (a-h,o-z)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CMS(NMG),CSMT,NMS,MTG(NMG,NMG)
      COMMON/SIMME1/RMS(3,NMA),RMST
      character*3 lbls,pt*7
      character S1,S2,Sn*2,tot*4
      COMMON/SIMME3/IPOGROU,lbls(NMG)
      character*5 pointgr(32)
      DATA (pointgr(i),i=1,32)
     */'C1 ','Ci ','C2 ','Cs ','C2h','D2 ','C2v','D2h','C4 ','S4 ','C4h'
     *,'D4 ','C4v','D2d','D4h','C3 ','C3i','D3 ','C3v','D3d','C6 ','C3h'
     *,'C6h','D6 ','C6v','D3h','D6h','T  ','Th ','O  ','Td ','Oh '/
      IPOGROU=0
      maxasp=0
      maxasi=0
      mplane=0
      invers=0
      nplane=0
      lbls(1)='E'
C_______________________________________________________
C gruppi lineari
C per questi gruppi non stampo le matrici di simmetria ma solo il gruppo
C
      if(MDEG.eq.5)then
        pt='C(inf)v'
        do i=1,NMS
          call ax_order(SIM,i,m,msign,inv)
          if(inv.eq.1)go to 1000
        enddo
        NMS=1
        go to 5200
 1000   pt='D(inf)h'
        NMS=1
        go to 5200
      endif
C_______________________________________________________
      if(NMS.eq.1)then
        pt='C1 '
        IPOGROU=1
        return
      endif
      do i=2,NMS
        tot='    '
        S1=' '
        S2=' '
        Sn='  '
        call ax_order(SIM,i,m,msign,inv)
        if(m.gt.maxasp.and.msign.eq.1)maxasp=m
        if(m.gt.maxasi.and.msign.eq.-1)maxasi=m
        if(m.eq.2.and.mplane.eq.0)mplane=1
c     write(out,'(a20,4i3)')'i,msign,m,inv',i,msign,m,inv
        S1='C'
        if(m.eq.2)then
Centro
          if(inv.eq.1)then
            invers=1
            Sn='i '
            go to 1100
          endif
C piano
          Sn='2 '
C asse 2
          if(msign.eq.-1)then
            Sn='s '
            nplane=nplane+1
          endif
          go to 1100
        endif
C asse improprio m>2
        write(Sn,'(i2)')m
        if(msign.eq.-1)S1='S'
 1100   tot=S1//Sn//S2
        call compatta(tot,4,k)
        lbls(i)=tot(1:3)
      enddo
c     write(out,*)'maxasp maxasi mplane nplane invers   MDEG    NMS'
c     write(out,'(7i7)')maxasp,maxasi,mplane,nplane,invers,MDEG,NMS
      pt='   '
      if(NMS.eq.48)pt='Oh '
      if(NMS.eq.60)pt='I  '
      if(NMS.eq.120)pt='Ih '
      if(pt.ne.'   ')go to 5100
      if(MDEG.eq.3)then
        if(NMS.eq.12)pt='T23'
        if(invers.eq.1)pt='Th '
        if(maxasi.eq.4)pt='Td '
        if(pt.eq.'   ')pt='O  '
        go to 5100
      endif
      if(NMS.eq.2)then
        pt=lbls(2)
        go to 5000
      endif
      S2=' '
      S1='C'
C Cn e S2n dove n=maxasp
      write(Sn,'(I2)')maxasp
      if(MOD(NMS,2).eq.0.and.maxasi.eq.NMS)then
        S1='S'
        write(Sn,'(I2)')maxasi
        if(NMS.ne.6)go to 5000
        pt='C3i'
        go to 5100
      endif
      if(NMS.eq.maxasp)go to 5000
      if(NMS.lt.maxasi)then
        write(Sn,'(I2)')maxasi
        S1='S'
        go to 5000
      endif
C Dn, Cnv e Cnh
      S1='D'
      if(NMS.eq.maxasp*2)then
        if(nplane.eq.0)go to 5000
        S1='C'
        S2='v'
        if(nplane.eq.maxasp)go to 5000
        S2='h'
        go to 5000
      endif
      S2='h'
      if(nplane.eq.maxasp)S2='d'
    1 format(/,' Schoenflies symbol = ',a7,'  CSM =',f8.4,5x,
     *'Molecular RMS =',f8.4)
    2 format(' CSM,see: Zabrodsky et al. (1993) JACS, 115, 8278-8298')
 5000 tot='    '
      tot=S1//Sn//S2
      call compatta(tot,4,k)
      pt=tot(1:3)
 5100 do m=1,32
        if(pt.eq.pointgr(m))IPOGROU=m
      enddo
 5200 write(out,1)pt,CSMT,RMST
      write(out,2)
      if(out.eq.9)then
         write(*,1)pt,CSMT,RMST
         write(*,2)
      endif
      RETURN
      END
      subroutine stasy
      implicit double precision (a-h,o-z)
C
C        STAMPA LE POSIZIONI EQUIVALENTI
C
c     character*120 linea
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      character*80 linea
      character*40 riga
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      dimension TRL(3,NMG)
      character*3 lbls
      COMMON/SIMME3/IPOGROU,lbls(NMG)
      call azzera(TRL,0.d0,3*NMS)
      i1=-39
      I=1
      if(NMS.gt.1)I=2
    2 format(2a40)
      write(out,2)(' Symmetry element its CSM and Max.Diff. ',k=1,I)
      DO 320 M=1,NMS
      mm=M
      call intpe(riga(1:39),SIM,TRL,mm,2,*310,*310)
  310 i1=i1+40
      i2=i1+39
      i3=I1+21
      i4=I1+37
      write(linea(i1:i2),3)M,lbls(M),riga
    3 format(I2,') [',a3,'] ',a30)
      write(linea(i3:i4),'(2f8.4)')CSM(M),DEV(M)
c     if(i1.eq.81)then
      if(i1.eq.41)then
          write(out,'(1x,a)')linea
      do 315 kk=1,80
c     do 315 kk=1,120
      linea(kk:kk)=' '
  315 continue
          i1=-39
      endif
  320 CONTINUE
      if(i1.ge.1)write(out,'(1x,a)')linea(1:i1+39)
      RETURN
      END
      subroutine trasfm(A,B,I,J)
      implicit double precision (a-h,o-z)
C TRASFERIMENTO DI UNA MATRICE
      DIMENSION A(3,3,1),B(3,3,1)
      DO 1 K=1,3
      DO 1 L=1,3
    1 B(K,L,J)=A(K,L,I)
      RETURN
      END
      subroutine trasfv(A,B,I,J)
      implicit double precision (a-h,o-z)
C TRASFERIMENTO DI UN VETTORE
      DIMENSION A(3,1),B(3,1)
      B(1,J)=A(1,I)
      B(2,J)=A(2,I)
      B(3,J)=A(3,I)
      RETURN
      END
      subroutine traspo(A,B,M,N)
      implicit double precision (a-h,o-z)
      DIMENSION A(3,3,1),B(3,3,1),C(3,3)
      DO 1 I=1,3
      DO 1 J=1,3
    1 C(I,J)=A(J,I,M)
      call trasfm(C,B,1,N)
      RETURN
      END
      subroutine verify(XO,A,MK,MN,MV,N)
      implicit double precision (a-h,o-z)
C VERIFICA SE A E' UNA MATRICE DI SIMMETRIA A MENO DI DXM
C IN CASO AFFERMATIVO PONE MV=1,  NMS=NMS + 1, SIM(I,J,NMS=A(I,J)
C  MK(I,NMS) CONTIENE IL NUMERO DELL'ATOMO OTTENUTO STASFORMANDO I CON
C  L'OPERAZIONE NMS
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      integer out
      COMMON/output/out
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      COMMON/AT2/SX(3,NMA),SIG(NMA),DXM(NMA),DCM,DCME,indwgh,indtol
      DIMENSION XO(3,1),MN(1),MK(NMA,1),A(3,3),V(3),W(3),iat(NMA)
      MV=1
Controlla per prima cosa che ls matrice A non sia gia' compresa in SIM
      do i=1,NMS
        if(ium(A,SIM,1.d-2,1,i).eq.1)go to 6
      enddo
      NMS=NMS+1
c     write(out,'(i3,3f10.6,i4)')(i,(XO(k,i),k=1,3),MN(i),i=1,N)
c     write(out,'(3f10.6)')((A(k,i),i=1,3),k=1,3)
      DO 4 I=1,N
      MK(I,NMS)=0
      JJ=0
      FMIN=DXM(I)
      call prodmv(A,XO,V,1,I,1)
c     write(out,'(3f10.6)')V
      DO 2 J=1,N
      IF(MN(I).NE.MN(J))GO TO 2
      call combli(XO,V,W,1.d0,-1.d0,J,1,1)
      call prods(W,W,F,1,1,2)
c     write(out,*)'I,J,NMS,F',I,J,NMS,F
      IF(F.GT.FMIN)GO TO 2
      FMIN=F
      JJ=J
    2 CONTINUE
c     write(out,*)'JJ, se JJ=0 rifiuta',JJ
      IF(JJ.EQ.0)GO TO 5
      MK(I,NMS)=JJ
    4 CONTINUE
c     write(out,*)'MK'
c     write(out,'(29i3)')(MK(k,NMS),k=1,N)
C analisi della matrice di moltiplicazione
      do k=1,N
      iat(k)=0
      enddo
      do 10 k=1,N
      if(iat(k).ne.0)go to 10
      iat(k)=1
      neq=1
      l=k
    7 m=MK(l,NMS)
      if(m.eq.k)go to 9
      iat(m)=1
      neq=neq+1
      if(neq.gt.N)go to 5
      l=m
      go to 7
    9 continue
      call ax_order(A,1,nord,msign,inv)
c     write(out,*)'nord,neq,inv',nord,neq,inv
C il numero di atomi equivalenti per l'operazione di ordine m puo' essere
C solo 1 o 2 (atomi giacenti sull'elemento di simmetria) m (caso normale)
C o 2m ( quando si ha un elemento riducibile 3/m 5/m etc.)
      if(neq.eq.1.or.neq.eq.2)go to 10
      if(neq.eq.nord)go to 10
      if(neq.eq.2*nord)go to 10
      go to 5
   10 continue
c     write(out,*)'accettata, nord=',nord
c     write(out,'(3f10.5)')((A(i,j),j=1,3),i=1,3)
      call trasfm(A,SIM,1,NMS)
      RETURN
    5 NMS=NMS-1
    6 MV=0
c     write(out,*)'return da verify,NMS=',NMS
      RETURN
      END
      subroutine work(MO)
C MO      = molecola da simmetrizzare
C X       = coordinate frazionarie originali
C NA      = numero totale di atomi
C MOL(I)  = molecola a cui appartiene l'atomo I
C MSP(I)  = indicatore della specie dell'atomo I
C AMAS(I) = massa dell'atomo I
C XO      = coordinate soli atomi da simmetrizzare
C N       = numero di atomi da simmetrizzare
C SIM     = matrici di simmetria del gruppo molecolare
C NMS     =  ordine del gruppo di simmetria trovato
      implicit double precision (a-h,o-z)
      PARAMETER (maxorder=8)
      PARAMETER (NMA=1000)
      PARAMETER (NMG=120)
      PARAMETER (NMV=NMA*NMG)
      common/rmsmin/RV(3,3),ppo(3,NMV),ppu(3,NMV),npu
      integer out
      COMMON/output/out
      CHARACTER*6 NAME,NNAME
      CHARACTER*52 comment
      COMMON/AT4/NAME(NMA)
      common/matcel/PC(7),PCR(7),O(3,3),OI(3,3),G(3,3),GI(3,3),CS(12)
      COMMON/AT0/NA,NMO
      COMMON/AT1/X(3,NMA),AMAS(NMA),RCOV(NMA),MSP(NMA)
      COMMON/AT2/SX(3,NMA),SIG(NMA),DXM(NMA),DCM,DCME,indwgh,indtol
      COMMON/AT3/MOL(NMA),MLG(NMA)
      COMMON/ORIE/OR(3,3),OT(3,3),OTI(3,3),BARC(3),BARO(3),RIN(3)
      COMMON/SIMME/SIM(3,3,NMG),DEV(NMG),CSM(NMG),CSMT,NMS,MTG(NMG,NMG)
      COMMON/SIMME1/RMS(3,NMA),RMST
      COMMON/SIMME2/AL(3,3),ALX(3,3)
      CHARACTER*3 lbls
      COMMON/SIMME3/IPOGROU,lbls(NMG)
      DIMENSION MK(NMA,NMG),MN(NMA),XO(3,NMA),XS(3,NMA),MD(NMA,2)
      DIMENSION D(NMA),AA(3,3),CO(3,NMA),IEQAT(NMA)
      DIMENSION PESO(NMA),DEL(4),DELM(4),IDM(4),delta(3,NMA)
      DIMENSION V(3),W(3),NNAME(4),ISU(NMG),IASU(NMA)
      CHARACTER type(2),car
      type(1)=' '
      type(2)='*'
      PIG=3.14159265358979D0
      deliner=0.1d0
      call azzera(SIM,0.d0,9*NMG)
      call azzera(AL,0.d0,9)
      SIM(1,1,1)=1.d0
      SIM(2,2,1)=1.d0
      SIM(3,3,1)=1.d0
      AL(1,1)=1.d0
      AL(2,2)=1.d0
      AL(3,3)=1.d0
      F=0.d0 
      IES=0 
      N=0 
      NMS=1
c selezione gruppo atomi da simmetrizzare
      write(out,1)MO
    1 format(//,38H           SYMMETRIZATION OF GROUP NR.,i2,//)
      if(out.eq.9)write(6,1)MO
      do 450 i=1,NA
      call azzera(AA,0.d0,9)
      AA(1,1)=SX(1,I)
      AA(2,2)=SX(2,I)
      AA(3,3)=SX(3,I)
      call prodmv(O,SX,AA,1,I,1)
      SIG(I)=0.d0
      den=0.d0
      do 435 k=1,3
        if(AA(K,K).le.0.d0)go to 435
        SIG(I)=SIG(I)+AA(k,k)
        den=den+1.d0
  435 continue
      if(den.ne.0.d0)then
       SIG(I)=SIG(I)/den
      else
       if(indtol.eq.3.or.indwgh.eq.3)then
        write(out,*)' WARNING: there is at least an atom with s.u.=0'
        write(out,*)
        write(out,*)'  CHANGE THE WEIGHTING AND TOLERANCE SCHEME!!!'
        STOP
       endif
      endif
      if(MOL(i).ne.MO)go to 450
      N=N+1
      IASU(i)=2
      IEQAT(N)=i
      do k=1,3
        XO(k,N)=X(k,i)
      enddo
      MN(N)=MSP(i)
      if(indwgh.eq.1) PESO(N)=AMAS(i)
      if(indwgh.eq.2) PESO(N)=1.0D0
      if(indwgh.eq.3)then
        if(SIG(i).lt..00001d0)then
          write(out,*)' ERROR: THE S.U. OF AN ATOM IS ZERO!!!'
          stop
        endif
        PESO(N)=1.d0/(SIG(I)*SIG(I))
      endif
      MK(N,1)=N
  450 continue
      porig=PESO(1)
      if(N.le.2)then
      write(out,*)' GROUP WITH LESS THAN THREE ATOMS'
      return
      endif
C
C     ORTOGONALIZZAZIONE DELLE COORDINATE DEGLI ATOMI DA SIMMETRIZZARE
C     CALCOLO DELLA DISTANZA MEDIA DAL BARICENTRO DEL GRUPPO
C
  460 call momin(XO,PESO,O,N)
c     do i=1,3
c     write(OUT,'(3f10.5,5x,3f10.5)')(OR(i,k),k=1,3),(OT(i,k),k=1,3)
c     enddo
C
C DM=maximum length of vectors XO(I)
C
      DM=0.d0
      DO 480 I=1,N
      call prods(XO,XO,D(I),I,I,2)
c     write(OUT,'(a9,i3,4f10.5)')'atom,D(I)',I,(XO(k,i),k=1,3),D(i)
      if(D(i).gt.DM)then
      DM=D(i)
      endif
  480 continue
C
      do 500 I=1,N
      DXM(I)=DCM
      if(indtol.eq.2)DXM(I)=DCM*D(I)/DM
      if(indtol.eq.3)DXM(I)=SIG(I)*DCM
  500 continue
      MDEG=0
      delr=DABS(RIN(1)-RIN(2))/RIN(1)
      IF(delr.LE.deliner)then
        MDEG=1
      endif
      delr=DABS(RIN(2)-RIN(3))/RIN(1)
      IF(delr.LE.deliner)then
        MDEG=MDEG+2
      endif
      NDEGA=MDEG+1
      nd=(MDEG+1)/2+1
      write(out,*)'PRINCIPAL INERTIA MOMENTS and DEGENERATION DEGREE'
      write(out,'(3f10.2,i10,/)')RIN,nd
C_______________________________________________________________________
C gruppi lineari
      if(MDEG.eq.1)then
        do i=1,N
          com=DSQRT(XO(1,i)*XO(1,i)+XO(2,i)*XO(2,i))
          if(com.gt.DXM(i))go to 990
        enddo
        call azzera(CO,0.d0,9)
        CO(1,1)=-1.d0
        CO(2,2)=-1.d0
        CO(3,3)= 1.d0
        call verify(XO,CO,MK,MN,MV,N)
        CO(3,3)=-1.d0
        call verify(XO,CO,MK,MN,MV,N)
        MDEG=5
        go to 2520
      endif
C_______________________________________________________________________
  990 GO TO (1020,1020,1000,2000),NDEGA
C
C ASSE UNICO=X  RUOTA LE COORDINATE IN MODO DA AVERE Z COME ASSE UNICO
C MODIFICA DI CONSEGUENZA ANCHE LA MATRICE DI ORIENTAZIONE
C
 1000 call azzera(CO,0.d0,9)
      MDEG=1
      CO(3,1)=1.d0 
      CO(1,3)=1.d0 
      CO(2,2)=-1.d0
      call r_frame(CO,XO,N)
      COM=RIN(3)
      RIN(3)=RIN(1) 
      RIN(1)=com
C ASSE UNICO Z
 1020 IU=3
      IB=2 
      IC=1
C RICERCA DELL'ASSE DI ORDINE MORD
 1100 MORD=maxorder+1
      IF(MDEG.EQ.0)MORD=3
      NASS=1
 1110 MORD=MORD-1
      IF(MORD.EQ.3.AND.MDEG.EQ.3)MORD=2
      IF(MORD.EQ.1)GO TO 1210
C GENERAZIONE DELLA MATRICE DI SIMMETRIA
 1120 COST=1.d0
      call azzera(CO,0.d0,9)
      ROT=2.d0*PIG/dfloat(MORD)
      CA=COS(ROT)
      CB=SIN(ROT)
 1140 CO(IB,IB)=CA*COST
      CO(IC,IC)=CO(IB,IB)
      CO(IB,IC)=CB*COST
      CO(IC,IB)=-CO(IB,IC)
      CO(IU,IU)=COST
      call verify(XO,CO,MK,MN,MV,N)
c     write(out,*)' 1140 tentativo con la matrice:',NINT(cost),MORD
c     write(out,'(3f9.5)')((CO(kkk,lll),lll=1,3),kkk=1,3)
c     write(out,*)'  NMS MORD NASS MV'
c     write(out,'(5i5)')NMS,MORD,NASS,MV
      IF(MV.EQ.1)GO TO 1150
      IF(COST.EQ.-1.d0)GO TO 1110
      COST=-1.d0
      mmez=MORD/2
      m2=2*mmez
      if(m2.ne.MORD)go to 1140
C se MORD=2*mmez ed mmez=dispari (assi -6=3/m -10=5/m)  si ha un
c elemento riducibile quindi e' inutile testarlo
      m2=2*mmez/2
      if(mmez-m2.eq.1)go to 1110
      GO TO 1140
C HA TROVATO UN ELEMENTO DI SIMMETRIA
 1150 IF(MORD.EQ.3.OR.MORD.EQ.6)IES=1
c     write(OUT,*)'MDEG,MORD',MDEG,MORD
      if(MDEG.eq.1.and.MORD.gt.2)go to 1300
C RICERCA ELEMENTI DI SIMMETRIA SU NUOVI ASSI
 1205 if(MORD.gt.2)go to 1110
 1210 NASS=NASS+1
C AL COMPLETAMENTO DEL GRUPPO
      IF(NASS.EQ.4)GO TO 2500
C     RICICLO SUL SECONDO O TERZO ASSE
C Buckyball: this write statement seems to fix a bug???
      write (6,*) NASS
      MORD=2
      I=IU 
      IU=IB 
      IB=IC 
      IC=I
      GO TO 1120
C RICERCA DELLA ROTAZIONE PER PORTARE EVENTUALI ELEMENTI DI SIMMETRIA
C A COINCIDERE CON GLI ASSI DI RIFERIMENTO
 1300 IF(N.le.1)then
        write(out,*)'Single atoms, not a molecule'
        return
      endif
C SELEZIONE DELL'ATOMO DA USARE PER DETERMINARE LA ROTAZIONE
      call azzera(CO,0.d0,18)
      CO(3,3)=1.d0
C divisione in gruppi equidistanti dall'asse unico
      do I=1,N
        D(I)=DSQRT(XO(1,I)*XO(1,I)+XO(2,I)*XO(2,I))
        MD(I,1)=0
        MD(I,2)=0
      enddo
      IAI=0
C IAI=NUMERO DI ENNUPLE DI ATOMI AVENTI DISTANZE COMPRESE IN UN RANGE
C         X< D <X+DXM*.5
C MD(I,1)=NUMERO DELLA ENNUPLA DI CUI L'ATOMO I FA PARTE
C MD(IAI,2)=N PER IL PRIMO ATOMO DELLA ENNUPLA DOVE N E' IL NUMERO DI
C         ATOMI DELLA STESSA
      N1=N-1
      do 1320 I=1,N1
      IF(MD(I,1).NE.0)GO TO 1320
      IAI=IAI+1
      MD(I,1)=IAI
      MD(IAI,2)=1
      J=I+1
      DO 1310 K=J,N
      IF(MN(I).NE.MN(K).or.K.eq.I.or.MD(K,1).ne.0)GO TO 1310
      IF(ABS(D(I)-D(K)).GT.0.5d0*DXM(I)) GO TO 1310
      MD(K,1)=IAI
      MD(IAI,2)=MD(IAI,2)+1
 1310 CONTINUE
 1320 CONTINUE
c ricerca del sottogruppo minimo
      JJ=10000
      do 1330 i=1,IAI
      if(MD(I,2).ge.JJ.or.MD(I,2).le.2)go to 1330
      JJ=MD(I,2)
 1330 continue
C fra i sottogruppi minimi viene scelto il piu' distante
C dall'asse unico
      dm=0.
      do 1350 I=1,IAI
      if(MD(I,2).ne.JJ)go to 1350
      do 1340 k=1,N
      if(MD(k,1).ne.I)go to 1340
        if(D(k).gt.dm)then
        dm=D(k)
        II=I
      endif
 1340 continue
 1350 continue
c     write(OUT,*)'al 1350 II,JJ',II,JJ
c     write(OUT,'(a5)')'MD(1)'
c     write(OUT,'(25i3)'),(MD(i,1),i=1,N)
c     write(OUT,'(a5,25i3)')'MD(2)',(MD(i,2),i=1,IAI)
      call azzera(AA,0.d0,9)
      N1=N-1
      DO 1370 I1=1,N1
      IF(MD(I1,1).NE.II)GO TO 1370
      IN1=I1+1
      DO 1360 I2=IN1,N
      IF(MD(I2,1).NE.II)GO TO 1360
      call combli(XO,XO,CO,1.d0,1.d0,I1,I2,1)    
      call prodv(CO,CO,CO,3,1,2)
      cam=CO(1,2)*CO(1,2)+CO(2,2)*CO(2,2)+CO(3,2)*CO(3,2)
      if(cam.le.0.00001)go to 1360
      call norm(CO,2)
c     write(OUT,'(4f10.5)')CO(1,2),CO(2,2),CO(3,2)
      call prodv(CO,CO,CO,2,3,1)
      call traspo(CO,CO,1,1)
      call r_frame(CO,XO,N)
      com=1.d0
 1355 AA(1,1)=-1
      AA(2,2)=-1
      AA(3,3)= com
      call verify(XO,AA,MK,MN,MV,N)
c     write(OUT,*)' 1350 tentativo con la matrice:'
c     write(OUT,'(3f9.5)')((AA(kkk,lll),lll=1,3),kkk=1,3)
c     write(OUT,*)'NMS,MV',NMS,MV
      AA(1,1)=-1
      AA(2,2)= 1
      AA(3,3)= com
      call verify(XO,AA,MK,MN,MV,N)
c     write(OUT,*)' 1350 tentativo con la matrice:'
c     write(OUT,'(3f9.5)')((AA(kkk,lll),lll=1,3),kkk=1,3)
c     write(OUT,*)'NMS,MV',NMS,MV
      AA(1,1)= 1
      AA(2,2)=-1
      AA(3,3)= com
      call verify(XO,AA,MK,MN,MV,N)
c     write(OUT,*)' 1350 tentativo con la matrice:'
c     write(OUT,'(3f9.5)')((AA(kkk,lll),lll=1,3),kkk=1,3)
c     write(OUT,*)'NMS,MV',NMS,MV
      AA(1,1)= 0
      AA(1,2)= 1
      AA(2,1)= 1
      AA(2,2)= 0
      AA(3,3)= com
      call verify(XO,AA,MK,MN,MV,N)
c     write(OUT,*)' 1350 tentativo con la matrice:'
c     write(OUT,'(3f9.5)')((AA(kkk,lll),lll=1,3),kkk=1,3)
c     write(OUT,*)'NMS,MV',NMS,MV
      if(NMS.ge.3)go to 1110
      if(com.le.0.0)go to 1360
      com=-com
      call azzera(AA,0.d0,9)
      go to 1355
cdi ognuna di queste quattro  matrici fare l'equivalente con -z!!!!!
 1360 continue
 1370 continue
      go to 2500
C 
CASO TRIPLAMENTE DEGENERE
 2000 N1=N-1
      N2=N-2
      DO 2010 I=1,N
      MD(I,1)=0
      MD(I,2)=0
      DO 2010 J=2,NMG
 2010 MK(I,J)=0
      IAI=0
C IAI=NUMERO DI ENNUPLE DI ATOMI AVENTI DISTANZE COMPRESE IN UN RANGE
C         X< D <X+DXM*.5
C MD(I,1)=NUMERO DELLA ENNUPLA DI CUI L'ATOMO I FA PARTE
C MD(IAI,2)=N PER IL PRIMO ATOMO DELLA ENNUPLA DOVE N E' IL NUMERO DI
C         ATOMI DELLA STESSA
      DO 2050 I=1,N1
      IF(MD(I,1).NE.0)GO TO 2050
      IAI=IAI+1
      MD(I,1)=IAI
      MD(IAI,2)=1
      J=I+1
      DO 2040 K=J,N
      IF(MN(I).NE.MN(K).or.K.eq.I.or.MD(K,1).ne.0)GO TO 2040
      IF(ABS(D(I)-D(K)).GT.0.5d0*DXM(I)) GO TO 2040
      MD(K,1)=IAI
      MD(IAI,2)=MD(IAI,2)+1
 2040 CONTINUE
 2050 CONTINUE
c     write(OUT,'(A5)')'MD(1)'
c     write(OUT,'(20i4)')(MD(kkk,1),kkk=1,N)
c     write(OUT,'(a5,/,20i4)')'MD(2)',(MD(kkk,2),kkk=1,IAI)
C SE CI SONO ATOMI SU POSIZIONI SPECIALI USA QUESTI PER TROVARE LA
C MATRICE DI ROTAZIONE
 2100 call azzera(CO,0.d0,9)
      mmd=1000
      do i=1,IAI
        if(MD(I,2).lt.mmd.and.MD(I,2).gt.1)then
          mmd=MD(I,2)
          II=i
        endif
      enddo
c     write(OUT,*)'II',II
c perche' ci sia simmetria I o Ih ci devono essere almeno 12 atomi equivalenti
      if(MD(II,2).eq.12.or.MD(II,2).eq.20.or.MD(II,2).eq.30.or.MD(II,2)
     *.gt.48)call icosahed(XO,PESO,N,MN,MK,MD,II,MDEG,*2780)
      if(MDEG.eq.1)go to 1020
      if(NMS.ge.5)then
        MDEG=1
        call azzera(AA,0.d0,9)
        write(out,'(i3,3f9.5)')(I,(XO(k,I),k=1,3),i=1,N)
        AA(1,1)= 1
        AA(2,2)= 1
        AA(3,3)=-1
        call verify(XO,AA,MK,MN,MV,N)
        AA(1,1)= 1
        AA(2,2)=-1
        AA(3,3)= 1
        call verify(XO,AA,MK,MN,MV,N)
        AA(1,1)=-1
        AA(2,2)= 1
        AA(3,3)= 1
        call verify(XO,AA,MK,MN,MV,N)
        go to 2500
      endif
C
C RICERCA ASSE 3 SU ATOMI IN POSIZIONE GENERALE
C NOTA: TRE ATOMI EQUIVALANTI PER UN ASSE 3 DEVONO FORMARE UN TRIANGOLO
C EQUILATERO
C II E' IL NUMERO DI ATOMI DEL SOTTOGRUPPO
C
 2110 CA=0.d0
C matrice di rotazione 3 coincidente con z
      call azzera(AA,0.d0,9)
      AA(3,3)=1.d0
      AA(1,1)=-.5d0
      AA(2,2)=-.5d0
      AA(2,1)=0.5d0*DSQRT(3.d0)
      AA(1,2)=-AA(2,1)
      NMS1=1
C
C     I1,I2,I3=INDICATORI DEI TRE ATOMI POSSIBILE GENERATORI DELL'ASSE 3
C
      N2=N-2
      N1=N-1
      DO 2200 I1=1,N2
      IF(MD(I1,1).NE.II)GO TO 2200
      IN1=I1+1
      DO 2190 I2=IN1,N1
      IF(MD(I2,1).NE.II)GO TO 2190
      call combli(XO,XO,CO,1.d0,-1.d0,I1,I2,4)
      call prods(CO,CO,CA,4,4,2)
      IN2=I2+1
      DO 2180 I3=IN2,N
      IF(MD(I3,1).NE.II)GO TO 2180
      call combli(XO,XO,CO,1.d0,-1.d0,I1,I3,5)
      call prods(CO,CO,CB,5,5,2)
      call combli(XO,XO,CO,1.d0,-1.d0,I2,I3,6)
      call prods(CO,CO,CC,6,6,2)
      IF(DABS(CA-CB).GT.DXM(I1))GO TO 2180
      IF(DABS(CA-CC).GT.DXM(I1))GO TO 2180
      IF(DABS(CB-CC).GT.DXM(I1))GO TO 2180
      CO(1,3)=XO(1,I1)+XO(1,I2)+XO(1,I3)
      CO(2,3)=XO(2,I1)+XO(2,I2)+XO(2,I3)
      CO(3,3)=XO(3,I1)+XO(3,I2)+XO(3,I3)
      call norm(CO,3)
c     write(out,'(a15,3i4)')'I1,I2,I3',I1,I2,I3
c     write(out,'(a5,3f9.5)')'asse',(CO(kkk,3),kkk=1,3)
C salvataggio del possibile asse 3
      CO(1,6+NMS1)=CO(1,3)
      CO(2,6+NMS1)=CO(2,3)
      CO(3,6+NMS1)=CO(3,3)
Controllo che l'asse trovato non coincida gia' con z cioe'
con il primo asse C3 trovato
      if(NMS1.gt.1)then
        call prods(CO,CO,pro,7,8,1)
        if(DABS(pro).gt.0.9d0)go to 2180
        if(pro.gt.0.d0)then
          CO(1,8)=-CO(1,8)
          CO(2,8)=-CO(2,8)
          CO(3,8)=-CO(3,8)
        endif
      endif
C definizione di un secondo asse
      do i=1,3
        call azzera(CO,0.d0,3)
        CO(I,1)=1.d0
        call prods(CO,CO,pro,1,3,1)
        if(DABS(pro).gt.0.5)go to 2150
      enddo
C verifica se l' asse trovato e' un C3
C rotazione per portare il C3 a coincidere con z
 2150 call prodv(CO,CO,CO,3,1,2)
      call norm(CO,2)
      call prodv(CO,CO,CO,2,3,1)
      call traspo(CO,CO,1,1)
      call r_frame(CO,XO,N)
      call verify(XO,AA,MK,MN,MV,N)
      NMS=1
      NMS1=NMS1+MV
c     write(OUT,*)' 2150 tentativo con la matrice:'
c     write(OUT,'(3f9.5)')((AA(kkk,lll),lll=1,3),kkk=1,3)
c     write(OUT,*)'NMS1,MV',NMS1,MV
      call traspo(CO,CO,1,1)
      call r_frame(CO,XO,N)
      if(NMS1.eq.3)GO TO 2350
 2180 CONTINUE
 2190 CONTINUE
 2200 CONTINUE
      write(out,2)
    2 format('********************** WARNING **********************',//,
     *'       INCREASING THE TOLERANCE COULD BE USEFUL',//,
     *'*****************************************************')
c     write(OUT,*)'NMS1,MV',NMS1,MV
      if(NMS1.eq.1)then
C se ha gia' modificato una volta i pesi non li modifica ulteriormente
C ma riduce la tolleranza accettata fra i momenti di inerzia per il
calcolo della degenerazione
        if(porig.ne.PESO(1))then
          deliner=deliner*0.1d0
          go to 460
        endif
C 
C pseudodegenerazione 3 senza assi C3. Esiste ancora la possibilita'
C che la pseudodegenerazione sia completa (MDEG=0) o che ci sia una
C asse 4,-4,5,-5,7,-7,8,-8
C 
        write(out,*)'Weights are changed'
        write(out,*)
        do k=1,N
          PESO(K)=PESO(K)*(D(k)/DM)**4
c         write(out,'(i3,g10.3)')k,PESO(K)
        enddo
        go to 460
      endif
C ha trovato un solo asse 3 [ = CO(7)] e si riporta in quel
C riferimento ripristinando NMS
C definizione di un secondo asse
      CO(1,3)=CO(1,7)
      CO(2,3)=CO(2,7)
      CO(3,3)=CO(3,7)
      do i=1,3
        call azzera(CO,0.d0,3)
        CO(I,1)=1.d0
        call prods(CO,CO,pro,1,3,1)
        if(DABS(pro).gt.0.5)go to 2310
      enddo
 2310 call prodv(CO,CO,CO,3,1,2)
      call norm(CO,2)
      call prodv(CO,CO,CO,2,3,1)
      call traspo(CO,CO,1,1)
      call r_frame(CO,XO,N)
      go to 2490
C ha trovato due assi C3 [CO(7) e CO(8)] e li usa per
C determinare il riferimento definitivo
 2350 call combli(CO,CO,CO,1.d0,1.d0,7,8,1)
      call norm(CO,1)
      call combli(CO,CO,CO,1.d0,-1.d0,7,8,3)
      call prodv(CO,CO,CO,3,1,2)
      call norm(CO,2)
      call prodv(CO,CO,CO,1,2,3)
      call traspo(CO,CO,1,1)
      call r_frame(CO,XO,N)
      call azzera(CO,0.d0,9)
      CO(1,1)=1.d0
      CO(2,2)=DSQRT(0.5d0)
      CO(3,3)=CO(2,2)
      CO(2,3)=CO(2,2)
      CO(3,2)=-CO(2,2)
c     write(out,*)'Nuove coordinate'
c     write(out,'(i2,3f9.5)')(lll,(XO(kkk,lll),kkk=1,3),lll=1,N)
      go to 2410
C ORTONORMALIZZAZIONE DELLA MATRICE DI ROTAZIONE
 2400 call norm(CO,5)
c     write(out,*)'MD'
c     write(out,'(30i4)')(MD(i,1),i=1,N)
c     write(out,'(30i4)')(MD(i,2),i=1,N)
      call prodv(CO,CO,CO,5,6,7)
      call norm(CO,7)
      call prodv(CO,CO,CO,7,5,6)
      call trasfv(CO,CO,5,1)
      call trasfv(CO,CO,6,2)
      call trasfv(CO,CO,7,3)
      call traspo(CO,CO,1,1)
C ROTAZIONE SU NUOVA TERNA
c     write(out,*)' CO prima di r_frame'
c     write(out,'(3f10.5)')((CO(kk,ii),ii=1,3),kk=1,3)
 2410 call r_frame(CO,XO,N)
C VERIFICA GENERALE ESISTENZA ASSE TERNARIO
c     write(out,*)' 2410 coordinate'
c     write(out,'(i2,3f9.5)')(lll,(XO(kkk,lll),kkk=1,3),lll=1,N)
      call azzera(AA,0.d0,9)
      AA(1,2)=1.d0
      AA(2,3)=1.d0
      AA(3,1)=1.d0
      call verify(XO,AA,MK,MN,MV,N)
c     write(out,*)' 2400 tentativo con la matrice:'
c     write(out,'(3f9.5)')((CO(kkk,lll),lll=1,3),kkk=1,3)
c     write(out,*)'NMS,MV',NMS,MV
      NASS=1
C ha torvato un asse 3 (diagonale) ora cerca un asse 2
      IU=1
      IB=2
      IC=3
      MORD=4
      call azzera(CO,0.d0,9)
      CO(1,1)=-1.d0
      CO(2,2)=-1.d0
      CO(3,3)=1.d0
      call verify(XO,CO,MK,MN,MV,N)
c     write(out,*)' 2400 tentativo con la matrice:'
c     write(out,'(3f9.5)')((CO(kkk,lll),lll=1,3),kkk=1,3)
c     write(out,*)'NMS,MV',NMS,MV
      if(MV.eq.1)GO TO 1120
C non esiste l'asse 2 allineato su x: si tratta di una pseudodegenerazione
C mi metto in modo che l'asse 3 sia allineato a z
      call azzera(CO,0.d0,9)
      call azzera(AA,0.d0,9)
      cost=dsqrt(0.5d0)
      CO(1,1)=cost
      CO(2,2)=cost
      CO(1,2)=-cost
      CO(2,1)= cost
      CO(3,3)= 1.d0
      ang=0.5d0*DACOS(-1.d0/3.d0)
      aaa=ang*180.d0/PIG
      cosa=DCOS(ang)
      sina=DSIN(ang)
      AA(2,2)=cosa
      AA(3,3)=cosa
      AA(2,3)=-sina
      AA(3,2)=sina
      AA(1,1)=1.d0
      call prodmm(AA,CO,CO,1,1,1)
      call r_frame(CO,XO,N)
 2490 NMS=1
      MORD=6
      MDEG=1
      IU=3
      IB=2
      IC=1
      write(6,5)
      if(out.eq.9)write(out,5)
    5 format(
     *'***********************************************************',/,
     *'WARNING: the degeneration degree is 3 but no cubic',/,
     *'         or icosahedral group can be found.',/,
     *'IF YOU SUSPECT THE EXISTENCE OF ONE OF THEM, PLEASE CHANGE:',/,
     *'1) the weighting scheme   OR',/,
     *'2) put MOL<0 to the atoms farest from the baricenter   OR',/,
     *'3) enlarge DCM',/,
     *'***********************************************************')
      go to 1120
C RICERCA DEL CENTRO DI SIMMETRIA
 2500 II=2
      call azzera(CO,0.d0,9)
      CO(1,1)=-1.d0
      CO(2,2)=-1.d0
      CO(3,3)=-1.d0
      call verify(XO,CO,MK,MN,MV,N)
c     write(out,*)' 2500 tentativo con la matrice:'
c     write(out,'(3f9.5)')((CO(kkk,lll),lll=1,3),kkk=1,3)
c     write(out,*)'NMS,MV',NMS,MV
      IF(NMS.NE.1)GO TO 2520
      write(out,3)
      RETURN 
complete the group
 2520 call compl_gr(MK,N,*4100)
C_______________________________________________________________________
C gruppi lineari
      if(MDEG.eq.5)go to 2790
C_______________________________________________________________________
      ntest=0
      do i=1,NMS
        if(MTG(i,i).ne.1)ntest=1
      enddo
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c variazione per far si che in un C2 o C2v l'asse 2 coincida con z
      if(NMS.eq.2.or.NMS.eq.4)then
        nax2=0
        lax2=0
        do i=2,NMS
          itrac=NINT(SIM(1,1,i)+SIM(2,2,i)+SIM(3,3,i))
          iabst=NINT(DABS(SIM(1,1,i))+DABS(SIM(2,2,i))+DABS(SIM(3,3,i)))
          if(iabst.ne.3)go to 2740
          if(itrac.eq.-1)then
            nax2=nax2+1
            lax2=i
          endif
        enddo
        if(nax2.eq.0.or.nax2.eq.3)go to 2780
        if(nax2.eq.1)then
          do l=1,3
            if(NINT(SIM(l,l,lax2)).eq.1)go to 2730
          enddo
 2730     if(l.eq.3)go to 2780
          call azzera(CO,0.d0,9)
          i=3-l
          k=6-i-l
          CO(k,l)=1.d0
          CO(l,k)=-1.d0
          CO(i,i)=1.d0
          call r_frame(CO,XO,N)
          com=RIN(3)
          RIN(3)=RIN(l)
          RIN(l)=com
        endif
        NMS=1
        go to 1020
      endif
 2740 if(ntest.eq.0)go to 2780
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Ottimizzazione del riferimento per MDEG>0
 2780 if(MDEG.eq.0)go to 2790
      do i=1,NMS
       ISU(i)=i
      enddo
      call s_coor(SIM,XO,XS,CO,D,DEL,DELM,CSMT,MK,IDM,ISU,NMS,N)
      npb=0
      do k=1,NMS
       do i=1,N
        npu=npb+i
        call trasfv(XO,ppu,i,npu)
        npu=npb+MK(i,k)
        call prodmv(SIM,XS,ppo,k,i,npu)
       enddo
       npb=npb+N
      enddo
      npu=npb
      call rms_min(V)
      arms=crms(V)
c     write(6,*)'v,RV'
c     write(6,'(3f10.5)')V,RV
      call r_frame(RV,XO,N)
Calcolo simmetria per sottogruppi di operazioni
 2790 CSM(1)=0.d0
      DEV(1)=0.d0
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c     write(out,*)'GROUP MULTIPLICATION TABLE'
c     write(out,*)
c     nst=NMS
c     nis=25
c     if(NMS.gt.24)nst=24
c     do i=1,NMS
c      write(out,'(24i3)')(MTG(i,j),j=1,nst)
c     enddo
c     write(out,*)
c     if(nst.ne.NMS)then
c      do i=1,NMS
c       write(out,'(24i3)')(MTG(i,j),j=nis,NMS)
c      enddo
c     endif
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c simmetrizzazione per sottogruppi di una solo elemento (e sue potenze)
      do 2950 I=2,NMS
      nes=1
      ISU(1)=nes
      k=i
 2800 nes=nes+1
      ISU(nes)=k
      l=MTG(i,k)
      if(l.eq.1)go to 2900
      k=l
      go to 2800
 2900 call s_coor(SIM,XO,XS,CO,D,DEL,DELM,comt,MK,IDM,ISU,nes,N)
      DEV(I)=DELM(4)
      CSM(I)=comt
 2950 continue
CALCOLA LE COORDINATEE SIMMETRIZZATE PER L'INTERO GRUPPO
 3000 do i=1,NMS
      ISU(i)=i
      enddo
      call s_coor(SIM,XO,XS,CO,D,DEL,DELM,CSMT,MK,IDM,ISU,NMS,N)
c     write(out,*)
c     write(out,*)'FRACTIONAL COORDINATES OF THE CENTRE OF MASS'
c     write(out,'(3f12.6,/)')BARC
c     write(out,*)'ORTHOGONAL COORDINATES OF THE CENTRE OF MASS'
c     write(out,'(3f12.6,/)')BARO
c     write(out,*)'PRINCIPAL INERTIA MOMENTS'
c     write(out,'(3f10.2,/)')RIN
      write(out,*)
      write(out,*)'   ORTHOGONALIZATION MATRIX'
      write(out,'(3f12.6)')((OT(k,i),i=1,3),k=1,3)
      write(out,*)
      write(out,*)'ATOM      ORTHOGONAL COORDINATES      VECTORS TO MAKE
     * SYMMETRICAL THE GROUP'
c     do I=1,N
c     write(out,'(24I4)')(MK(I,K),K=1,NMS)
c     enddo
C Recupero degli atomi con MOL<0
      N1=N
      do 3100 i=1,NA
      if(MOL(i).ne.-MO)go to 3100
      N1=N1+1
      IEQAT(N1)=i
      MN(N1)=MSP(i)
      MK(N1,1)=-N1
      call combli(X,BARC,XO,1.d0,-1.d0,i,1,N1)
      call prodmv(OT,XO,XO,1,N1,N1)
      XS(1,N1)=XO(1,N1)
      XS(2,N1)=XO(2,N1)
      XS(3,N1)=XO(3,N1)
 3100 continue
      NSTART=N+1
      IF(N1.eq.N)go to 3400
      do 3300 I=NSTART,N1
      do 3200 L=2,NMS
      MK(i,l)=0
      FMIN=10000.
      call prodmv(SIM,XO,V,L,I,1)
      DO 3150 K=NSTART,N1
      IF(MN(I).NE.MN(K))GO TO 3150
      call combli(XO,V,W,1.d0,-1.d0,K,1,1)
      call prods(W,W,F,1,1,2)
      if(F.le.FMIN)then
        KK=K
        FMIN=F
      endif
 3150 continue
      if(FMIN.LT.2.*DCME)then
      MK(I,L)=KK
      else
      MK(I,1)=0
      endif
 3200 continue
 3300 continue
c     do I=1,N1
c     write(out,'(24I4)')(MK(I,K),K=1,NMS)
c     enddo
      do 3350 i=NSTART,N1
      if(MK(I,1).eq.0)go to 3350
      do k=2,NMS
        call prodmv(SIM,XO,V,k,i,1)
        l=MK(i,k)
        call combli(XS,V,XS,1.d0,1.d0,l,1,l)
      enddo
 3350 continue
      do 3360 I=NSTART,N1
      if(MK(I,1).eq.0)go to 3360
        do k=1,3
          XS(k,I)=XS(k,I)/dfloat(NMS)
        enddo
 3360 continue
 3400 write(out,*)
      call asymunit(MK,IASU,N1,NMS)
      do i=1,N1
        call combli(XS,XO,delta,1.d0,-1.d0,i,i,i)
        write(out,'(a6,2x,2(3f10.5,5x))')NAME(IEQAT(I)),
     *  (XO(k,I),k=1,3),(delta(k,i),k=1,3)
      enddo
      write(out,*)
      write(out,'(5x,a34,12x,a)')'SYMMETRIZED ORTHOGONAL COORDINATES',
     *'ATOMIC R.M.S.'
      ind1=0
      ind2=0
      do I=1,N1
      comment=' '
      car=' '
      if(MK(i,1))3410,3430,3420
 3410 comment='#'
      ind1=1
 3420 car=type(IASU(I))
      go to 3440
 3430 comment='$'
      ind2=2
 3440 if(MK(i,1).gt.0)then
        write(out,'(a6,i2,2(3f9.5,2x),a1)')NAME(IEQAT(I)),
     *  MOL(IEQAT(I)),(XS(k,I),k=1,3),(RMS(k,I),k=1,3),car
      else
        write(out,'(a6,i2,3f9.5,29x,a1,1x,a1)')NAME(IEQAT(I)),
     *  MOL(IEQAT(I)),(XS(k,I),k=1,3),comment,car
      endif
C IASU(I)=2 se l'atomo rappresenta l'unita' asimmetrica altrimenti IASU(I)=1
      enddo
      ind1=ind1+ind2
      write(out,*)
      write(out,*)'* Atom defining the asymmetric unit for the found sym
     *metry group.'
      if(ind1.eq.0)go to 3460
      if(ind1.eq.2)go to 3450
      write(out,*)
      write(out,*)'# This atom was symmetrized but NOT used to find the
     *symmetry group and to calculate CMS, RMS and so on.'
 3450 if(ind1.eq.1)go to 3460
       write(out,*)
       write(out,*)'$ It was IMPOSSIBLE to symmetrize this atom accordin
     *g to the found symmetry group and within the given tolerance.'
 3460 write(out,*)
      do i=1,4
        DEL(i)=DEL(I)/dfloat(N)
        id1=idm(i)
        if(IDM(i).ne.0) then
        ie1=IEQAT(id1)
          NNAME(i)=NAME(ie1)
        else
          NNAME(i)=' --- '
        endif
      enddo
      write(out,*)
      write(out,'(1x,a29,4F10.5,/)')'AVERAGE DIFFERENCE ON X,Y,Z,D',DEL
      write(out,'(1x,a29,4F10.5)')'MAXIMUM DIFFERENCE ON X,Y,Z,D',DELM
      write(out,'(1x,a29,4a10,/)')'DUE TO THE ATOMS',
     *(NNAME(i),i=1,4)
      if(out.eq.9)then
       write(*,*)
       write(*,'(1x,a29,4F10.5,/)')'AVERAGE DIFFERENCE ON X,Y,Z,D',DEL
       write(*,'(1x,a29,4F10.5)')'MAXIMUM DIFFERENCE ON X,Y,Z,D',DELM
       write(*,'(1x,a29,4a10,/)')'DUE TO THE ATOMS',
     * (NNAME(i),i=1,4)
      endif
C
      do i=1,N 
      do j=1,NMS
      K=MK(i,j)
      call prodmv(SIM,XS,V,j,i,1)
      do l=1,3
      diff=abs(XS(l,k)-V(l))
      if(diff.gt.0.001d0)then
       write(out,*)
       write(out,*)
     * '                    ========> CAUTION!!! <=======' 
       write(out,*)
       write(out,*)  ' THE SYMMETRY OPERATIONS ARE NOT CONSISTENT WITH THE
     * EQUIVALENTS POINTS.'
       write(out,'(a27,f10.5)')' DIFFERENCE FOUND =',diff
       write(out,*)
     * '     ========> TRY TO REDUCE THE CONSTANT OF TOLERANCE <=======' 
c      write(out,'(3(3f11.7,2x))')(XS(kk,i),kk=1,3),V,(XS(kk,k),kk=1,3)
       go to 3500
c      RETURN
      endif
      enddo
      enddo
      enddo
C 
 3500 write(out,*)'Bond lengths and bond angles after symmetrization'
      call bond(XS,IASU,N1)
      call out_bond(XS,IEQAT,IASU,N)
      if(PC(7).gt.2.d0)then
      write(out,'(7x,a)')' SYMMETRIZED FRACTIONAL COORDINATES'
      do 3600 i=1,N1
      k=IEQAT(i)
      call prodmv(OTI,XS,X,1,i,k)
      call combli(X,BARC,X,1.d0,1.d0,k,1,k)
      write(out,'(1x,a6,i2,6f9.5)')NAME(k),MOL(k),(X(j,k),j=1,3),
     *(SX(j,k),j=1,3)
 3600 continue
      endif
      call schoenfl(MDEG)
      if(IPOGROU.ne.0.and.(IPOGROU.lt.16.or.IPOGROU.gt.27))then
      call stasy
      else
      write(out,6)
      do 3800 k=1,NMS
      write(out,8)k,CSM(k),DEV(k),lbls(K)
    8 format(/,i3,' CSM =',f7.2,5x,'MAX. DIFF. (Angstrom)=',F6.4,
     *5x,'TYPE ',a3)
      write(out,7)((SIM(I,J,K),J=1,3),I=1,3)
 3800 continue
      endif
      IF(IES.NE.1)RETURN
      if(IPOGROU.le.0.or.IPOGROU.gt.32)RETURN
      AL(2,2)=SQRT(3.D0)*0.5D0
      AL(1,2)=-0.5D0
      AL(1,1)=1.d0
      AL(3,3)=1.d0
      call inv33(AL,ALX,1,1)
      write(out,*)' SYMMETRY OPERATIONS IN HEXAGONAL COORDINATES'
      DO 3900 I=1,NMS
      call prodmm(ALX,SIM,AA,1,I,1)
      call prodmm(AA,AL,SIM,1,1,I)
 3900 continue
      call stasy
      write(out,*)
      write(out,11)
   11 format(' OBLIQUE COORDINATES (HEXAGONAL SYSTEM)')
      do 4000 I=1,N1
      call prodmv(ALX,XS,CO,1,I,I)
      write(out,'(1x,a8,3f10.5)')NAME(IEQAT(I)),(CO(J,I),J=1,3)
 4000 continue
 4100 RETURN 
    3 format(' NO SYMMETRY EXISTS WITHIN GIVEN TOLERANCE')
    6 format(/,' SYMMETRY GROUP MATRICES',/)
    7 format(3(1X,3F16.10,/))
      END
