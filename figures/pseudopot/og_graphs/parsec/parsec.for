      program readqe
      implicit none
      integer NMX
      parameter (NMX=5000)
      integer i,ic1,ic2,ic3,irows
      real*8 rr1,rr2,rr3,rr4
      real*8 r,pr,rho

      dimension r(NMX),pr(NMX),rho(NMX)

      open(unit=10,file='pp_r',status='unknown')
      open(unit=20,file='pp_local',status='unknown')
      open(unit=30,file='pp_wave',status='unknown')

      ic1=1
      irows=NMX/4
      print*,irows
      do 100 i=1,irows
        read(10,*,end=101) rr1,rr2,rr3,rr4
        r(ic1)=rr1
        r(ic1+1)=rr2
        r(ic1+2)=rr3
        r(ic1+3)=rr4
        ic1=ic1+4
100   continue
101   ic1 = ic1-1
      print*,'number of grid points',ic1

      ic2=1
      do 200 i=1,irows
        read(20,*,end=201) rr1,rr2,rr3,rr4
        pr(ic2)=rr1
        pr(ic2+1)=rr2
        pr(ic2+2)=rr3
        pr(ic2+3)=rr4
        ic2=ic2+4
200   continue
201   ic2 = ic2-1
      print*,'number of grid points',ic2

      ic3=1
      do 300 i=1,irows
        read(30,*,end=301) rr1,rr2,rr3,rr4
        rho(ic3)=rr1
        rho(ic3+1)=rr2
        rho(ic3+2)=rr3
        rho(ic3+3)=rr4
        ic3=ic3+4
300   continue
301   ic3 = ic3-1
      print*,'number of grid points',ic3


      do 400 i=1,ic3
        write(15,*) r(i),pr(i)
400   continue

      do 500 i=1,ic3
        write(25,*) r(i),rho(i)
500   continue

      stop
      end
