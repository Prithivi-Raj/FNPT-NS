program fnptMlpg3D
implicit none

  integer(kind=4)::i,j,k,l
  integer(kind=4)::nx,ny,nt,ntt,nxy
  real(kind=8)::dt,t1,t2
  real(kind=8)::data(8)

  nx=81
  ny=21
  dt=0.01d0
  t1=0d0
  t2=18.75d0
  ntt=5300
  nt=(t2-t1)/dt;
  nxy=nx*ny

  open(11,file='Output PHIT.dat')
  open(21,file='Output PHIT2.dat')

  write(*,*)'nt  = ',nt
  write(*,*)'ntt = ',ntt

  do i=1,nt
    do j=1,nxy
      read(11,*)data(1:8)
    enddo
  enddo

  do i=nt+1,ntt
    do j=1,nxy
      read(11,*)data(1:8)
      write(21,18)data(1)-t2,data(2:8)
    enddo
  enddo

  18 format(F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5)

  close(11)
  close(21)



end program fnptMlpg3D