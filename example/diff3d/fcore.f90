!
!   core.f90
!   diff3dhp
!
!   Created by keno on 2016/06/26.
!   Copyright 2016年 keno. All rights reserved.
!

!> ********************************************************************
!! @brief Initialize arrays
!! @param [in]     sz   array length
!! @param [in]     g    guide cell
!! @param [in,out] q,w  variables
!<
subroutine initialize (sz, g, q, w)
implicit none
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  q, w

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL FIRSTPRIVATE(ix, jx, kx, g)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-g, kx+g
do j=1-g, jx+g
do i=1-g, ix+g
  q(i,j,k) = 0.0
  w(i,j,k) = 0.0
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine initialize


!> ************************************
!! @brief Euler explicit scheme
!! @param [in]     sz      array length
!! @param [in]     g       guide cell
!! @param [in,out] q       variable
!! @param [in]     w       temporary work array
!! @param [in]     dh      mesh width
!! @param [in]     dt      time increment
!! @param [in]     alpha   coefficient
!! @param [out]    res     residual
!<
subroutine euler_explicit (sz, g, q, w, dh, dt, alpha, res)
implicit none
integer                                                :: i, j, k, ix, jx, kx, g
integer, dimension(3)                                  :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: q, w
real                                                   :: dh, dt, alpha, c, res, delta, q0

ix = sz(1)
jx = sz(2)
kx = sz(3)

c = alpha*dt/(dh*dh)
res = 0.0

!$OMP PARALLEL FIRSTPRIVATE(ix, jx, kx, c) &
!$OMP PRIVATE(q0, delta) &
!$OMP REDUCTION(+:res)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1, kx
do j=1, jx
do i=1, ix
  q0 = q(i,j,k)
  delta = c * ( q(i-1,j,k) + q(i+1,j,k) &
              + q(i,j-1,k) + q(i,j+1,k) &
              + q(i,j,k-1) + q(i,j,k+1) &
              - 6.0 * q0 )
  w(i,j,k) = q0 + delta
  res = res + delta * delta
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL FIRSTPRIVATE(ix, jx, kx)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1, kx
do j=1, jx
do i=1, ix
  q(i,j,k) = w(i,j,k)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine euler_explicit


!> ************************************
!! @brief Boundary condition
!! @param [in]     sz   array length
!! @param [in]     g    guide cell
!! @param [in,out] p    variable
!! @param [in]     dh   mesh width
!! @param [in]     org  origin of subdomain
!! @param [in]     nID  table
!<
subroutine bc (sz, g, p, dh, org, nID)
implicit none
integer                                                :: i, j, k, ix, jx, kx, g
integer, dimension(3)                                  :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: p
real                                                   :: pi, x, y, dh, ox, oy, oz
real, dimension(3)                                     :: org
integer, dimension(6)                                  :: nID

ix = sz(1)
jx = sz(2)
kx = sz(3)

pi = 2.0*asin(1.0)

ox = org(1)
oy = org(2)
oz = org(3)

!$OMP PARALLEL FIRSTPRIVATE(ix, jx, kx, g, pi, dh, ox, oy)

! ZMINUS Dirichlet
if ( nID(5).lt.0 ) then

!$OMP DO SCHEDULE(static) PRIVATE(i, j, x, y)
do j=1-g,jx+g
do i=1-g,ix+g
  x = ox + dh*(real(i-1)+0.5)
  y = oy + dh*(real(j-1)+0.5)
  p(i,j,0) = sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO

end if


! ZPLUS Dirichlet
if ( nID(6).lt.0 ) then

!$OMP DO SCHEDULE(static) PRIVATE(i, j, x, y)
do j=1-g,jx+g
do i=1-g,ix+g
  x = ox + dh*(real(i-1)+0.5)
  y = oy + dh*(real(j-1)+0.5)
  p(i,j,kx+1) = sin(pi*x)*sin(pi*y)
end do
end do
!$OMP END DO

end if


! XMINUS
if ( nID(1).lt.0 ) then

!$OMP DO SCHEDULE(static)
do k=1,kx
do j=1,jx
  p(0,   j,k) = -p(1, j,k)
end do
end do
!$OMP END DO

end if


! XPLUS
if ( nID(2).lt.0 ) then

!$OMP DO SCHEDULE(static)
do k=1,kx
do j=1,jx
  p(ix+1,j,k) = -p(ix,j,k)
end do
end do
!$OMP END DO

end if


! YMINUS
if ( nID(3).lt.0 ) then

!$OMP DO SCHEDULE(static)
do k=1,kx
do i=1,ix
  p(i,0,   k) = -p(i,1, k)
end do
end do
!$OMP END DO

end if


! YPLUS
if ( nID(4).lt.0 ) then

!$OMP DO SCHEDULE(static)
do k=1,kx
do i=1,ix
  p(i,jx+1,k) = -p(i,jx,k)
end do
end do
!$OMP END DO

end if

!$OMP END PARALLEL

return
end subroutine bc


!> ************************************
!! @brief fileout by sph format
!! @param [in]     sz     array length
!! @param [in]     g      guide cell
!! @param [in]     stp    step
!! @param [in]     tm     time
!! @param [in]     dh     mesh width
!! @param [in]     org    origin
!! @param [in]     fname  file name
!! @param [in]     s      variable
!<
subroutine write_sph (sz, g, stp, tm, dh, org, fname, s)
implicit none
integer                                                  :: ix, jx, kx, i, j, k, g, stp
integer, dimension(3)                                    :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   :: s
real                                                     :: dh, tm, ox, oy, oz
real,dimension(3)                                        :: org
character*20                                             :: fname

ix = sz(1)
jx = sz(2)
kx = sz(3)

ox = org(1)-dh
oy = org(2)-dh
oz = org(3)-dh

open (unit=22, file=fname, form='unformatted')
write (22) 1, 1
write (22) ix+2*g, jx+2*g, kx+2*g
write (22) ox, oy, oz
write (22) dh, dh, dh
write (22) stp, tm
write (22) (((s(i,j,k),i=1-g,ix+g),j=1-g,jx+g),k=1-g,kx+g)
close (unit=22)

return
end subroutine write_sph
