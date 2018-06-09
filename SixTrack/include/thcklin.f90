#ifdef DEBUG
!     if (i.ge.673.and.j.le.2) then
!       call warr('bl30 i j',xv(1,2),i,j,0,0)
!       call warr('xv12,i,ktrack ',xv(1,2),i,ktrack(i),0,0)
!     endif
#endif
  puxve1=xv(1,j)
  puzve1=yv(1,j)
  puxve2=xv(2,j)
  puzve2=yv(2,j)
  sigmv(j)=(((((sigmv(j)+as(1,1,j,jx))+puxve1*((as(2,1,j,jx)+ as &!hr03
  (4,1,j,jx)*puzve1)+as(5,1,j,jx)*puxve1))+ puzve1*(as           &!hr03
  (3,1,j,jx)+as(6,1,j,jx)*puzve1))                               &!hr03
  +as(1,2,j,jx))+puxve2*(as(2,2,j,jx)+ as                        &!hr03
  (4,2,j,jx)*puzve2+as(5,2,j,jx)*puxve2))+ puzve2*(as            &!hr03
  (3,2,j,jx)+as(6,2,j,jx)*puzve2)                                 !hr03
    xv(1,j)=(al(1,1,j,jx)*puxve1+ al(2,1,j,jx)*puzve1)+          &
  real(idz1,fPrec)*al(5,1,j,jx)                                   !hr03
    xv(2,j)=(al(1,2,j,jx)*puxve2+ al(2,2,j,jx)*puzve2)+          &
  real(idz2,fPrec)*al(5,2,j,jx)                                   !hr03
    yv(1,j)=(al(3,1,j,jx)*puxve1+ al(4,1,j,jx)*puzve1)+          &
  real(idz1,fPrec)*al(6,1,j,jx)                                   !hr03
    yv(2,j)=(al(3,2,j,jx)*puxve2+ al(4,2,j,jx)*puzve2)+          &
  real(idz2,fPrec)*al(6,2,j,jx)                                   !hr03
#ifdef DEBUG
!     if (i.ge.673.and.j.le.2) then
!       call warr('al30 i j',xv(1,2),i,j,0,0)
!       call warr('xv12,i,ktrack ',xv(1,2),i,ktrack(i),0,0)
!     endif
#endif
  ! A.Mereghetti and P.Garcia Ortega, for the FLUKA Team
  ! last modified: 07-03-2018
  ! store old particle coordinates
  if (lbacktracking) call aperture_saveLastCoordinates(i,ix,kz(jx))