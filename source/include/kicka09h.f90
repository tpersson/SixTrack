! start include/kicka09h.f90
#ifndef TILT
  mpe=9
  mx=7
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  ab1(7)=(28.0_fPrec*ekk)*cxzyr                                          !hr02
  ab2(7)=(-28.0_fPrec*ekk)*cxzyi                                         !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(6)=(56.0_fPrec*ekk)*cxzyr                                          !hr02
  ab2(6)=(-56.0_fPrec*ekk)*cxzyi                                         !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(5)=(70.0_fPrec*ekk)*cxzyr                                          !hr02
  ab2(5)=(-70.0_fPrec*ekk)*cxzyi                                         !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(4)=(56.0_fPrec*ekk)*cxzyr                                          !hr02
  ab2(4)=(-56.0_fPrec*ekk)*cxzyi                                         !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(3)=(28.0_fPrec*ekk)*cxzyr                                          !hr02
  ab2(3)=(-28.0_fPrec*ekk)*cxzyi                                         !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  qu=(eight*ekk)*cxzyr                                               !hr02
  qv=(eight*ekk)*cxzyi                                               !hr02
  ab2(2)=-one*qv                                                   !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  dyy1=ekk*cxzyr
  dyy2=(-one*ekk)*cxzyi                                            !hr02
  ab1(8)=(eight*ekk)*xl                                              !hr02
  ab2(8)=(-eight*ekk)*zl                                             !hr02
  ab1(9)=ekk
#else
  mpe=9
  mx=7
  cxzr=xl
  cxzi=zl
  cxzyr=cxzr**2-cxzi**2                                            !hr08
  cxzyi=cxzr*cxzi+cxzi*cxzr
  tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
  tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk1=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck1=tiltckuk
  tiltckuk=tiltck1*tiltc(k)-tiltsk1*tilts(k)
  tiltsk2=tiltck1*tilts(k)+tiltsk1*tiltc(k)
  tiltck2=tiltckuk
  tiltckuk=tiltck2*tiltc(k)-tiltsk2*tilts(k)
  tiltsk3=tiltck2*tilts(k)+tiltsk2*tiltc(k)
  tiltck3=tiltckuk
  tiltckuk=tiltck3*tiltc(k)-tiltsk3*tilts(k)
  tiltsk4=tiltck3*tilts(k)+tiltsk3*tiltc(k)
  tiltck4=tiltckuk
  tiltckuk=tiltck4*tiltc(k)-tiltsk4*tilts(k)
  tiltsk5=tiltck4*tilts(k)+tiltsk4*tiltc(k)
  tiltck5=tiltckuk
  ab1(7)=(28.0_fPrec*ekk)*(tiltck5*cxzyr+tiltsk5*cxzyi)                  !hr02
  ab2(7)=(28.0_fPrec*ekk)*(tiltsk5*cxzyr-tiltck5*cxzyi)                  !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(6)=(56.0_fPrec*ekk)*(tiltck4*cxzyr+tiltsk4*cxzyi)                  !hr02
  ab2(6)=(56.0_fPrec*ekk)*(tiltsk4*cxzyr-tiltck4*cxzyi)                  !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(5)=(70.0_fPrec*ekk)*(tiltck3*cxzyr+tiltsk3*cxzyi)                  !hr02
  ab2(5)=(70.0_fPrec*ekk)*(tiltsk3*cxzyr-tiltck3*cxzyi)                  !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(4)=(56.0_fPrec*ekk)*(tiltck2*cxzyr+tiltsk2*cxzyi)                  !hr02
  ab2(4)=(56.0_fPrec*ekk)*(tiltsk2*cxzyr-tiltck2*cxzyi)                  !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  ab1(3)=(28.0_fPrec*ekk)*(tiltck1*cxzyr+tiltsk1*cxzyi)                  !hr02
  ab2(3)=(28.0_fPrec*ekk)*(tiltsk1*cxzyr-tiltck1*cxzyi)                  !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  qu=(eight*ekk)*(tiltck*cxzyr+tiltsk*cxzyi)                         !hr02
  qv=(eight*ekk)*(tiltck*cxzyi-tiltsk*cxzyr)                         !hr02
  ab1(2)=qu
  ab2(2)=-one*qv                                                   !hr02
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  dyy1=ekk*(tiltc(k)*cxzyr+tilts(k)*cxzyi)
  dyy2=ekk*(tilts(k)*cxzyr-tiltc(k)*cxzyi)                         !hr02
  tiltckuk=tiltck5*tiltc(k)-tiltsk5*tilts(k)
  tiltsk=tiltck5*tilts(k)+tiltsk5*tiltc(k)
  tiltck=tiltckuk
  ab1(8)=(eight*ekk)*(tiltck*xl+tiltsk*zl)                           !hr02
  ab2(8)=(eight*ekk)*(tiltsk*xl-tiltck*zl)                           !hr02
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  ab1(9)=ekk*tiltck
  ab2(9)=ekk*tiltsk
#endif
! end include/kicka09h.f90
