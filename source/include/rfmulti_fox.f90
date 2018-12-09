
  !---- Zero the arrays

fact=1   
do j=1,nordm-1   
  fact=fact*j   
end do

  normal = zero
  skew = zero
  pnl = zero
  psl = zero

  do j=1,napx 
  
  crabamp=ed(ix)*nzz(j)

  if(isSkew .eq. 0) then 
    if(nordm .eq. 1) then
      pnl = pi/2 - crabpase_t
      normal = crabamp/e0f
    else
      pnl =  -crabpase_t
      normal = -crabamp*fact
    endif
  else
    if(nordm .eq. 1) then
      psl = -pi/2 - crabpase_t
      skew = crabamp/e0f
    else
      psl = - crabpase_t
      skew = crabamp*fact
    endif
  endif

  krf = (((one/(clight*(e0f/e0)))*crabfreq)*two)*pi
  

!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT) ;
    
!!FOX  DXI=X(1)*C1M3;
!!FOX  DYI=X(2)*C1M3;

!!FOX  FIELDC(1) = (NORMAL * COS(PNL  - KCRABDA))
!!FOX  FIELDS(1) = (NORMAL * sin(PNL  - KCRABDA))
!!FOX  FIELDC(2) = (SKEW   * cos(PSL  - KCRABDA))
!!FOX  FIELDS(2) = (SKEW   * sin(PSL  - KCRABDA))
   


CP_RE = zero;
CP_IM = zero;

SP_RE = zero;
SP_IM = zero;

    do iord = nordm, 1, -1
      if(iord .eq. nordm) then
!! FOX CP_RE = FIELDC(1);
!! FOX CP_IM = FIELDC(2);
!! FOX SP_RE = FIELDS(1);
!! FOX SP_IM = FIELDS(2);

        Sp1 = FIELDS(1)+imag*FIELDS(2);
      else 
        Cp0 = Cp0 * (x_t+imag*y_t) / (iord)  
        Sp1 = Sp1 * (x_t+imag*y_t) / (iord+1) 

      endif
    enddo
    Sp1 = Sp1 * (x_t+imag*y_t);
    

    dpx = -REAL(Cp0)*c1e3*moidpsv(j);
    dpy = AIMAG(Cp0)*c1e3*moidpsv(j);
    dpt = - krf * REAL(Sp1)*c1e3*e0f;
    
    yv1(j) = yv1(j) + dpx
    yv2(j) = yv2(j) + dpy
    ejv(j) = ejv(j) + dpt
    

    ejf0v(j)=ejfv(j)
    ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
    dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
    yv1(j)=(ejf0v(j)/ejfv(j))*yv1(j)
    yv2(j)=(ejf0v(j)/ejfv(j))*yv2(j)
    if(ithick.eq.1) call envarsv(dpsv,moidpsv,rvv,ekv)

    print *, "aaaaa", dpx, dpy, dpt, Sp1, imag
  end do



