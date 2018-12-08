
  !---- Zero the arrays
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero


fact=1   
do j=1,nordm-1   
  fact=fact*j   
end do


  do j=1,napx 
  
  crabamp=ed(ix)*nzz(j)
  
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero
  

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
  


  print * , "krf",  "crabfreq", "crabamp",  "phase", pnl, normal
  print * , krf, crabfreq, crabamp, crabph(ix), pnl, normal, "ddddd"
  print *, psl, skew, "skeeewww11"


    x_t = xv1(j)*c1m3 
    y_t = xv2(j)*c1m3

      field_cos(1) = (normal * cos(pnl  - krf * sigmv(j)))
      field_sin(1) = (normal * sin(pnl  - krf * sigmv(j)))
      field_cos(2) = (skew   * cos(psl  - krf * sigmv(j)))
      field_sin(2) = (skew   * sin(psl  - krf * sigmv(j)))
   

    Cm2 = zero; Sm2a = zero; Cm1 = zero; Sm1a = zero;
    Cp0 = zero; Sp0 = zero; Cp1 = zero; Sp1 = zero;

    do iord = nordm-1, 0, -1
      if(iord .eq. nordm-1) then
        Cp0 = field_cos(1)+imag*field_cos(2);
        Sp1 = field_sin(1)+imag*field_sin(2);
      else 
        Cp0 = Cp0 * (x_t+imag*y_t) / (iord+1)  
        Sp1 = Sp1 * (x_t+imag*y_t) / (iord+2) 
        print *, "ooord",  (iord+1)
        print *, "ooorf2",  (iord+2)
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



