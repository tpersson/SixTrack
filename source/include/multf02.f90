!Here it is in angles because the transformation is afterwards

!FOX  YV1J=BBI(I,1)+BBI(I,2)*XL+AAI(I,2)*ZL+ (DKI(I,1)/DKI(I,3))*(XL*XL-0.5*ZL*ZL)*C1M3 ;
!FOX  YV2J=AAI(I,1)-BBI(I,2)*ZL+AAI(I,2)*XL-
!FOX  (DKI(I,1)/DKI(I,3))*(XL*ZL*C1M3+C1M6*(DKI(I,1)/DKI(I,3))*(ZL*ZL*ZL/6.0));

!yv1j=bbiv(1,1,i)+bbiv(2,1,i)*(xlvj+((strack(i)*(10e-3*xlvj**2-0.5*zlvj**2)*10e-3)))+aaiv(2,1,i)*zlvj     !hr03
!yv2j=(aaiv(1,1,i)-bbiv(2,1,i)*zlvj)+aaiv(2,1,i)*xlvj-((strack(i)**2*bbiv(2,1,i)*zlvj**3)/6.0)*10e-6 &
!& -strack(i)*xlvj*zlvj*10e-3*bbiv(2,1,i) 

!FOX  CRKVE=XL ;
!FOX  CIKVE=ZL ;
