subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, cch1, &
     fx, fxl1, fxh1, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(1),hi(1), coord_type
  integer          :: ccl1, cch1
  integer          :: fxl1, fxh1
  integer          :: fyl1, fyh1
  double precision :: cc(ccl1:cch1)
  double precision :: fx(fxl1:fxh1)
  double precision :: dx(1), problo(1)

  ! Local variables
  integer          :: i
  double precision :: rlo,rhi,rcen

  if (coord_type .eq. 0) then ! Cartesian

     do i=lo(1),hi(1)
        cc(i) = 0.5d0 * ( fx(i) + fx(i+1) )
     enddo
  
  else if (coord_type .eq. 1) then ! Cylindrical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo*fx(i) + rhi*fx(i+1) ) / rcen
     enddo
  
  else  ! Spherical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo**2 * fx(i) + rhi**2 * fx(i+1) ) / rcen**2
     enddo

  end if

end subroutine bl_avg_fc_to_cc

