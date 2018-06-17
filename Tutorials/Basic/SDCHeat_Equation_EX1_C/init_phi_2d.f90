subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,r2

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        !        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
        x = prob_lo(1) + (dble(i)) * dx(1)

        r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0

        phi(i,j) = 1.d0 + exp(-r2)
        phi(i,j) = sin(x*3.14159265358979323846d0*2d0)

     end do
  end do

end subroutine init_phi

subroutine err_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="err_phi")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,r2,sym,Tfin,tupi

  Tfin=0.01d0
  tupi=3.14159265358979323846d0*2d0
  sym=(-2.0d0+2.0d0*cos(tupi*dx(1)))/(dx(1)*dx(1))
!  print *,Tfin,tupi,sym
  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        !        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
        x = prob_lo(1) + (dble(i)) * dx(1)

        r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0

        phi(i,j) = phi(i,j)-sin(x*tupi)*exp(Tfin*sym)
     end do
  end do

end subroutine err_phi
