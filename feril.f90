module feril
!****************************************************************************
! feril = fortran electron repulsion integral (ERI) library
!
! The main function to be used is coulomb_rep4, a 4-center ERI calculator
! for Gaussian type orbitals (GTOs). Also includes a 2-center ERI, 
! coulomb_rep2, which is a wrapper for coulomb_rep4. 
! This work is based off of the pyquante2 code.
!   https://github.com/rpmuller/pyquante2/blob/master/pyquante2/ints/two.py
! which is a python implementation of:
!    'Gaussian Expansion Methods for Molecular Orbitals' H. Taketa,
!    S. Huzinaga & K. O-ohata. H. Pys. Soc. Japan, 21, 2313, 1966.
! Special functions are based off of 
!  -Gamma & Incomplete Gamma Functions:
!    'Computation of Special Functions' by S. Zhang & J. Jin
!  -Boys Function:
!    'Molecular Electronic-Structure Theory' by T. Helgaker, P. Jorgensen
!                                               & J. Olsen
!
! Author:          Date:           Version:
! Aaron Mahler     04Nov2018            0.0
!****************************************************************************

! @TODO: add timing to see what takes the longest

  use kinds, only : dp
  use constants, only: pi
  implicit none

  ! Hard coded limit so b_arrays can be on the stack.
  ! INT(KIND=8) limit is @ vmax_am_tot=29, because fact_ratio2 is maxed 
  ! @ (28,14). It also disagrees with pyquante2 when vmax_am_tot = 25 by 
  ! ~0.2% for maxed out input angular momentum.
  integer, parameter :: vmax_am_tot = 21 ! 4*5 + 1 (4 h-type orbitals)
  real(dp), dimension(vmax_am_tot)   :: bax, bay, baz
  real(dp), dimension(0:vmax_am_tot) :: fbar ! array of fboys evaluations
  real(dp) :: ericoef = 2.0_dp * (pi**2.5_dp)
  real(dp), parameter :: taytol = 0.01_dp  ! tolerance to eval fboys with Taylor series
  real(dp), parameter :: anatol = 150.0_dp ! tol to analytical eval fboys for big arg

  contains
function fboys(nin, xin) result(ret)
!****************************************************************************
! fboys evaluates the Boys function(n,x) := \int_0^1 exp(-xt^2)t^(2n) dt
!
!  for small x<0.01, Taylor series at x=0
!  for large x>150, analytic form of \int_0^\infty exp(-xt^2)t^(2n) dt
!  for intermediate range use F_n(x) = \gamma(n+1/2,x) / ( 2x^(n+1/2) ),
!    where \gamma is the lower incomplete gamma function  
!
!  Parameters:
!    Input, real(dp) :: nin, xin ! both inputs should be >= 0
!    Output, real(dp) :: ret ! Boys func, range is 1/(2n+1) >= F_n(x) > 0
!****************************************************************************
  use kinds, only : dp
  use constants, only : pi
  implicit none
  ! Params
  real(dp), intent(in) :: nin, xin
  ! Locals
  real(dp) :: r1
  integer  :: ii
  ! Return
  real(dp) :: ret

  ! input error checking, change nin to int type?
  if ( (xin .lt. 0.0_dp) .or. (nin .lt. 0.0_dp) ) then
    write(6,'("ERR in fboys, negative input not allowed")')
    ret = -1.0_dp
    return
  end if

  if ( xin .lt. taytol ) then
  ! Use truncated series expansion at x=0 for small x
  ! F_n(x) = \sum_{k=0}^{\infty} (-x)^k / ( k! (2n+2k+1) )

    ret = 1.0_dp / ( 2.0_dp*nin + 1.0_dp )
    !r1 = 1.0_dp
    !do ii = 1, 5 ! loop method for testing different truncation lengths
    !  r1 = real(ii, dp)*r1
    !  ret = ret + (-xin)**ii / (r1 * (2.0_dp*nin + 2.0_dp*ii + 1.0_dp))
    !end do
    ! manually unroll loop
    ! truncation at k=5 tested against Mathematica up to n=50 for x=0.01
    ret = ret - xin / (2.0_dp*nin + 3.0_dp)
    ret = ret + xin**2 / (4.0_dp*nin + 10.0_dp)
    ret = ret - xin**3 / (12.0_dp*nin + 42.0_dp)
    ret = ret + xin**4 / (48.0_dp*nin + 216.0_dp)
    ret = ret - xin**5 / (240.0_dp*nin + 1320.0_dp)
    !write(6,*) 'fboys taylor'
    return

  else if ( xin .gt. anatol ) then
  ! F_n(x) \approx (2n-1)!! / 2^(n+1) * Sqrt( \pi / x^(2n+1) ); for large x
  ! tested against Mathematica up to n=50 for x=150
    r1 = 1.0_dp
    if ( nin .ge. 2.0_dp ) then ! double factorial
        do ii = 2*int(nin)-1, 1, -2
            r1 = r1 * real(ii, dp)
        end do
    end if
    ret = r1 * 2.0_dp**(-nin-1) * sqrt( pi * xin**(-2*nin-1) )
   ! write(6,*) 'fboys large'
    return

  end if

  ! else use incomplete gamma func
  r1 = nin + 0.5_dp
  ret = incgaml(r1, xin) / ( 2.0_dp * xin**r1 )
 ! write(6,*) 'fboys incgaml'

  return
end function fboys


function gam_zj(xin) result(ga)
!****************************************************************************
! gam_zj evaluates the Gamma Function(n) := 
!    \int_0^{\infty} exp(-t) t^{n-1) dt
!
!  Reference:
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
!
!  Parameters:
!    Input, real (dp) Xin, the argument.
!    Xin must not be 0, or any negative integer.
!    Output, real (dp) GA, the value of the Gamma function.
!****************************************************************************
  use kinds, only : dp
  use constants, only : pi

  implicit none
  ! Params
  real(dp), intent(in) :: xin
  ! Locals
  real(dp), dimension (26) :: gra = (/ &
    1.0E+00_dp, &
    0.5772156649015329E+00_dp, &
   -0.6558780715202538E+00_dp, &
   -0.420026350340952E-01_dp, &
    0.1665386113822915E+00_dp, &
   -0.421977345555443E-01_dp, &
   -0.96219715278770E-02_dp, &
    0.72189432466630E-02_dp, &
   -0.11651675918591E-02_dp, &
   -0.2152416741149E-03_dp, &
    0.1280502823882E-03_dp, & 
   -0.201348547807E-04_dp, &
   -0.12504934821E-05_dp, &
    0.11330272320E-05_dp, &
   -0.2056338417E-06_dp, & 
    0.61160950E-08_dp, &
    0.50020075E-08_dp, &
   -0.11812746E-08_dp, &
    0.1043427E-09_dp, & 
    0.77823E-11_dp, &
   -0.36968E-11_dp, &
    0.51E-12_dp, &
   -0.206E-13_dp, &
   -0.54E-14_dp, &
    0.14E-14_dp, &
    0.1E-15_dp /)

  real(dp) :: gr, rr, rz
  integer  :: ik, im, im1
  ! Return
  real(dp) :: ga

  ifgam1: if ( xin == aint(xin) ) then

    if ( 0.0_dp < xin ) then
      ga = 1.0_dp
      im1 = int(xin) - 1
      do ik = 2, im1
        ga = ga * real(ik, dp)
      end do
    else
      ga = 1.0E+300_dp
    end if

  else ifgam1

    if ( 1.0_dp < abs(xin) ) then
      rz = abs(xin)
      im = int(rz)
      rr = 1.0_dp
      do ik = 1, im
        rr = rr * ( rz - real(ik, kind=dp) )
      end do
      rz = rz - real(im, kind=dp)
    else
      rz = xin
    end if

    gr = gra(26)
    do ik = 25, 1, -1
      gr = gr * rz + gra(ik)
    end do

    ga = 1.0_dp / ( gr * rz )

    if ( 1.0_dp < abs(xin) ) then
      ga = ga * rr
      if ( xin < 0.0_dp ) then
        ga = - pi / ( xin * ga * sin(pi*xin) )
      end if
    end if

  end if ifgam1

  return
end function gam_zj

function incgaml(ain, xin) result(gout)
!****************************************************************************
! incgaml evaluates the Lower Incomplete Gamma Function(a,x) := 
!   \int_0^x exp(-t) t^{n-1) dt
! 
! Reference:
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996, ISBN: 0-471-11963-6, LC: QA351.C45.
! https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
!
!  Parameters:
!    Input, real ( kind = 8 ) ain, the parameter.
!    Input, real ( kind = 8 ) xin, the argument.
!    Output, real(dp) gout
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  real(dp), intent(in) :: ain, xin
  ! Locals
  real(dp) :: rr, rs
  real(dp) :: t0, xam
  integer :: ik
  ! Return
  real(dp) :: gout

  xam = -xin + ain*log(xin)

  if ( 700.0_dp < xam .or. 170.0_dp < ain ) then
    write(6, '(/,a)') 'ERR (incgaml), arg(s) too large! returning 0.0'
    gout = 0.0_dp
    return
  end if

  if ( xin .lt. 1.0e-50_dp ) then
  ! Tested against Mathematica, Gamma[1,0,10^-50] = 10e-50
  ! which is small enough to cutoff for this module
    gout = 0.0_dp

  else if ( xin .le. (1.0_dp+ain) ) then
    rs = 1.0_dp / ain
    rr = rs
    dok1: do ik = 1, 60
      rr = rr * xin / ( ain + real(ik, dp) )
      rs = rs + rr
      if ( abs(rr/rs) .lt. 1.0E-15_dp ) then
          exit dok1
      end if
    end do dok1
    gout = exp(xam) * rs

  else if ( (1.0_dp+ain) .lt. xin ) then
    t0 = 0.0_dp
    do ik = 60, 1, -1
      t0 = ( real(ik,dp) - ain ) / ( 1.0_dp + real(ik, dp)/(xin+t0) )
    end do 
    gout = gam_zj(ain) - ( exp(xam) / (xin+t0) )
    !gout = gamma(ain) - ( exp(xam) / (xin+t0) )

  else ! usually a NaN
    write(*,*) 'ERR (incgaml) undefined for requested input' 
    write(*,*) 'xin: ', xin
    write(*,*) 'ain: ', ain
    stop 4 ! SIGILL: illegal instruction
  end if

  return
end function incgaml


function fact_ratio2(ain, bin) result(ret)
!****************************************************************************
! Factorial ratio 2 := a! / ( b! (a-2b)! )
!****************************************************************************
  implicit none
  ! Params
  integer, intent(in) :: ain, bin ! ain >= 2bin
  ! Locals
  integer :: i1, imx, imn
  ! Return
  integer(kind=8) :: ret

  if ( (ain.eq.0) .or. (ain.eq.1) .or. &
       (bin.eq.0) ) then
      ret =  1
      return
  else if (ain.gt.28) then
      write(6,'("ERR, fact_ratio2 cannot handle numbers > 28")')
      ret = 0
      return
  end if

  imx = max(bin, ain - 2*bin)
  imn = min(bin, ain - 2*bin)

  ret = 1
  do i1 = imx+1, ain
      ret = ret * int(i1, 8)
  end do
  if ( imn .gt. 1) then
      do i1 = 2, imn
          ret = ret / int(i1, 8)
      end do
  end if

  return
end function fact_ratio2


function binomcof(ain, bin) result(ret)
!****************************************************************************
! binomial coefficient := a! / ( b! (a-b)! )
!****************************************************************************
  implicit none
  ! Params
  integer, intent(in) :: ain, bin
  ! Locals
  integer :: i1, imx, imn
  ! Return
  integer :: ret

  if ( (ain.eq.0) .or. (ain.eq.1) .or. &
       (bin.eq.0) ) then
      ret =  1
      return
  end if

  imx = max(bin, ain - bin)
  imn = min(bin, ain - bin)

  ret = 1
  do i1 = imx+1, ain
      ret = ret * i1
  end do
  if ( imn .gt. 1) then
      do i1 = 2, imn
          ret = ret / i1
      end do
  end if

  return
end function binomcof


function binom_prefac_v2(si, iai) result(ret)
!****************************************************************************
! integral prefactor with the binom coefficient from Augspruger & Dykstra
! adapted from 
! https://github.com/rpmuller/pyquante2/blob/master/pyquante2/ints/one.py
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: si, iai
  ! Return
  real(dp) :: ret

  ret = real(binomcof(iai, si), dp)
  return
end function binom_prefac_v2


function binom_prefac(si, iai, ibi, xpai, xpbi) result(ret)
!****************************************************************************
! integral prefactor with the binom coefficient from Augspruger & Dykstra
! adapted from 
! https://github.com/rpmuller/pyquante2/blob/master/pyquante2/ints/one.py
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: si, iai, ibi
  real(dp), intent(in) :: xpai, xpbi
  ! Locals
  integer :: i1
  ! Return
  real(dp) :: ret

  ret = 0.0_dp
  do i1 = 0, si
    if ( ((si-iai).le.i1) .and. (i1.le.ibi) ) then
      ret = ret + real(binomcof(iai, si-i1), dp) * &
                  real(binomcof(ibi, i1), dp) * &
                  xpai ** (iai - si + i1) * &
                  xpbi ** (ibi - i1)
    end if
  end do

  return
end function binom_prefac

function fb_v2(l1in, rin, gin) result(ret)
!****************************************************************************
! @TODO: combine with b0, rename to something more relevant
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: l1in, rin
  real(dp), intent(in) :: gin
  ! Locals
  real(dp) :: r1
  ! Return
  real(dp) :: ret

  ret = 0.0_dp
      !write(6,*) 'yee', i1, iin, l1in
  !if (iin.le.l1in) then !.and. (i1.le.0) ) then
  !if (iin.eq.l1in) then !.and. (i1.le.0) ) then
    !r1 = real(binomcof(l1in, iin), dp) * &
    !       0.0_dp ** (l1in - iin) !* &
    ret = real(fact_ratio2(l1in, rin), dp) * (4.0_dp*gin)**(rin-l1in)
  !end if

  return 
end function fb_v2


function b0(iin, rin, gin) result(ret)
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: iin, rin
  real(dp), intent(in) :: gin
  ! Return
  real(dp) :: ret

  ret = real(fact_ratio2(iin, rin), dp) * (4.0_dp*gin)**(rin-iin)

  return
end function b0


function fb(iin, l1in, l2in, pin, ain, bin, rin, gin) result(ret)
!****************************************************************************
! @TODO: combine with b0, rename to something more relevant
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: iin, l1in, l2in, rin
  real(dp), intent(in) :: pin, ain, bin, gin
  ! Return
  real(dp) :: ret

  ret = b0(iin, rin, gin) * &
        binom_prefac(iin, l1in, l2in, pin-ain, pin-bin)

  return 
end function fb


function b_term(i1i, i2i, r1i, r2i, ui, l1i, l2i, l3i, l4i, pxi, axi, &
                bxi, qxi, cxi, dxi, gam1, gam2, deli) result(ret)
!****************************************************************************
! Logic for filling elements of b_array
!****************************************************************************
   use kinds, only : dp
   implicit none
   ! Params
   integer,  intent(in) :: i1i, i2i, r1i, r2i, ui
   integer,  intent(in) :: l1i, l2i, l3i, l4i
   real(dp), intent(in) :: pxi, axi, bxi
   real(dp), intent(in) :: qxi, cxi, dxi
   real(dp), intent(in) :: gam1, gam2, deli
   ! Locals
   real(dp) :: sgn1, p1, p2, fb1, fb2, frat
   ! Return
   real(dp) :: ret

   sgn1 = real( (-1)**(i2i+ui), dp)
   p1 = (qxi-pxi) ** (i1i + i2i - 2*(r1i+r2i+ui))
   p2 = deli ** (i1i + i2i - ui - 2*(r1i+r2i))
   fb1 = fb(i1i, l1i, l2i, pxi, axi, bxi, r1i, gam1)
   fb2 = fb(i2i, l3i, l4i, qxi, cxi, dxi, r2i, gam2)
   frat = real( fact_ratio2( i1i + i2i - 2*(r1i+r2i), ui), dp)

   ret = fb1 * fb2 * frat * sgn1 * p1 / p2

   return
end function b_term


subroutine b_array(l1i, l2i, l3i, l4i, pin, ain, bin, qin, cin, din, &
                   g1i, g2i, deli, bout)
!****************************************************************************
! @TODO: combine with b_term
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: l1i, l2i, l3i, l4i
  real(dp), intent(in) :: pin, ain, bin
  real(dp), intent(in) :: qin, cin, din
  real(dp), intent(in) :: g1i, g2i, deli
!  real(dp), allocatable, intent(out) :: bout(:)
  real(dp), intent(out) :: bout(vmax_am_tot)
  ! Locals
  integer :: i1, i2, ir1, ir2, iu, ij
  integer :: imax

!  if ( allocated(bout) ) deallocate(bout)
  imax = l1i + l2i + l3i + l4i + 1
!  allocate( bout(imax) )

  bout(:) = 0.0_dp
  do i1  = 1, l1i + l2i + 1
  do i2  = 1, l3i + l4i + 1
  do ir1 = 1, (i1-1)/2 + 1
  do ir2 = 1, (i2-1)/2 + 1
  do iu  = 1, (i1+i2-2)/2 - ir1 - ir2 + 3
    ij = i1 + i2 - 2*(ir1+ir2-2) - iu ! +4 for 1-based indexing
    bout(ij) = bout(ij) + b_term(i1-1, i2-1, ir1-1, ir2-1, iu-1, &
                                 l1i, l2i, l3i, l4i, &
                                 pin, ain, bin, qin, cin, din, &
                                 g1i, g2i, deli)
!    write(6,'(2x,6(i2,x,","))') i1-1, i2-1, ir1-1, ir2-1, iu-1, ij-1
  end do
  end do
  end do
  end do
  end do

end subroutine b_array


subroutine gauss_prod_center(alpa, acen, alpb, bcen, ocen)
!****************************************************************************
! fill ocen with the center of the combined gaussians
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  real(dp), intent(in)  :: alpa, alpb       ! exponents
  real(dp), intent(in)  :: acen(3), bcen(3) ! centers
  real(dp), intent(out) :: ocen(3)          ! output center
  ! Locals
  real(dp) :: alpab

  alpab = 1.0_dp / (alpa + alpb)
  ocen(1) = (alpa*acen(1) + alpb*bcen(1)) * alpab
  ocen(2) = (alpa*acen(2) + alpb*bcen(2)) * alpab
  ocen(3) = (alpa*acen(3) + alpb*bcen(3)) * alpab

end subroutine gauss_prod_center

  
function diff3_norm2(vec1, vec2) result(ret)
!****************************************************************************
! find norm^2 of the difference of two 3d vectors
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  real(dp), intent(in) :: vec1(3), vec2(3)
  ! Return
  real(dp) :: ret

  ret = ( vec1(1) - vec2(1) )**2 + &
        ( vec1(2) - vec2(2) )**2 + &
        ( vec1(3) - vec2(3) )**2

  return
end function diff3_norm2


subroutine b_array_v2(l1i, l3i, cmain, g1i, g2i, deli, bout)
!****************************************************************************
! b_array just for 2 center, all functions inlined for speed
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  integer,  intent(in) :: l1i, l3i
  real(dp), intent(in) :: cmain ! cin - ain
  real(dp), intent(in) :: g1i, g2i, deli
  real(dp), intent(out) :: bout(vmax_am_tot)
  ! Locals
  integer :: ir1, ir2, iu, ij
  integer :: l1p3
  real(dp) :: sgn1, p1, p2, fb1, fb2, frat
  real(dp) :: fg1, fg2

  bout(:) = 0.0_dp

  l1p3 = l1i + l3i
  fg1 = 4.0_dp * g1i
  fg2 = 4.0_dp * g2i

  do ir1 = 0, l1i/2
  do ir2 = 0, l3i/2
  do iu  = 0, (l1i+l3i)/2 - ir1 - ir2
    ij = l1i + l3i - 2*(ir1+ir2) - iu + 1 

    sgn1 = real( (-1)**(l3i+iu), dp )
    p1 = cmain ** (l1p3 - 2*(ir1+ir2+iu))
    p2 = deli ** (l1p3 - iu - 2*(ir1+ir2))
    fb1 = real(fact_ratio2(l1i, ir1), dp) * fg1**(ir1-l1i)
    fb2 = real(fact_ratio2(l3i, ir2), dp) * fg2**(ir2-l3i)
    frat = real( fact_ratio2( l1p3 - 2*(ir1+ir2), iu), dp)

    bout(ij) = bout(ij) + ( fb1 * fb2 * frat * sgn1 * p1 / p2 )
  end do
  end do
  end do

  return
end subroutine b_array_v2


function coulomb_rep2_v2( acen, amom, alpa, &
                          ccen, cmom, alpc ) result(ret)
!****************************************************************************
! 2-center ERI, rewritten from coulomb_rep4, not just a wrapper
! take advantage of 2-center to speed up b_array
! switch to downward recursion of fboys so only eval once
!
! Parameters:
!  Input, real(dp), dim(3) :: acen, ccen ! origins of the GTOs
!  Input, real(dp) :: alpa, alpc         ! alpha(abcd), exponents
!  Input, integer, dim(3) :: amom, cmom  ! ang. mom. of GTOs
!  Output, real(dp) :: ret               ! the value of the ERI
!****************************************************************************
  use kinds, only : dp
  use constants, only : pi
  implicit none
  ! Params
  real(dp), intent(in) :: acen(3), ccen(3) ! origins
  real(dp), intent(in) :: alpa, alpc       ! exponents
  integer,  intent(in) :: amom(3), cmom(3) ! ang. mom.
  ! Locals
  real(dp), parameter :: rtol = 1.0e-8_dp
  integer :: ii, jj, kk
  real(dp) :: rpq2
  real(dp) :: delt, fbarg, fbnum
!  real(dp), dimension(vmax_am_tot) :: bax, bay, baz
!  real(dp), dimension(0:vmax_am_tot) :: fbar
  integer  :: fbarsz
  ! Return
  real(dp) :: ret

  call start_clock('crep2v2')

  ret = 0.0_dp
  rpq2 = diff3_norm2(acen, ccen)                 ! find p & q diff norm
  delt = 0.25_dp * ( 1.0_dp/alpa + 1.0_dp/alpc ) ! b_array arg
  fbarg = 0.25_dp * rpq2 / delt                  ! fboys arg

  if ( (alpa.lt.rtol) .or. (alpc.lt.rtol) .or. (delt.lt.rtol) ) then
      write(6,'("ERR (coulomb_rep2), invalid input, returning 0")')
      ret = 0.0_dp
      return
  end if

  !call start_clock('barrays')
  call b_array_v2(amom(1), cmom(1), ccen(1)-acen(1), alpa, alpc, delt, bax)
  call b_array_v2(amom(2), cmom(2), ccen(2)-acen(2), alpa, alpc, delt, bay)
  call b_array_v2(amom(3), cmom(3), ccen(3)-acen(3), alpa, alpc, delt, baz)
  !call stop_clock('barrays')

  !call start_clock('fboys comp')
  fbarsz = sum( amom + cmom ) 
  if ( (fbarg.lt.taytol) .or. (fbarg.gt.anatol) ) then
      ! taylor series and analytical are fast, go ahead and calculate
      do ii = 0, fbarsz
          fbar(ii) = fboys( real(ii,dp), fbarg )
      end do
  else ! use downward recursion
      fbar(fbarsz) = fboys( real(fbarsz,dp), fbarg )
      if ( fbarsz .gt. 0 ) then
          do ii = fbarsz-1, 0, -1
              fbar(ii) = 2*fbarg*fbar(ii+1) + exp(-fbarg)
              fbar(ii) = fbar(ii) / real( 2*ii + 1, dp )
          end do
      end if
  end if
  !call stop_clock('fboys comp')

  !call start_clock('fboys contract')
  if ( fbarsz .gt. 0 ) then
      do ii = 1, amom(1) + cmom(1) + 1
      do jj = 1, amom(2) + cmom(2) + 1
      do kk = 1, amom(3) + cmom(3) + 1
        ret = ret + bax(ii)*bay(jj)*baz(kk)*fbar(ii+jj+kk-3)
      end do
      end do
      end do
  else ! s-types, no need to contract
      ! I can't tell from run times if eval'ing s-types separately is any faster
      ret = bax(1)*bay(1)*baz(1)*fbar(0)
  end if
  !call stop_clock('fboys contract')

  ret = ret * ericoef / ( alpa*alpc*sqrt(alpa+alpc) )

  call stop_clock('crep2v2')
  return
end function coulomb_rep2_v2


function coulomb_rep4(acen, amom, alpa, bcen, bmom, alpb, &
                      ccen, cmom, alpc, dcen, dmom, alpd) result(ret)
!****************************************************************************
! 4-center ERI, adapted from pyquante2/ints/two.py
!
! Parameters:
!  Input, real(dp), dim(3) :: acen, bcen, ccen, dcen ! origins of the GTOs
!  Input, real(dp) :: alpa, alpb, alpc, alpd ! alpha(abcd), exponents
!  Input, integer, dim(3) :: amom, bmom, cmom, dmom ! ang. mom. of GTOs
!  Output, real(dp) :: ret ! the value of the ERI
!****************************************************************************
  use kinds, only : dp
  use constants, only : pi
  implicit none
  ! Params
  real(dp), intent(in) :: acen(3), bcen(3), ccen(3), dcen(3) ! origins
  real(dp), intent(in) :: alpa, alpb, alpc, alpd             ! exponents
  integer,  intent(in) :: amom(3), bmom(3), cmom(3), dmom(3) ! ang. mom.
  ! Locals
  real(dp), parameter :: rtol = 1.0e-8_dp
  integer :: ii, jj, kk
  real(dp) :: rab2, rcd2, rpq2
  real(dp) :: gam1, gam2, delt, fbarg, fbnum
  real(dp) :: pcen(3), qcen(3)
!  real(dp), allocatable :: bax(:), bay(:), baz(:)

!  real(dp), dimension(vmax_am_tot) :: bax, bay, baz
!  real(dp), dimension(0:vmax_am_tot) :: fbar
  integer  :: fbarsz
  ! Return
  real(dp) :: ret

  call start_clock('coulomb_rep4')

  ret = 0.0_dp

  ! Combine a & b => p
  rab2 = diff3_norm2(acen, bcen)
  gam1 = alpa + alpb
  call gauss_prod_center(alpa, acen, alpb, bcen, pcen)

  ! Combine c & d => q
  rcd2 = diff3_norm2(ccen, dcen)
  gam2 = alpc + alpd
  call gauss_prod_center(alpc, ccen, alpd, dcen, qcen)

  rpq2 = diff3_norm2(pcen, qcen)                 ! find p & q diff norm
  delt = 0.25_dp * ( 1.0_dp/gam1 + 1.0_dp/gam2 ) ! b_array arg
  ! @TODO?
!  delt = 0.25_dp * ( gam1**-1 + gam2**-1 ) ! b_array arg
  fbarg = 0.25_dp * rpq2 / delt                  ! fboys arg

  if ( (gam1.lt.rtol) .or. (gam2.lt.rtol) .or. (delt.lt.rtol) ) then
      write(6,'("ERR (coulomb_rep4), invalid input, returning 0")')
      ret = 0.0_dp
      return
  end if

  call start_clock('barrays')
  call b_array(amom(1), bmom(1), cmom(1), dmom(1), &
               pcen(1), acen(1), bcen(1), qcen(1), ccen(1), dcen(1), &
               gam1, gam2, delt, bax)
  call b_array(amom(2), bmom(2), cmom(2), dmom(2), &
               pcen(2), acen(2), bcen(2), qcen(2), ccen(2), dcen(2), &
               gam1, gam2, delt, bay)
  call b_array(amom(3), bmom(3), cmom(3), dmom(3), &
               pcen(3), acen(3), bcen(3), qcen(3), ccen(3), dcen(3), &
               gam1, gam2, delt, baz)
  call stop_clock('barrays')

  call start_clock('fboys contract')

  fbarsz = sum( amom + bmom + cmom + dmom ) 
  do ii = 0, fbarsz
      fbar(ii) = fboys( real(ii,dp), fbarg )
  end do

  do ii = 1, amom(1) + bmom(1) + cmom(1) + dmom(1) + 1
  do jj = 1, amom(2) + bmom(2) + cmom(2) + dmom(2) + 1
  do kk = 1, amom(3) + bmom(3) + cmom(3) + dmom(3) + 1
    !fbnum = real( ii + jj + kk - 3, dp )
    !ret = ret + bax(ii)*bay(jj)*baz(kk)*fboys(fbnum, fbarg)
    ret = ret + bax(ii)*bay(jj)*baz(kk)*fbar(ii+jj+kk-3)
  end do
  end do
  end do
  call stop_clock('fboys contract')

  ret = ret * 2.0_dp * (pi**2.5_dp) / ( gam1*gam2*sqrt(gam1+gam2) ) * &
        exp( -alpa * alpb * rab2 / gam1 ) * &
        exp( -alpc * alpd * rcd2 / gam2 )

  call stop_clock('coulomb_rep4')
  return
end function coulomb_rep4


function coulomb_rep2(acen, amom, alpa, bcen, bmom, alpb) result(ret)
!****************************************************************************
! 2-center ERI, wrapper for coulomb_rep4
!****************************************************************************
  use kinds, only : dp
  implicit none
  ! Params
  real(dp), intent(in) :: acen(3), bcen(3) ! origins
  real(dp), intent(in) :: alpa, alpb       ! exponents
  integer,  intent(in) :: amom(3), bmom(3) ! ang. mom.
  ! Locals
  integer :: am0(3)
  real(dp) :: alp0
  ! Return
  real(dp) :: ret

  am0 = (/ 0, 0, 0 /)
  alp0 = 0.0_dp

  ret = coulomb_rep4(acen, amom, alpa, acen, am0, alp0, &
                     bcen, bmom, alpb, bcen, am0, alp0)

  return
end function coulomb_rep2


end module feril
