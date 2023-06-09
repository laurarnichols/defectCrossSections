module capcs
  
  use constants, only: dp, HartreeToJ

  implicit none 
  real(kind=dp),parameter :: Kb =  1.38064852d-23
  real(kind=dp),parameter :: pi= 3.14159265358979
  real(kind=dp),parameter :: tpi = 6.2831853071795864769 
  real(kind=dp),parameter :: Thz = 1.0d12
  real(kind=dp),parameter :: hbar = 1.0545718d-34
  real(kind=dp),parameter :: eV = 1.6021766208d-19
  real(kind=dp),parameter :: meV =  1.6021766208d-22
  real(kind=dp),parameter :: q_comvert = 3.920205055d50
  real(kind=dp),allocatable :: ipfreq(:),Sj(:),Mj(:),E(:)
  real(kind=dp) :: temperature,beta,ematrixif
  real(kind=dp) :: limit,alpha,dstep,bin,gamma0, Eif, Ei, Ef, ematrix_real, ematrix_img, lambda
  complex(kind=dp) ::  G1t, wif,wif0,wif1
  character(len=256) :: dummy
  character(len=256) :: Sjinput, Mjinput, Einput
  integer :: nstep, nw, nfreq, nn, nE
   namelist /capcsconf/ ematrixif,temperature,Sjinput, Mjinput, nn,&
                       limit,gamma0,alpha, dstep,nw, bin, Einput, ematrix_real, ematrix_img

contains
subroutine init()
implicit none
character(len=256) :: dummy
open(13,file='input.in')
!read input
read(13,capcsconf)
beta = 1.0d0/Kb/temperature
!Eif = Ef - Ei !in Hartree
!Eif = Eif * eV !in J
!Eif = -0.4*eV
bin = bin*meV
ematrixif = ematrixif*HartreeToJ**2 ! Joule^2
!write(*,*) "gamma0 = ", gamma0
!transform ematrixif from Hartree unit to joule^2
end subroutine

  subroutine readEnergy(E)
  
    implicit none

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: E(:)
      !! Energy for delta function

    ! Local variables:
    integer :: iDum
      !! Dummy integer
    integer :: iE
      !! Loop index


    open(12,file=trim(Einput))

    read(12,*)
    read(12,*) nE

    allocate(E(nE))

    do iE = 1, nE

      read(12,*) iDum, iDum, E(iE) ! in Hartree

    end do

    E(:) = E(:)*HartreeToJ

    close(12)

    return

  end subroutine readEnergy

subroutine readphonon()
implicit none
integer :: ifreq, modeIndex
open(12,file=trim(Sjinput))
read(12,*)nfreq
!read mode number
allocate(Sj(1:nfreq))
allocate(ipfreq(1:nfreq))

do ifreq=1,nfreq
   read(12,*)modeIndex,Sj(ifreq),ipfreq(ifreq) !freq read from Sj.out is f(in Thz)*2pi
end do
!read frequency and Sj
ipfreq=ipfreq*Thz!*tpi convert to Hz*2pi
!transform frequency to Hz
close(12)
!write(*,*) "Frequency Read"
end subroutine

subroutine readelectron()
implicit none
integer :: ifreq
open(12,file=trim(Mjinput))
read(12,*)nfreq
!read mode number
allocate(Mj(1:nfreq))

do ifreq=1,nfreq
   read(12,*) dummy,dummy,dummy,Mj(ifreq) !Mj is in Hartree^2
end do
Mj=Mj*HartreeToJ**2*q_comvert!convert to Joule^2/(kg*m^2)

close(12)
!write(*,*) "Frequency Read"
end subroutine

subroutine testelectron()
implicit none
integer :: ifreq
open(12,file="electron_test.txt")
write(12,*) ematrixif/HartreeToJ**2
do ifreq=1,nfreq
   write(12,*) ifreq," ",Mj(ifreq)/HartreeToJ**2 !Mj is in Hartree^2
end do
close(12)
!write(*,*) "Frequency Read"
end subroutine
end module capcs



function G0_t(inputt) result(G0t) 
use capcs
  integer :: ifreq
  real(kind=dp) :: inputt, nj, omega, e_factor
  complex(kind=dp) G0t, tmp1, tmp2
  tmp1 = 0.0_dp
  tmp2 = 0.0_dp
  do ifreq=1, nfreq
   omega=ipfreq(ifreq)
   nj=1/(exp(hbar*omega*beta)-1)
   tmp1=tmp1+(nj+1)*Sj(ifreq)*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
   !tmp2=tmp2+(nj+1)*Sj(ifreq)/e_factor*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*e_factor*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
  end do
  G0t= tmp1! - exp(tmp2)
!time-dependent line-shape-function, expand it to the first order of Sj, you will get familiar form with Huang-Rhys factor
  !!G0t=G0t*exp(-0.25*(alpha*inputt)**2)
  !narrow Gaussian to simulate delta function, what we care about is the area, we make it narrow so neighboring LSF do not interfere
end function

function sinc(x) result (sincx)
use capcs
    real(kind=dp) :: x, sincx
    if ( x < 1e-20 ) then 
        sincx = 1.0_dp
    else
        sincx = sin(x)/x
    end if
end function

function Aj(t1) result (A)
use capcs
    real(kind=dp) :: t1, wt1, omega
    complex(kind=dp) :: t2, t3, wt2, wt3, A, A1, A2, tmp
    integer :: ifreq
    A = 0.0_dp
    t2 = dcmplx(-t1, -hbar*beta)
    t3 = dcmplx(0.0_dp, -hbar*beta)
    do ifreq=1, nfreq
        omega=ipfreq(ifreq)
        wt1 = omega*t1*0.5_dp
        wt2 = omega*t2*0.5_dp
        wt3 = omega*t3*0.5_dp
        A1 = 2.0_dp*hbar/omega*Sj(ifreq)
        tmp = cos(wt1)*sin(wt2)/sin(wt3)
        A1 = A1*4.0_dp*tmp*tmp
        !(*,*) "A1 = ", A1
        A2 = dcmplx(0.0_dp, -hbar/omega)*cos(wt1-wt2)/sin(wt3)
        !write(*,*) "A2 = ", A2
        A = A+Mj(ifreq)*(A1+A2)*0.25_dp
        !write(*,*) "A2 = ", A2
        !write(*,*) "t = ", t1
        !write(*,*) "omega = ", omega
        !write(*,*) "Sj = ", Sj(ifreq)
    end do
end function

function intg(a, b, c, d, t, dt, tmpa2, tmpb2) result (res)
use capcs
    real(kind=dp) :: t, dt!t1, wt1, omega
    complex(kind=dp) :: res, a, b, c, d, tmpa2, tmpb2!,t2, t3, wt2, wt3, A, A1, A2, tmp
    !res = c*exp(b*t)+d/b*(b*t-1)*exp(b*t)
    res = (c+d/b*(b*t-1))*tmpb2
    !write(*,*) c,b,t
    t = t-dt
    !res = res - (c*exp(b*t)+d/b*(b*t-1)*exp(b*t))
    res = res - (c+d/b*(b*t-1))*tmpa2
    res = res/b
    !res = exp(a-zlog(b))
    !res = res*exp(a)/b 
end function

program captureCS
use capcs
use OMP_LIB
implicit none
integer :: i,j,k,numofcores
real(kind=dp) :: rangemax, dtime, dw, inputt, t1, t2, inta, intb, S1, S2, sinc
complex(kind=dp) :: temp,tmpa1,tmpa2,G0_t,Aj,tmpb1,tmpb2, intg, fa, fb, fc, fd, tmpa3, tmpb3
real(kind=dp)::omega_tmp
real(kind=dp),allocatable :: S(:)
call init()
call readphonon()
call readelectron()
call readEnergy(E)

allocate(S(1:nE))
S = 0.0d0
!do i=1, nE
!  S(i) = 0.0d0
!end do
!call testelectron()
!write(*,*) Aj(0.000001_dp)

!!limit = 1e-4
dtime = 1/Thz
write(*,*)"simps"
write(*,*)dstep
dstep = dtime*dstep
!steplength for time zone integration
!!alpha = alpha*meV/hbar
gamma0 = gamma0*meV/hbar
!gamma0 = gamma0*gamma0 ! Lorentz
!transform alpha to be in Hz
!!rangemax= -log(limit*1.0)/gamma0

nstep=200000000!Ceiling(rangemax/dstep)
write(*,*) nstep*dstep
!maximum time and steps
!!bin = bin*meV
!!Eif = Eif*eV
!energy in joule
!!dw=bin/hbar/float(nw)
!frequency step for energy zone integration
!!inta=(Eif-0.5*bin)/hbar
!!intb=(Eif+0.5*bin)/hbar
!two ends of  frequency interval
!!write(*,*)"hello",nstep/10000,"dstep",dstep
!!write(*,*)"inta",inta,"intb", intb, "dw",dw
!inputt=nstep*dstep*nn
S1=0.0d0!+Real(ematrixif*exp(G0_t(inputt)+cmplx(0.0,inputt*EE/hbar,dp)-gamma0*inputt*inputt))*dstep
!$OMP PARALLEL DO default(shared) private(inputt,t1,t2,S2,i,j,k,omega_tmp,tmpa1,tmpa2,tmpa3,tmpb1,tmpb2,tmpb3)&
& firstprivate(Eif,inta,dstep,nw,dw,nstep) reduction(+:S1,S) 
DO i=1, nstep-1, 2!nstep-1, 0, -1
  inputt=(float(i))*dstep
  !write(*,*) inputt
  t1=inputt + nstep*dstep*nn
   !tmp1=(ematrixif+Aj(inputt))*G0_t(inputt)*exp(cmplx(0.0,inputt*EE/hbar,dp))*exp(-gamma0*inputt*inputt)*sinc(bin/hbar*inputt*0.5)
  !tmp1 =
  !Aj(inputt)*G0_t(inputt)*exp(cmplx(0.0,inputt*EE/hbar,dp))*exp(-gamma0*inputt*inputt)*sinc(bin/hbar*inputt*0.5)
  tmpa1 = (ematrixif)!*sinc(bin/hbar*inputt*0.5)
  !tmpa2 = G0_t(inputt)+cmplx(0.0,inputt*Eif/hbar,dp)-gamma0*inputt*inputt
  tmpa2 = G0_t(t1)-gamma0*t1 ! Lorentz
  !write(*,*) "tmpa1: ",tmpa1
  !write(*,*) "tmpa2: ",tmpa2
  t2 = t1 + dstep
  tmpb1 = (ematrixif)!*sinc(bin/hbar*inputt*0.5)
  !tmpb2 = G0_t(inputt)+cmplx(0.0,inputt*Eif/hbar,dp)-gamma0*inputt*inputt
  tmpb2 = G0_t(t2)-gamma0*t2 ! Lorentz
  !write(*,*) i, " ", t1, " ", t2
  DO j=1, nE
    Eif = E(j)
    tmpa3 = tmpa2+cmplx(0.0,t1*Eif/hbar,dp)
    tmpb3 = tmpb2+cmplx(0.0,t2*Eif/hbar,dp)
    S(j) = S(j) + Real(4d0*tmpa1*exp(tmpa3)+2d0*tmpb1*exp(tmpb3))
    !write(*,*) i, " ", j, " ", Real(4d0*tmpa1*exp(tmpa2+cmplx(0.0,t1*Eif/hbar,dp))+2d0*tmpb1*exp(tmpb2+cmplx(0.0,t2*Eif/hbar,dp)))*dstep 
    !write(*,'(2i4,4ES25.14E2)') i, j, tmpa3, exp(tmpa3)
    ! write(*,'(i4,8ES25.14E2)') i, tmpa1, tmpa2,tmpb1, tmpb2
  END DO
 !time zone integration
END DO
!write(*,*) "nw",nw,"dstep",dstep 
!write(*,*)"energy bin",bin/meV
!write(*,*)(S1*2.0d0/tpi/hbar)*tpi/hbar*ematrixif
!write(*,*) nE, " ", EE, "     ", S1 
!write(*,*) "ematrixif", ematrixif
!write(*,*) S1!*2.0d0/hbar/hbar*ematrixif
!wif0, final result we expect
DO i=1, nE
   Eif = E(i)
   t1 = nstep*dstep*nn
   t2 = (float(nstep))*dstep+nstep*dstep*nn
   S(i) = S(i) + Real(ematrixif*exp(G0_t(t1)+cmplx(0.0,t1*Eif/hbar,dp)-gamma0*t1)) ! Lorentz
   S(i) = S(i) - Real(ematrixif*exp(G0_t(t2)+cmplx(0.0,t2*Eif/hbar,dp)-gamma0*t2)) ! Lorentz
   S(i) = S(i)*2.0d0/3.0d0*dstep
   S(i) = S(i)/hbar/hbar
   write(*,'(f10.5,f7.1,i5,ES35.14E3)') E(i)/eV, temperature, nn, S(i)
END DO
end program
