
! This file defines the Fortran routines required for AUTO-07P to perform bifurcation
! analysis of our immune model.


! Detailed Breakdown of the Routines:

!     FUNC: Defines the system of differential equations.

!     STPNT: Sets the initial conditions for the variables and the parameters.

!     BCND: Defines the boundary conditions. Since this is an ordinary differential 
!           equation (ODE) problem, there are no boundary conditions, so this is left empty.

!     ICND: Defines any integral constraints. For this problem, there are no integral 
!           constraints, so this routine remains empty as well.

!     FOPT: Defines an optional function that can be used for optimization or additional 
!           output. This is not needed for our system, so we leave it as a placeholder.

!     PVLS: This routine can be used to monitor certain variables or parameters during the 
!           continuation process.
!           Here, we output the current values of SS and II during continuation. This can help 
!           track how these variables evolve as we vary the parameters (e.g., ββ).

! Additional Notes:

!     par array: This is the parameter array p(1) = beta1_scaled, p(2) = omega_scaled
!     u array: This represents the state variables:
!              u(1)= M0, u(2)= M1, u(3)= M2, 
!              u(4)= K0, u(5)= KA, 
!              u(6)= E0, u(7)= EF, u(8)= EA

! -------------------------------------------------------------------------------------------------

! FUNC: Defines the system of differential equations.

! Input arguments :
!      ndim   :   Dimension of the algebraic or ODE system 
!      u      :   State variables
!      icp    :   Array indicating the free parameter(s)
!      par    :   Equation parameters

! Values to be returned :
!      f      :   Equation or ODE right hand side values

! Normally unused Jacobian arguments : ijac, dfdu, dfdp (see manual)

subroutine FUNC(ndim, u, icp, par, ijac, f, dfdu, dfdp)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndim, ijac, icp(*)
    DOUBLE PRECISION, INTENT(IN) :: u(ndim), par(*)
    DOUBLE PRECISION, INTENT(OUT) :: f(ndim)
    DOUBLE PRECISION, INTENT(INOUT) :: dfdu(ndim,ndim), dfdp(ndim,*)

    ! Fixed parameters
    DOUBLE PRECISION :: muM_scaled = 0.01 ! muM/MX, muM = 1.0e4, MC = 1.0e6
    DOUBLE PRECISION :: muK_scaled = 3.5e-3 ! muK/KC, muK = 3.5e3, KC = 1.0e6
    DOUBLE PRECISION :: etaE = 1.0
    DOUBLE PRECISION :: etaK = 0.02
    DOUBLE PRECISION :: etaM = 0.2
    DOUBLE PRECISION :: deltaM = 0.02
    DOUBLE PRECISION :: deltaK = 0.02
    DOUBLE PRECISION :: deltaE0 = 1.0
    DOUBLE PRECISION :: deltaE = 0.14
    DOUBLE PRECISION :: sigma_scaled = 1.0e4 ! sigma*EC, sigma = 1e-5, EC = 1.0e9
    DOUBLE PRECISION :: CM1_scaled = 0.1 ! CM1/MC, CM1 = 1.0e5, MC = 1.0e6
    DOUBLE PRECISION :: CM2_scaled = 1.0e-4 ! CM2/MC, CM2 = 100.0, MC = 1.0e6
    DOUBLE PRECISION :: CKA_scaled = 1.0e-1 ! CKA/KC, CKA = 1.0e5, KC = 1.0e6
    DOUBLE PRECISION :: beta2 = 1.0e-3
    DOUBLE PRECISION :: beta12 = 5.0e-5
    DOUBLE PRECISION :: beta21 = 5.0e-5
    DOUBLE PRECISION :: thetaM = 5.0e-8
    DOUBLE PRECISION :: thetaK = 0.7
    DOUBLE PRECISION :: gamma = 0.8
    DOUBLE PRECISION :: rhoF = 0.1
    DOUBLE PRECISION :: muE_scaled = 1.26e-6 ! muE/EC, muE = 1260.0, EC = 1.0e9
    DOUBLE PRECISION :: rho0 = 0.2
  
    ! Extract the state variables and continuation parameters
    DOUBLE PRECISION M0, M1, M2, K0, KA, E0, EF, EA ! State variables
    DOUBLE PRECISION beta1_scaled, omega_scaled ! Parameters (of interest)
    ! p(1) = beta1_scaled = beta1*EC, p(2) = omega_scaled = omega*KC (or omega*MC)
    beta1_scaled = 10**(par(1) +9 )
    omega_scaled = 10**(par(2) +6 )
    
    M0 = u(1)
    M1 = u(2)
    M2 = u(3)
    K0 = u(4)
    KA = u(5)
    E0 = u(6)
    EF = u(7)
    EA = u(8)

    ! System equations ---------------------------------------------------------------------
    f(1) = muM_scaled + etaM*(1.0D0-(M0+M1+M2))*(M0+M1+M2)-beta1_scaled*(EF+EA)*M0-thetaM*(KA/(CKA_scaled+KA))*M0-beta2*M0-deltaM*M0
    f(2) = beta1_scaled*(EF+EA)*M0+thetaM*(KA/(CKA_scaled+KA))*M0+beta21*M2-beta12*M1-deltaM*M1                           
    f(3) = beta2 * M0 + beta12 * M1 - beta21 * M2 - deltaM * M2                                                                  
    f(4) = muK_scaled-thetaK*K0*(M1/(CM1_scaled+M1))-deltaK*K0                                                                   
    f(5) = thetaK*K0*(M1/(CM1_scaled+M1))+etaK*(1.0D0-(KA+K0))*KA-sigma_scaled*(EF+EA)*KA-deltaK*KA                              
    f(6) = muE_scaled-deltaE0*(rho0*E0+(1-rho0)*E0)                                                                       
    f(7) = deltaE0*rho0*E0-omega_scaled*(gamma*KA+(1.0D0-gamma)*M1)*EF-rhoF*EF-deltaE*EF                            
    f(8) = rhoF*EF+etaE*(1.0D0-EA)*(M2/(CM2_scaled+M2))*EA-omega_scaled*(gamma*KA+(1.0D0-gamma)*M1)*EA-deltaE*EA

    IF(ijac.EQ.0)RETURN
  
    ! The jacobian -------------------------------------------------------------------------
    ! M0
    dfdu(1,1) = -1.0D0 * beta1_scaled*(EA + EF) - beta2 - deltaM + etaM*(1.0D0 - 2.0D0*(M0 + M1 + M2)) - (thetaM * KA)/(CKA_scaled + KA)
    dfdu(1,2) = etaM*(1.0D0 - 2.0D0*(M0 + M1 + M2))
    dfdu(1,3) = etaM*(1.0D0 - 2.0D0*(M0 + M1 + M2))
    dfdu(1,4) = 0.0D0
    dfdu(1,5) = -1.0D0 * (thetaM * CKA_scaled * M0)/((CKA_scaled + KA) * (CKA_scaled + KA))
    dfdu(1,6) = 0.0D0
    dfdu(1,7) = -1.0D0 * beta1_scaled * M0
    dfdu(1,8) = -1.0D0 * beta1_scaled * M0
    ! M1
    dfdu(2,1) = beta1_scaled*(EA + EF) + (thetaM * KA)/(CKA_scaled + KA)
    dfdu(2,2) = -1.0D0 * beta12 - deltaM
    dfdu(2,3) = beta21
    dfdu(2,4) = 0.0D0 
    dfdu(2,5) = (thetaM * CKA_scaled * M0)/((CKA_scaled + KA) * (CKA_scaled + KA))
    dfdu(2,6) = 0.0D0
    dfdu(2,7) = beta1_scaled * M0
    dfdu(2,8) = beta1_scaled * M0
    ! M2
    dfdu(3,1) = beta2
    dfdu(3,2) = beta12
    dfdu(3,3) = -1.0D0 * beta21 - deltaM
    dfdu(3,4) = 0.0D0
    dfdu(3,5) = 0.0D0
    dfdu(3,6) = 0.0D0
    dfdu(3,7) = 0.0D0
    dfdu(3,8) = 0.0D0
    ! K0
    dfdu(4,1) = 0.0D0
    dfdu(4,2) = -1.0D0 * (thetaK * CM1_scaled * K0)/((CM1_scaled + M1) * (CM1_scaled + M1))
    dfdu(4,3) = 0.0D0
    dfdu(4,4) = -1.0D0 * deltaK - (thetaK * M1)/(CM1_scaled + M1)
    dfdu(4,5) = 0.0D0
    dfdu(4,6) = 0.0D0
    dfdu(4,7) = 0.0D0
    dfdu(4,8) = 0.0D0
    ! KA
    dfdu(5,1) = 0.0D0
    dfdu(5,2) = (thetaK * CM1_scaled * K0)/((CM1_scaled + M1) * (CM1_scaled + M1))
    dfdu(5,3) = 0.0D0
    dfdu(5,4) = -1.0D0 * KA * etaK + (thetaK * M1)/(CM1_scaled + M1)
    dfdu(5,5) = -1.0D0 * sigma_scaled * (EA + EF) - deltaK + etaK*(1 - K0 - 2*KA)
    dfdu(5,6) = 0.0D0
    dfdu(5,7) = -1.0D0 * KA * sigma_scaled
    dfdu(5,8) = -1.0D0 * KA * sigma_scaled
    ! E0
    dfdu(6,1) = 0.0D0
    dfdu(6,2) = 0.0D0
    dfdu(6,3) = 0.0D0
    dfdu(6,4) = 0.0D0
    dfdu(6,5) = 0.0D0
    dfdu(6,6) = -1.0D0 * deltaE0
    dfdu(6,7) = 0.0D0
    dfdu(6,8) = 0.0D0
    ! EF
    dfdu(7,1) = 0.0D0
    dfdu(7,2) = -1.0D0 * omega_scaled * (1-gamma) * EF
    dfdu(7,3) = 0.0D0
    dfdu(7,4) = 0.0D0
    dfdu(7,5) = -1.0D0 * gamma * omega_scaled * EF
    dfdu(7,6) = rho0 * deltaE0
    dfdu(7,7) = -1.0D0 * omega_scaled * ((1-gamma) * M1 + gamma * KA) - deltaE - rhoF
    dfdu(7,8) = 0.0D0
    ! EA
    dfdu(8,1) = 0.0D0
    dfdu(8,2) = -1.0D0 * omega_scaled * (1-gamma) * EA
    dfdu(8,3) = (CM2_scaled * etaE * EA * (1.0D0 - EA)) / ((CM2_scaled + M2) * (CM2_scaled + M2))
    dfdu(8,4) = 0.0D0
    dfdu(8,5) = -1.0D0 * EA * gamma * omega_scaled
    dfdu(8,6) = 0.0D0
    dfdu(8,7) = rhoF
    dfdu(8,8) = -1.0D0 * omega_scaled * ((1-gamma) * M1 + gamma * KA) - deltaE + (M2 * etaE * (1.0D0 - 2.0D0 * EA)) / (CM2_scaled + M2)

    IF(ijac.EQ.1)RETURN

end subroutine FUNC

! -------------------------------------------------------------------------------------------------

! STPNT: Sets the initial conditions for the variables and the parameters.

! Input arguments :
!      ndim   :   Dimension of the algebraic or ODE system 

! Values to be returned :
!      u      :   A starting solution vector
!      par    :   The corresponding equation-parameter values

! Note : For time- or space-dependent solutions this subroutine has
!        the scalar input parameter t contains the varying time or space
!        variable value.

subroutine STPNT(ndim, u, par, t)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndim
    DOUBLE PRECISION, INTENT(INOUT) :: u(ndim),par(*)
    DOUBLE PRECISION, INTENT(IN) :: t

    ! For rho0 = 0.1, stable high
    ! par(1) = -7.5     ! Initial value for beta1_scaled
    ! par(2) = -5.5     ! Initial value for omega_scaled
  
    ! u(1) = 0.0025987    ! M0
    ! u(2) = 0.9474       ! M1
    ! u(3) = 0.0024922    ! M2 
    ! u(4) = 0.0053585    ! K0
    ! u(5) = 1.4678e-06   ! KA 
    ! u(6) = 1.26e-6      ! E0
    ! u(7) = 1.5014e-07   ! EF
    ! u(8) = 0.23115      ! EA

    ! For rho0 = 0.2, stable high
    par(1) = -7.5     ! Initial value for beta1_scaled
    par(2) = -5.5     ! Initial value for omega_scaled
  
    u(1) = 0.0025987    ! M0
    u(2) = 0.9474       ! M1
    u(3) = 0.0024922    ! M2 
    u(4) = 0.0053585    ! K0
    u(5) = 1.4678e-06   ! KA 
    u(6) = 1.26e-6      ! E0
    u(7) = 3.0029e-07   ! EF
    u(8) = 0.23115      ! EA

    ! For rho0 = 0.01, stable high
    ! par(1) = -7.5     ! Initial value for beta1_scaled
    ! par(2) = -5.5     ! Initial value for omega_scaled
  
    ! u(1) = 0.0025987    ! M0
    ! u(2) = 0.9474       ! M1
    ! u(3) = 0.0024922    ! M2 
    ! u(4) = 0.0053585    ! K0
    ! u(5) = 1.4678e-06   ! KA 
    ! u(6) = 1.26e-6      ! E0
    ! u(7) = 1.5014e-08   ! EF
    ! u(8) = 0.23115      ! EA

    ! For rho0 = 0.5, stable high
    ! par(1) = -7.5     ! Initial value for beta1_scaled
    ! par(2) = -5.5     ! Initial value for omega_scaled
  
    ! u(1) = 0.0025987    ! M0
    ! u(2) = 0.9474       ! M1
    ! u(3) = 0.0024922    ! M2 
    ! u(4) = 0.0053585    ! K0
    ! u(5) = 1.4678e-06   ! KA 
    ! u(6) = 1.26e-6      ! E0
    ! u(7) = 7.5072e-07   ! EF
    ! u(8) = 0.23115      ! EA




end subroutine STPNT

! -------------------------------------------------------------------------------------------------

!     BCND: Defines the boundary conditions. Since this is an ordinary differential 
!           equation (ODE) problem, there are no boundary conditions, so this is left empty.

!     ICND: Defines any integral constraints. For this problem, there are no integral 
!           constraints, so this routine remains empty as well.

!     FOPT: Defines an optional function that can be used for optimization or additional 
!           output. This is not needed for our system, so we leave it as a placeholder.

!     PVLS: This routine can be used to monitor certain variables or parameters during the 
!           continuation process.

subroutine BCND(ndim, par, t, u0, u1, ijac, fb, dbc)
end subroutine BCND
  
subroutine ICND(ndim, par, t, u, udot, ijac, fi, dfi)
end subroutine ICND
  
subroutine FOPT(ndim, u, icp, par, ijac, fs, dfdu, dfdp)
end subroutine FOPT
  
subroutine PVLS(ndim, u, par)
    implicit none
    integer ndim, i
    double precision u(ndim), par(*)
    integer :: unit_number
    character(len=30) :: filename

    RETURN

    ! Specify the file name
    filename = 'output.txt'

    ! Assign a unit number for the file
    unit_number = 10  ! You can use any valid number not already used

    ! Open the file for writing
    open(unit=unit_number, file=filename, status='old', position='append', action='write')

    ! Write the two parameters
    write(unit_number, '(ES12.5, A, ES12.5, A)', advance='no') par(1), ', ', par(2), ', '

    ! Loop over the array u and write values to the file, each followed by a comma
    do i = 1, ndim
        if (i < ndim) then
            write(unit_number, '(ES12.5, A)', advance='no') u(i), ', '
        else
            write(unit_number, '(ES12.5)') u(i)  ! No comma after the last value
        end if
    end do


    ! Close the file
    close(unit_number)

    ! print *, 'Text has been written to ', filename

    ! Output current values of the state variables S and I for monitoring
    ! write(*,*) 'Current value: EA = ', u(8), 'beta1_scaled = ', par(1), 'omega_scaled = ', par(2)
end subroutine PVLS
  