module Axial_Asymmetry_PES_mod
    use iso_fortran_env, only: real64
    use constant_mod
    use common_func_mod, only: LegendreP, elliptical_integral_1st, conf_hyper_func
    use nucleus_mod, only: nucleus_property
    implicit none
    private

    integer, parameter :: dp = real64
    ! integer, parameter, public :: Z_num = 10
    ! integer, parameter, public :: N_num = 10
    ! integer, public :: A_num = Z_num + N_num
    ! real(dp), parameter, public :: pi = 3.1415926535897932384626433832795_dp
    ! real(dp), parameter, public :: hbar_c = 197.3269804_dp ! MeV*fm
    ! real(dp), parameter, public :: hbar = 6.582119569e-22_dp  ! MeV*s
    ! real(dp), parameter, public :: m_nucleon = 938.90595_dp   ! MeV/c^2
    ! real(dp), parameter, public :: e2 = 1.4399764_dp          ! MeV*fm
    ! real(dp), parameter, public :: r0 = 1.16_dp ! fm
    type :: axial_asymmetry_PES_type
        real(dp) :: epsilon
        real(dp) :: epsilon2
        real(dp) :: epsilon3
        real(dp) :: epsilon22
        real(dp) :: epsilon4
        real(dp) :: epsilon6
        real(dp) :: gamma
        real(dp) :: energy
        real(dp) :: k
        real(dp) :: phi_0
        real(dp) :: d_lambda
        real(dp) :: hbar_omega_x
        real(dp) :: hbar_omega_y
        real(dp) :: hbar_omega_z
        real(dp) :: hbar_0omega0 
        real(dp) :: hbar_omega0
        real(dp) :: radius
        real(dp) :: semi_a, semi_b, semi_c

        character(len=20) :: particle_type = "proton" ! or "neutron"

        ! Nilsson model parameters
        ! real(dp) :: hbar_omega_z
        real(dp) :: hbar_omega_perp
        real(dp) :: xi
        real(dp) :: eta
        real(dp) :: zeta
        real(dp) :: rho
        real(dp) :: cos_theta_t

        real(dp), allocatable :: wave_func(:,:,:)

        

        
        contains
            procedure :: set
            procedure :: potential
    end type axial_asymmetry_PES_type

    public :: axial_asymmetry_PES_type


    contains
        

        subroutine set(this,nucleus)
            class(axial_asymmetry_PES_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp) :: temp, d_x, sum, x
            integer :: n, n_max
            n_max = 100000
            ! This subroutine can be used to generate a grid of deformation parameters and calculate the corresponding energy values.
            ! The actual implementation will depend on the specific model being used for the energy calculation.
            this%hbar_0omega0 = 41.0_dp*nucleus%A**(-1.0_dp/3.0_dp)
            this%epsilon2 = this%epsilon*cos(this%gamma)
            this%epsilon22 = this%epsilon*sin(this%gamma)*2.0_dp

            ! integral from -1 to 1
            sum = 0.0_dp
            do n = 0, n_max - 1
                d_x = 2.0_dp / n_max
                x = -1.0_dp + n * d_x
                temp = 1 - this%epsilon*LegendreP(2,x) + 2.0_dp*this%epsilon3*LegendreP(3,x) &
                + 2.0_dp*this%epsilon4*LegendreP(4,x) + 2.0_dp*this%epsilon6*LegendreP(6,x)
                sum = sum + temp**(-3.0_dp/2.0_dp)*d_x*0.5_dp
            end do


            if (this%particle_type == "neutron") then
                this%hbar_0omega0 = this%hbar_0omega0 * (1.0_dp + 1.0_dp/3.0_dp* (nucleus%N - nucleus%Z)/nucleus%A)
            else if (this%particle_type == "proton") then
                this%hbar_0omega0 = this%hbar_0omega0 * (1.0_dp - 1.0_dp/3.0_dp* (nucleus%N - nucleus%Z)/nucleus%A)
            end if

            this%hbar_omega0 = this%hbar_0omega0 * sum**(1.0_dp/3.0_dp)/(1.0_dp + this%epsilon/3.0_dp)&
            /(1.0_dp - 2.0_dp*this%epsilon/3.0_dp )**(1.0_dp/2.0_dp)

            this%hbar_omega_x = this%hbar_omega0 * (1 + this%epsilon2/3.0_dp + 0.5_dp*this%epsilon22)
            this%hbar_omega_y = this%hbar_omega0 * (1 + this%epsilon2/3.0_dp - 0.5_dp*this%epsilon22)
            this%hbar_omega_z = this%hbar_omega0 * (1 - 2.0_dp*this%epsilon2/3.0_dp)

            ! this%radius = nucleus%R_a
            ! this%semi_a = this%radius * hbar / (m_nucleon * this%hbar_omega_x)**0.5_dp&
            ! *(m_nucleon * this%hbar_omega0 / hbar**2)**(1.0_dp/2.0_dp)
            ! this%semi_b = this%radius * hbar / (m_nucleon * this%hbar_omega_y)**0.5_dp&
            ! *(m_nucleon * this%hbar_omega0 / hbar**2)**(1.0_dp/2.0_dp)
            ! this%semi_c = this%radius * hbar / (m_nucleon * this%hbar_omega_z)**0.5_dp&
            ! *(m_nucleon * this%hbar_omega0 / hbar**2)**(1.0_dp/2.0_dp)

        end subroutine set

        subroutine potential(this,nucleus)
            implicit none
            class(axial_asymmetry_PES_type), intent(inout) :: this
            type(nucleus_property), intent(in) :: nucleus
            real(dp) :: k, phi_0
            integer :: n, n_max
            n_max = 100000
            k = sqrt((this%semi_b**2 - this%semi_a**2) / (this%semi_c**2 - this%semi_a**2))
            print *, "k: ", k
            phi_0 = asin(sqrt((this%semi_c**2 - this%semi_a**2) )/this%semi_c)
            print *, "phi_0: ", phi_0

            this%energy = 0.6_dp*(nucleus%Z**2*e2)/sqrt(this%semi_c**2 - this%semi_a**2)&
            *elliptical_integral_1st(k, phi_0)
            
        end subroutine potential
    

end module Axial_Asymmetry_PES_mod

program name
    use iso_fortran_env, only: real64
    use constant_mod
    use common_func_mod, only: LegendreP, elliptical_integral_1st, conf_hyper_func
    use Axial_Asymmetry_PES_mod
    implicit none
    type(axial_asymmetry_PES_type) :: pes
    integer, parameter :: dp = real64
    real(dp) :: energy_old

    pes%epsilon = 0.1_dp
    pes%gamma = 10.0_dp * pi / 180.0_dp !

    ! call pes%set()
    ! print *, pes%semi_a*pes%semi_b*pes%semi_c
    ! print *, pes%radius**3
    ! call pes%potential()
    ! print *, "Energy: ", pes%energy, " MeV"
    ! energy_old = pes%energy
    ! pes%epsilon = 0.0_dp
    ! pes%gamma = 10.0_dp * pi / 180.0_dp !
    ! call pes%set()
    ! call pes%potential()
    ! print *, "Energy: ", pes%energy, " MeV"
    ! print *, "Ratio-1: ", (energy_old - pes%energy)/pes%energy

    
end program name