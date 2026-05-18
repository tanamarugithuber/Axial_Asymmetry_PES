module common_func_mod
    use iso_fortran_env, only: real64
    implicit none
    private
    integer, parameter :: dp = real64
    public :: LegendreP, elliptical_integral_1st, conf_hyper_func

contains

function LegendreP(n,x) result(solution)
            implicit none
            integer, intent(in) :: n
            real(dp), intent(in) :: x
            real(dp) :: solution, theta
            integer :: i
            if (n == 0) then
                solution = 1.0_dp
            else if (n == 1) then
                solution = x
            else if (n == 2) then
                solution = 0.5_dp*(3.0_dp*x**2 - 1.0_dp)
            else if (n == 3) then
                solution = 0.5_dp*(5.0_dp*x**3 - 3.0_dp*x)
            else if (n == 4) then
                solution = (35.0_dp*x**4 - 30.0_dp*x**2 + 3.0_dp)*0.125_dp
            else if (n == 5) then
                solution = (63.0_dp*x**5 - 70.0_dp*x**3 + 15.0_dp*x)*0.125_dp
            else if (n == 6) then
                solution = (231.0_dp*x**6 - 315.0_dp*x**4 + 105.0_dp*x**2 - 5.0_dp)*0.0625_dp
            else 
                theta = acos(x)
                solution = 0.0_dp
                do i = 0, n
                    solution = solution + cos((n - 2*i)*theta)*gamma(real(2*n -2*i, dp ))*gamma(real(2*i, dp )) &
                    / (gamma(real(n - i, dp ))**2*gamma(real(i, dp ))**2*2.0_dp**(2*n))
                end do
            end if
        end function LegendreP

        function elliptical_integral_1st(x,k) result(solution)
            implicit none
            real(dp), intent(in) :: x, k
            real(dp) :: solution
            integer :: n, n_max
            real(dp) ::  d_x, t
            n_max = 100000
            d_x = x / n_max
            ! This function can be used to calculate the elliptic integral of the first kind, which is often used in nuclear physics calculations.
            ! The actual implementation of the function will depend on the specific model being used for the energy calculation.
            solution = 0.0_dp
            do n = 0, n_max - 1
                t = n * d_x
                solution = solution + 1.0_dp / sqrt(1.0_dp - t**2)&
                /sqrt(1.0_dp - k**2 * t**2) * d_x
            end do
        end function elliptical_integral_1st

        function conf_hyper_func(a,c,x) result(solution)
            implicit none
            real(dp), intent(in) :: a, c, x
            real(dp) :: solution
            real(dp) :: pi_a, pi_c, n_factorial, temp
            integer :: i, n
            ! This function can be used to calculate the energy of the nucleus based on the deformation parameters.
            ! The actual implementation of the function will depend on the specific model being used for the energy calculation.
            solution = 1.0_dp ! Placeholder for the actual energy calculation
            n_factorial = 1.0_dp
            
            pi_a = 1.0_dp
            pi_c = 1.0_dp
            do n = 1, 100000
                n_factorial = n_factorial * n
                pi_a = pi_a * (a + n - 1.0_dp) 
                pi_c = pi_c * (c + n - 1.0_dp)
                temp = pi_c * n_factorial
                solution = solution + pi_a * x**n / temp
            end do
        end function conf_hyper_func

end module common_func_mod
