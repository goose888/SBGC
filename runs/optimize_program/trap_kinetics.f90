program trap_kinetics
   implicit none

   integer :: i
   real, parameter :: D_dot = 9.51294E-11
   real, parameter :: D0 = 3000.
 !  real, dimension(6), parameter :: Temp = (/315.15, 319.15, 326.15, 330.15, 340.15, 345.15/)
   real, dimension(5), parameter :: Temp = (/313.75, 326.15, 330.15, 340.15, 345.15/)

   real :: lamda_prime
   real :: lamda
   real, parameter :: t = 3.1536E13
   real, parameter :: k = 0.00008617
 !  real, dimension(6), parameter :: t_a_obs = (/3.58497E+12, 5.90848E+12, 2.79899E+12, 1.96274E+12, 1.10705E+12, 6.63284E+11/)
   real, dimension(5), parameter :: t_a_obs = (/5.90848E+12, 2.79899E+12, 1.96274E+12, 1.10705E+12, 6.63284E+11/)

 !  real, dimension(6) :: t_a   ! Apparent age
   real, dimension(5) :: t_a   ! Apparent age
   real :: obj  ! Objective function
   real :: E    ! Paramter 1
   real :: S    ! Paramter 2
   real :: accu_sqrterr

   ! Read in parameters
   open(52, file='tkpar', status='unknown')
      read(52,*) E, S
   close(52)

   print *, 'E=', E, 'S=', S

   ! Model
   ! Calculate the apparent age
   accu_sqrterr = 0.
   do i = 1, 5
      lamda = 1. / S * exp(E/(k*Temp(i)))
     ! t_a(i) = lamda * (1 - exp((-t)/lamda))
      lamda_prime = 1. / (D_dot / D0 + 1/lamda)
      t_a(i) = - D0 / D_dot * log(1 - D_dot/D0 * lamda_prime * (1 - exp(-t/lamda_prime)))
      accu_sqrterr = accu_sqrterr + (t_a(i) - t_a_obs(i)) * (t_a(i) - t_a_obs(i))/5.
   end do

   print *, 't_a:', t_a(1)
  ! print *, 't_a_obs:', t_a_obs

   ! Calculate the objective function
   obj = sqrt(accu_sqrterr)

   print *, 'obj:', obj

   ! Write out the objective function result
   open(53, file='tk_obj.dat', status='unknown')
      write(53,*) obj
   close(53)

end program trap_kinetics
