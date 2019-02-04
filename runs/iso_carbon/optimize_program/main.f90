
  program opt_prog

  USE userfunc_mod
  
  implicit none

  external :: grobfd, grcnfd

  bl(1) = 2.0d-1
  bl(2) = 5.754d-7
  bl(3) = 1
  bl(4) = 1.5
  bl(5) = 2

  bu(1) = 3.0d0
  bu(2) = 3.18d-5
  bu(3) = 40
  bu(4) = 2.1
  bu(5) = 8

  ! ! For site 110
  ! x(1)  = 1.53d0
  ! x(2)  = 1.918d-6
  ! x(3)  = 3.86
  ! x(4)  = 1.6
  ! x(5)  = 4.0

  ! ! For site 43
  ! x(1)  = 2.984839
  ! x(2)  = 5.7540001E-07
  ! x(3)  = 19.91537
  ! x(4)  = 2.080639
  ! x(5)  = 2.036407

  ! ! For site 143
  ! x(1)  = 0.2700000
  ! x(2)  = 5.7540001E-06
  ! x(3)  = 1.410000
  ! x(4)  = 1.600000
  ! x(5)  = 4.009000

  ! ! For site 146
  ! x(1)  = 1.530000
  ! x(2)  = 2.7325099E-05
  ! x(3)  = 1.410000
  ! x(4)  = 1.500003
  ! x(5)  = 2.001000

  ! For site 197
  x(1)  = 1.530000
  x(2)  = 5.7540001E-07
  x(3)  = 14.00000
  x(4)  = 1.600000
  x(5)  = 2.001000

  call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,     &
             miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,   &
             f,g,iw,iwsize,w,nwsize,obj,constr,grobfd,grcnfd)

 
  end program opt_prog
