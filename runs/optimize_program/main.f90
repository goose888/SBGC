
  program opt_prog

  USE userfunc_mod
  
  implicit none

  external :: grobfd, grcnfd

  bl(1) = 1.0d0
  bl(2) = 1.0d8
  bl(3) = 3.1536d12

  bu(1) = 2.0d0
  bu(2) = 1.0d16
  bu(3) = 3.1536d13

  ! For site 197
  x(1)  = 1.5d0
  x(2)  = 1.0d13
  x(3)  = 1.0d13

  call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,     &
             miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,   &
             f,g,iw,iwsize,w,nwsize,obj,constr,grobfd,grcnfd)

 
  end program opt_prog
