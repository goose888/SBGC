! This is file to set inputs varibles for FFSQP and interface with program 
! 
! Min Xu @ISWS
! Shijie Shu for ISAM-1DSBGC calibration
!

  module userfunc_mod
    implicit none

    integer, parameter :: r8=kind(1.0d0)

!
! user should change the following parameters to fits their application
!
    integer, parameter, public ::     & !
                     nparam  =    5,  & !
                     nf      =    1,  & !
                     neqn    =    0,  & !
                     nineqn  =    0,  & !
                     nineq   =    0,  & !
                     neq     =    0,  & !
                     mode    =  100,  & !
                     iprint  =    2,  & !
                     miter   =  1000     !

    !integer, public ::     & !
    !                 iprint  =    0     !

    integer, parameter, public ::     & !
                     nfsize  = max(1,nf), &
                     ngsize  = max(1,nineq+neq), &
                     iwsize  = 6*nparam+8*max(1,nineq+neq)+7*max(1,nf)+30, & !
                     nwsize  = 4*nparam*nparam+5*max(1,nineq+neq)*nparam+3*max(1,nf)*nparam &
                               +26*(nparam+max(1,nf))+45*max(1,nineq+neq)+100


    real(r8), public, save ::             & !
                     bigbnd = 1.0d+20,   & !
                     eps    = 1.0d-15,   & !
!                     eps    = 1.0d+4,   & !
                     epseqn = 1.0d-15,   & ! not effects in this case
                     udelta = 1.0d-1       ! have no idea to choose

! ---------------------------------------------------------------------------------
! defined in ffsqp

   integer, public :: iw(iwsize)

   real(r8), public  ::  x(nparam), & !
                        bl(nparam), & !
                        bu(nparam), & !
                         f(nfsize), & !
                         g(ngsize), & !
                         w(nwsize)

!   external obj,constr,gradob,gradcn
! ---------------------------------------------------------------------------------

!
! the following needn't be changed, the values of them will be set in run time
!
    integer, public  :: inform  

!
! following is set to get weight function from different schemes
!

    real(r8), public :: var_pkm,var_pgm,cov_pkg

  contains

! ---------------------------------
    subroutine obj(nparam,j,x,fj)
    implicit none

    integer  :: nparam, j, lb, ip, ik
    integer  :: p
    real(r8) :: x(nparam),fj
    real     :: diff_soc
    real     :: diff_d14c
    real     :: max_cryo_depth, d_const, cryo_tao, defac_t_q10, defac_z_s, f_damp_bd, biome_bb, biome_alpha
    character(len=80) :: cmd

    fj = 0.0_r8

    cmd = './isam < namelist.hist > /dev/null'

    ! Set the new parameters and overwrite the parameter input file
    open(17, file='param_sbgc.in', status='old')
    read(17, *) max_cryo_depth, d_const, cryo_tao, defac_t_q10, defac_z_s, f_damp_bd, biome_bb, biome_alpha
    close(17)

    call system('rm param_sbgc.in')

    max_cryo_depth = x(1)
    d_const = x(2)
    cryo_tao = x(3)
    defac_t_q10 = x(4)
    defac_z_s = x(5)

    open(17, file='param_sbgc.in', status='new')
    write(17, *) max_cryo_depth
    write(17, *) d_const
    write(17, *) cryo_tao
    write(17, *) defac_t_q10
    write(17, *) defac_z_s
    write(17, *) f_damp_bd
    write(17, *) biome_bb
    write(17, *) biome_alpha
    close(17)

    print*, 'Running proc...', trim(cmd)
    print*, 'Parameters:', x(1), x(2), x(3), x(4), x(5)
    call system(trim(cmd))
    open(20, file='difference_d14c.dat', status='old')
       read(20,*) fj
    close(20)

    print *, 'obj func1 (RMSE):', fj !diff, stn(ik), spn(ip)

    return
    end subroutine obj

!------------------------------------
    subroutine constr(nparam,j,x,gj)
    integer :: nparam,j
    real(r8) ::  x(nparam), gj

    return
    end subroutine constr

!-----------------------------------------------

  end module userfunc_mod
