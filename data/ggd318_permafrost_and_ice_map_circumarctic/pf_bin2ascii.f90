program pf_bin2ascii
!     sample program to read permafrost arrays

      implicit none

      integer, parameter :: width=1441, height=1441

!     Nh grid dimensions
!!!      parameter (width = 1441, height = 1441)

!     Nl grid dimensions
!!!      parameter (width = 721, height = 721)

!     1/2 x 1/2 deg grid dimensions
!!!      parameter (width = 720, height = 360)

      character*1 :: pf(width, height)
      integer :: row, col  !, x, y

!     open file and read data
!     open(unit=1, file='NhIPA.byte', form='system')
      open(unit=1, file='nhipa.byte', form='binary')
      read(1) pf
      close(1)
      
      print *, 'test:', pf(22,22)

!     query grid
! 100  continue

!     x,y start at 0,0 in lower left corner
!     rdpix is an IDL routine used for testing

!      print *, 'enter x y from rdpix'
!      read *, x, y
!      col = height-y
!      row = x+1

!     write out ascii file
      open(unit=12, file='pf.ascii', status='unknown')
      do col=1,height
         do row=1,width
           ! write(12, "(i6,i6,i6)")  col, row, ichar(pf(row, col))
            write(12, "(i6)") ichar(pf(row, col))
         end do
      enddo
      close(10)

!      goto 100

end
