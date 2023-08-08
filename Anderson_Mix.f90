module Anderson_Mix
  ! This file was provided by Dr. Bill Atkinson, Trent U.
  implicit none
  PRIVATE
  PUBLIC :: AM_startup,AM,AM_shutdown
  integer, parameter :: M=4
  integer, parameter :: dbl=8  
  integer, allocatable :: IPIV(:)
  real (dbl), allocatable, save :: x_old(:,:),F_old(:,:),F(:)
  contains

    subroutine AM_startup(N)
      integer, intent (in) :: N

      if (.not.allocated(x_old)) then
         print *, "allocating workspace"
         allocate(x_old(N,M),F_old(N,M),F(N),IPIV(N))
      else
         print *,"Workspace Already Allocated.  Doing Nothing."
      end if
    end subroutine AM_startup

      
    subroutine AM_shutdown()
      if (allocated(x_old)) then
         print *, "deallocating workspace"
         deallocate(x_old,F_old,F,IPIV)
      else
         print*,"Workspace Alread Cleared.  Doing Nothing."
      end if
    end subroutine AM_shutdown

    subroutine AM(x,y,x_new)
      real (dbl), intent (in) :: x(:),y(:)
      real (dbl), intent (out) :: x_new(:)

      integer :: N
      integer, save :: iteration = 0, head = 0
      integer :: ir,ic,is,it
      real (dbl) :: simple_mix1=0.001_dbl, simple_mix2=0.001_dbl
      real (dbl) :: A(M,M),B(M)

      integer :: info
      
      N=size(x)
      if (.not.allocated(x_old)) then
         call AM_startup(N)
      end if

      iteration=iteration+1
      F = y-x
      head = mod(head,M)+1
      if (iteration<=M+1) then
         ! use simple mixing
         x_old(:,head)=x
         F_old(:,head)=F
         x_new = x + simple_mix1*F
      else
         do ir=1,M
            is=mod(head+ir-2,M)+1
            B(ir) = dot_product(F-F_old(:,is),F)
            do ic=1,M
               it=mod(head+ic-2,M)+1
               A(ir,ic) = dot_product( (F-F_old(:,is)),(F-F_old(:,it)) )
            end do
            A(ir,ir) = 1.001*A(ir,ir) ! do this to increase stability
         end do
         call dgesv(M,1,A,M,IPIV,B,M,INFO)
         if (info/=0) then
            print *,"Error with DGESV in subroutine AM"
            print *,"Info = ",info
            stop
         end if

         x_new = x + simple_mix2*F 
         do ir=1,M
            is=mod(head+ir-2,M)+1
            x_new = x_new + B(ir)*(x_old(:,is)-x + simple_mix2*(F_old(:,is)-F))
         end do
         F_old(:,head) = F
         x_old(:,head) = x
      end if

    end subroutine AM
    


end module Anderson_Mix
