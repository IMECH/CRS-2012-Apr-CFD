! ==============================================================
!   One-dimensional shock tube solved by Lax and leapflog scheme
! --------------------------------------------------------------
! At t = 0:
!       +----------------------
!       |                     |
!       |    rho1  u1  p1     |    rho2  u2  p2     
!       |                     -----------------------
!       +---------------------+---------------------+--------> X
!       0                    L/2                    L
!
! Solve the Riemann problem. At both ends, the tube is closed
! so u = 0. Equations:
!       d(rho)  /dt + d(rho*u)      /dx = 0
!       d(rho*u)/dt + d(rho*u*u + p)/dx = 0
!       d(E)    /dt + d(u*(rho*E+p))/dx = 0
!
! ==============================================================
! !!!!!!!!!!!!!!!!!!!!!!!!!! Output !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! "solution.dat" = Data file containing for each grid point,
!                  coordinate, density, velocity, pressure, 
!                  in the following format:
!
!              x(1)        rho(1)        u(1)        p(1)
!              x(2)        rho(2)        u(2)        p(2)
!               .            .            .           .
!               .            .            .           .
!              x(n)        rho(n)        u(n)        p(n)
! 
! Use the matlab program, data2figure.m, to plot the solutions.
! ==============================================================
! !!!!!!!!!!!!!!  READ THIS BEFORE YOU START !!!!!!!!!!!!!!!!!!!
!  
! W in the program is the matrix of conservative varibles: rho, 
! rho*u and E stored as columns in the matrix W
!        _                                              _
!       |    rho_1        (rho*u)_1           E_1        |
!       |    rho_2        (rho*u)_2           E_2        |
!  W =  |      .              .                .         |  
!       |      .              .                .         |
!       |_   rho_n        (rho*u)_n           E_n       _|
!
! F in the program is the matrix of flux: rho*u, rho*u**2+p and
! u*(rho*E+p) stored as columns in the matrix F
!        _                                              _
!       |  (rho*u)_1   (rho*u**2+p)_1   (u*(rho*E+p))_1  |
!       |  (rho*u)_2   (rho*u**2+p)_2   (u*(rho*E+p))_2  |
!  F =  |    .               .               .           |
!       |    .               .               .           |
!       |_ (rho*u)_n   (rho*u**2+p)_n   (u*(rho*E+p))_n _|
!
! rho: Density,  u: Velocity,  p : Pressure
!
! ==============================================================
! !!!!!!!!!!!!!!!!!  Record of revisions !!!!!!!!!!!!!!!!!!!!!!!
! 
!        Date        Programmer        Description of change
!       4/02/12      zhou lvwen            original code
!---------------------------------------------------------------
!

Program TubeShock

implicit none
! ==============================================================
! !!!!!!!!!!!!!! Declare constant parameters !!!!!!!!!!!!!!!!!!!
integer, parameter    :: o     = 8    ! Double Precision
real(o), parameter    :: gamma = 1.4  ! Ratio of specific heats 
real(o), parameter    :: L     = 1.0  ! Length of domain
real(o), parameter    :: CFL   = 0.5  ! CFL number for stability
real(o), parameter    :: Tend  = 0.2  ! final time
integer, parameter    :: n     = 10001! Number of grid points
real(o), parameter    :: dx = L/(n-1) ! Spatial step size
! ............... Define intial conditions .....................
! Right conditions:
real(o), parameter    :: p1    = 0.1  ! Pressure in right side
real(o), parameter    :: rho1  = 0.125! Density  at right side
real(o), parameter    :: u1    = 0.0  ! Velocity in right side
! Left conditions:
real(o), parameter    :: p2    = 1.0  ! Pressure in left  side
real(o), parameter    :: rho2  = 1.0  ! Density  at left  side
real(o), parameter    :: u2    = 0.0  ! Velocity in left  side

! ==============================================================
! !!!!!!!!!! Declare variable types & definitions !!!!!!!!!!!!!!
integer               :: i            ! loop index
integer               :: nsteps= 0    ! Number of time steps
real(o)               :: t     = 0.0  ! current time
real(o)               :: dt           ! time step
real(o),dimension(n)  :: x     = [(i*dx, i = 0, n-1)]
integer,dimension(n-2):: nc    = [(i, i= 2, n-1)]
! center point index  :: nc    = [2, 3, ..., n-1]
! left   point index  :: nc-1  = [1, 2, ..., n-2]
! right  point index  :: nc+1  = [3, 4, ..., n  ]
real(o),dimension(n)  :: p            ! Pressure corresponding x
real(o),dimension(n)  :: rho          ! Density  corresponding x
real(o),dimension(n)  :: u            ! Velocity corresponding x
real(o),dimension(n)  :: E            ! energy   corresponding x
real(o),dimension(n,3):: W            ! conservative varibles
real(o),dimension(n,3):: Wbefore      ! for leapfrog
real(o),dimension(n,3):: F            ! flux
real(o)               :: eta   = 1.0  ! Art_visc coefficient
character             :: scheme*10    ! Difference scheme
integer               :: process = 0  ! process index

! Initial conditions
p  (1:(n+1)/2) = p2  ;  p  ((n+3)/2:n) = p1  ;
rho(1:(n+1)/2) = rho2;  rho((n+3)/2:n) = rho1;
u  (1:(n+1)/2) = u2  ;  u  ((n+3)/2:n) = u1  ;
E = p/((gamma-1)*rho) + 0.5*u**2

! Prompted to select a difference scheme in the terminal
call Initialization(scheme)

Do while (t<Tend)

   ! Determine time step by CFL condition
   dt = timestep(CFL, dx)

   ! code the variables
   call u2WF(rho, u, p, E, W, F)

   !Add artificial viscosity
   W = W + artificial_visc(rho, W, eta)

   ! Use lax or leapfrog to update the solution
   if (scheme =='lax') then
      call Lax(W, F, dt)
   else if (scheme=='leapfrog') then
      call leapfrog(W, Wbefore, F, dt, nsteps)
   end if

   ! boundary conditions: u(0,t) = u(L,t) = 0
   call boundary(u)
   
   ! Decode the variables
   call W2u(W, rho, u, p, E)

   ! update the time and time steps
   t = t + dt
   nsteps = nsteps + 1

   ! draw the program running process bar
   call processBar(t, Tend, process)

End do



call output(t)

! End of program
! **************************************************************
! **************************************************************

contains
! Below is a list of subroutines and functions used in the main 
! program.

! **************************************************************
subroutine Initialization(scheme)
! Print some information on the termianl about the program and
! Prompted to select a difference scheme in the terminal
implicit none
character,intent(inout):: scheme*10   ! Difference scheme

write(*,*), "------- One-dimensional shock tube problem -------" 
write(*,*), " Copyrhigt by Zhou Lvwen: zhou.lv.wen[at]gmail.com"
write(*,*), "---------------- intial conditions ---------------" 
write(*,*), '-Right conditions: '
write(*, '("   rho = ",f6.3, "  u = ",f6.3, "  p = ",f6.3)')&
                    rho1,            u1,            p1
write(*,*), '-Left  conditions: '
write(*, '("   rho = ",f6.3, "  u = ",f6.3, "  p = ",f6.3)')&
                    rho2,            u2,            p2
write(*,'(" -Length of domain:", f6.2)') , L 
write(*,*)
write(*,*), "------------- computational parameters -----------" 
write(*,*), 'The recommended parameters will be used'
write(*,'(" -CFL number for stability :", f5.2)') CFL
write(*,'(" -Number of grid points    :", i5  )') n
write(*,'(" -Final time               :", f5.2)') Tend
write(*,*)
write(*,*), "----------------- Difference scheme --------------" 
1 write(*,'(A,$)'),  &
             ' Please select Difference scheme (lax/leapfrog): '
read(*,'(A)') scheme
select case (scheme)
   case ('lax')
       write(*,*), 'You select lax scheme'
   case ('leapfrog')
       write(*,*), 'You select leap-forg scheme'
   case default
       write(*,*), 'please input "lax" or "leapfrog"'; goto 1
end select
write(*,*)

end subroutine Initialization
! --------------------------------------------------------------

! **************************************************************
function timestep(CFL, dx) result(dt)
! Determine dt by CFL condition: dt = CFL * dx / max_wavespeed
!                                     CFL <= 1.0
!
implicit none
real(o), intent(in)   :: CFL          ! CFL number for stability
real(o), intent(in)   :: dx           ! Spatial step size
real(o)               :: dt           ! time step
real(o)               :: umax         ! maximal absolute speed
real(o),dimension(n)  :: c            ! speed of sound

! Find the speed of sound
  c = sqrt(gamma*p/rho);

! Compute maximal absolute value of the characteristic speeds
  umax = maxval(abs(u+c));

! Determine the time step using the CFL condition
  dt = CFL*dx/umax
end function timestep

! --------------------------------------------------------------


! **************************************************************
function artificial_visc(rho, W, eta) result(Q)
! Artificial dissipative flux
!     Q = 1/2 * eta * sw * [ W(i+1) - 2W(i) + W(i-1) ] 
!
! where sw is density switch: 
!         |  |rho(i+1) - rho(i)| - |rho(i) - rho(i-1)|  |
!    sw = | ------------------------------------------- | 
!         |  |rho(i+1) - rho(i)| + |rho(i) - rho(i-1)|  |    
!

implicit none
real(o),intent(in),dimension(n)  :: rho 
real(o),intent(in),dimension(n,3):: W 
real(o),intent(in)               :: eta! Art_visc coefficient
real(o),dimension(n,3)           :: Q  ! Artificial viscidity
real(o),dimension(n-2)           :: dur, dul
real(o),dimension(n-2,3)         :: sw ! density switch


dur = dabs( rho(nc+1)-rho(nc  ) )
dul = dabs( rho(nc  )-rho(nc-1) )
sw(:,1) = dabs((dur - dul) / (dur+dul + 1e-5))
sw(:,2) = sw(:,1);  sw(:,3) = sw(:,1)
Q = 0
Q(nc, :) = 0.5*sw*eta*(W(nc+1,:)-2*W(nc,:)+W(nc-1, :))
end function artificial_visc

! --------------------------------------------------------------


! **************************************************************
subroutine u2WF(rho, u, p, E, W, F)
! Code the variables: 
!                     rho, u, p, E    ==>    W, F
implicit none

real(o),intent(in)    :: p(n)         ! Pressure corresponding x
real(o),intent(in)    :: rho(n)       ! Density  corresponding x
real(o),intent(in)    :: u(n)         ! Velocity corresponding x
real(o),intent(in)    :: E(n)         ! energy   corresponding x
real(o),intent(out)   :: W(n,3)       ! conservative varibles
real(o),intent(out)   :: F(n,3)       ! flux

! Code the variables
  W(:,1) = rho  ;  W(:,2) = rho*u       ;  W(:,3) = rho*E
  F(:,1) = rho*u;  F(:,2) = rho*u**2 + p;  F(:,3) = u*(rho*E+p)

end subroutine u2WF
! --------------------------------------------------------------

! **************************************************************
subroutine W2u(W, rho, u, p, E)
! Decode the variables:
!                       W    ==>   rho, u, p, E
implicit none
real(o),intent(in)    :: W(n,3)       ! conservative varibles
real(o),intent(out)   :: p(n)         ! Pressure corresponding x
real(o),intent(out)   :: rho(n)       ! Density  corresponding x
real(o),intent(out)   :: u(n)         ! Velocity corresponding x
real(o),intent(out)   :: E(n)         ! energy   corresponding x

  rho = W(:,1)
  u   = W(:,2)/rho
  E   = W(:,3)/rho
  p   = (gamma-1)*rho*(E-0.5*u**2)
end subroutine W2u
! --------------------------------------------------------------

! **************************************************************
subroutine Lax(W, F, dt)
! Lax scheme for system of 1D conservation laws:
!   
!  W(j,n+1) = + 1/2*  [W(j+1,n) + W(j-1,n)] 
!             - 1/2*r*[F(j+1,n) + F(j-1,n)]   
!  r = dt/dx
!
implicit none
real(o),intent(in)    :: F(n,3)       ! flux
real(o),intent(in)    :: dt           ! time step
real(o),intent(inout) :: W(n,3)       ! conservative varibles


  W(nc,:) = 0.5*( W(nc+1,:)+W(nc-1,:) )                     &
                           - dt/(2*dx) * ( F(nc+1,:)-F(nc-1,:) )

end subroutine Lax
! --------------------------------------------------------------

! **************************************************************
subroutine leapfrog(W, Wbefore, F, dt, nsteps)
! Lax scheme for system of 1D conservation laws:
!   
!  W(j,n+1) = W(j, n-1) - r*[F(j+1,n) - F(j-1,n)]   
!  r = dt/dx
!
implicit none
real(o),intent(in)    :: F(n,3)       ! flux
real(o),intent(in)    :: dt           ! time step
integer,intent(in)    :: nsteps       ! Number of time steps
real(o),intent(inout) :: Wbefore(n,3) ! conservative varibles
real(o),intent(inout) :: W(n,3)       ! conservative varibles
real(o)               :: temp(n,3)    ! temp of W_{t}

if (nsteps == 0) then 
   Wbefore = W
   call Lax(W, F, dt)
else
  temp = W
  W(nc,:) = Wbefore(nc,:) - dt/dx * ( F(nc+1,:)-F(nc-1,:) )
  Wbefore = temp
end if

end subroutine leapfrog
! --------------------------------------------------------------

! **************************************************************
subroutine boundary(u)
! boundary conditions for the shock tube. should be implemented
! The boundary conditions are: 
!        u(0,t)=u(L,t) = 0  => rho*u = 0 at both ends
! rho should be extrapolated at both ends
implicit none
real(o),intent(inout) :: u(n)         ! Velocity corresponding x

  u([1, n]) = 0;

end subroutine boundary
! --------------------------------------------------------------


! **************************************************************
subroutine output(t)
!  Write output data file
!
!  Output:  Data file "solution.dat" containing for each point 
!           the following:
!           coordinate, density, velocity, pressure
implicit none
real(o), intent(in)   :: t            ! current time

open(unit = 8, file='solution.dat', status = 'replace')
write(8, '("%%%%%%%%", " 1D shock tube problem ", "%%%%%%%%%")')
write(8,*) '%'
write(8, '("% % % % % % % ", " conditions ", " % % % % % % %")')
write(8,'("% Right conditions            :")')
write(8,'("%   rho = ", f5.2, "  u = ", f5.2, "  p = ", f5.2)')&
                   rho1,            u1,            p1
write(8,'("% Left  conditions            :")')
write(8,'("%   rho = ", f5.2, "  u = ", f5.2, "  p = ", f5.2)')&
                   rho2,            u2,            p2
write(8,'("% Length of domain            :", f5.2)') L
write(8,*) '%'
write(8, '("% % % % % % % ", " parameter ", "% % % % % % % %")')
write(8,'("% CFL number for stability    :", f5.2)') CFL
write(8,'("% Number of grid points       :", i5  )') n
write(8,'("% time                        :", f5.2)') t
write(8,'("% Difference scheme           :", A   )') scheme
write(8,*) '%'
write(8, '("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
write(8,*) 
write(8, '("%", 6x, "x", 8x, "rho", 8x, "u", 8x, "p")')
write(8, '(4f10.4)') (x(i), rho(i), u(i), p(i), i=1, n) 
close(8)

write(*,*),  "Done! The flow data are writen to 'solution.dat' "
end subroutine output
! --------------------------------------------------------------

! **************************************************************
subroutine processBar(t, Tend, process)
!  draw the program running process bar
implicit none
real(o), intent(in)   :: t            ! current time
real(o), intent(in)   :: Tend         ! final time
integer, intent(inout):: process

if (nsteps ==1 ) then
   write(*,*), 'Running process:'
   write(*,*), '-------20%-------40%-------60%-------80%------100%'
   write(*,'(A,$)') ' '
elseif (int(t/(Tend/50)) == 50) then
   write(*,'(A,$)') '>|'
   write(*,*) 
elseif (int(t/(Tend/50)) > process) then
   write(*,'(A,$)') '=' 
   process = int(t/(Tend/50))
end if
end subroutine processBar
! --------------------------------------------------------------

End program Tubeshock
