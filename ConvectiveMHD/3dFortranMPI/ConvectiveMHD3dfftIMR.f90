PROGRAM main	
	!-----------------------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program numerically solves the 3D convective incompressible magnetohydrodynamic 
	! equations on a Cubic Domain [0,2pi]x[0,2pi]x[0,2pi] using fourier pseudo-spectral 
	! methods and Implicit Midpoint rule timestepping. 
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  Nt				= number of timesteps to take
	!  Tmax				= maximum simulation time
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD     = value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Re				= Reynolds number
  !  Rem      = Magnetic reynolds number
  !  Ra       = Rayleigh number
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  count			= keep track of information written to disk
	!  iol				= size of array to write to disk
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfxyz			= Forward 3d fft plan 
	!  planbxyz			= Backward 3d fft plan
	!  dt				    = timestep
	! .. Arrays ..
	!  u				  = velocity in x direction
	!  v				  = velocity in y direction
	!  w				  = velocity in z direction
	!  uold				= velocity in x direction at previous timestep
	!  vold				= velocity in y direction at previous timestep
	!  wold				= velocity in z direction at previous timestep
	!  ux				  = x derivative of velocity in x direction
	!  uy				  = y derivative of velocity in x direction
	!  uz				  = z derivative of velocity in x direction
	!  vx				  = x derivative of velocity in y direction
	!  vy				  = y derivative of velocity in y direction
	!  vz				  = z derivative of velocity in y direction
	!  wx				  = x derivative of velocity in z direction
	!  wy				  = y derivative of velocity in z direction
	!  wz				  = z derivative of velocity in z direction
	!  uxold			= x derivative of velocity in x direction
	!  uyold			= y derivative of velocity in x direction
	!  uzold			= z derivative of velocity in x direction
	!  vxold			= x derivative of velocity in y direction
	!  vyold			= y derivative of velocity in y direction
	!  vzold			= z derivative of velocity in y direction
	!  wxold			= x derivative of velocity in z direction
	!  wyold			= y derivative of velocity in z direction
	!  wzold			= z derivative of velocity in z direction
	!  utemp			= temporary storage of u to check convergence
	!  vtemp			= temporary storage of v to check convergence
	!  wtemp			= temporary storage of w to check convergence
	!  temp_r			= temporary storage for untransformed variables
	!  uhat				= Fourier transform of u
	!  vhat				= Fourier transform of v
	!  what				= Fourier transform of w
	!  rhsuhatfix	= Fourier transform of righthand side for u for timestepping
	!  rhsvhatfix	= Fourier transform of righthand side for v for timestepping
	!  rhswhatfix	= Fourier transform of righthand side for w for timestepping
	!  bx				  = x component of magnetic field
	!  by				  = y component of magnetic field
	!  bz				  = z component of magnetic field
	!  bxold			= x component of magnetic field at previous timestep
	!  byold			= y component of magnetic field at previous timestep
	!  bzold			= z component of magnetic field at previous timestep
	!  bxx				= x derivative of x component of magnetic field
	!  bxy				= y derivative of x component of magnetic field
	!  bxz				= z derivative of x component of magnetic field
	!  byx				= x derivative of y component of magnetic field
	!  byy				= y derivative of y component of magnetic field
	!  byz				= z derivative of y component of magnetic field
	!  bzx				= x derivative of z component of magnetic field
	!  bzy				= y derivative of z component of magnetic field
	!  bzz				= z derivative of z component of magnetic field
	!  bxxold			= x derivative of x component of magnetic field at previous timestep
	!  bxyold			= y derivative of x component of magnetic field at previous timestep
	!  bxzold			= z derivative of x component of magnetic field at previous timestep
	!  byxold			= x derivative of y component of magnetic field at previous timestep
	!  byyold			= y derivative of y component of magnetic field at previous timestep
	!  byzold			= z derivative of y component of magnetic field at previous timestep
	!  bzxold			= x derivative of z component of magnetic field at previous timestep
	!  bzyold			= y derivative of z component of magnetic field at previous timestep
	!  bzzold			= z derivative of z component of magnetic field at previous timestep
	!  bxtemp			= temporary storage of bx to check convergence
	!  bytemp			= temporary storage of by to check convergence
	!  bztemp			= temporary storage of bz to check convergence
	!  bxhat			= Fourier transform of bx
	!  byhat			= Fourier transform of by
	!  bzhat			= Fourier transform of bz
  !  theta      = Temperature
  !  thetax     = x derivative of temperature
  !  thetay     = y derivative of temperature
  !  thetaz     = z derivative of temperature
  !  thetaxold  = x derivative of temperature  at previous timestep
  !  thetayold  = y derivative of temperature  at previous timestep
  !  thetazold  = z derivative of temperature  at previous timestep
  !  thetatemp  = temporary storage of temperature to check convergence
  !  thetahat   = Fourier transform of temperature 
  !  rhsthetafix    = Fourier transform of righthand side for temperature for timestepping
	!  rhsbxhatfix		= Fourier transform of righthand side for bx for timestepping
	!  rhsbyhatfix		= Fourier transform of righthand side for by for timestepping
	!  rhsbzhatfix		= Fourier transform of righthand side for bz for timestepping	
	!  nonlinuhat			= Fourier transform of nonlinear term for u
	!  nonlinvhat  		= Fourier transform of nonlinear term for v
	!  nonlinwhat			= Fourier transform of nonlinear term for w
  !  nonlinthetahat = Fourier transform of nonlinear term for theta
	!  phat				    = Fourier transform of nonlinear term for pressure, p
	!  temp_c			    = temporary storage for Fourier transforms
	!  realtemp			  = Real storage
	!
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	!  x				= x locations
	!  y				= y locations
	!  z				= y locations
	!  time				= times at which save data
	!  name_config			= array to store filename for data to be saved
	!    		
	! REFERENCES
	! Draws on work on "Infinte prandtl convection on the 2D torus" by Doering, Muite and Whitehead (forthcoming) 
	!
	! ACKNOWLEDGEMENTS
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	!
	! This program has not been optimized to use the least amount of memory
	! but is intended as an example only for which all states can be saved
	!
	!--------------------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! 2DECOMP&FFT -- Fast Fourier Transform in the West Library
	!			(http://2decomp.org/)
				 
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	USE MPI
	IMPLICIT NONE	
	! declare variables
   	INTEGER(kind=4), PARAMETER 		:: Nx=64
	INTEGER(kind=4), PARAMETER 		:: Ny=64
	INTEGER(kind=4), PARAMETER 		:: Nz=64
   	INTEGER(kind=4), PARAMETER 		:: Lx=1
	INTEGER(kind=4), PARAMETER 		:: Ly=1
	INTEGER(kind=4), PARAMETER 		:: Lz=1
	INTEGER(kind=4), PARAMETER		:: Nt=5
	REAL(kind=8), PARAMETER			:: dt=0.05d0/Nt
	REAL(kind=8), PARAMETER			:: Re=1.0d0	
	REAL(kind=8), PARAMETER			:: Rem=1.0d0	
	REAL(kind=8), PARAMETER			:: Ra=1.0d0	
	REAL(kind=8), PARAMETER			:: tol=0.1d0**10
	REAL(kind=8), PARAMETER			:: theta=0.0d0

	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::	ReInv=1.0d0/REAL(Re,kind(0d0))
	REAL(kind=8), PARAMETER		::	RemInv=1.0d0/REAL(Rem,kind(0d0))
	REAL(kind=8), PARAMETER		::  	dtInv=1.0d0/REAL(dt,kind(0d0)) 
	REAL(kind=8)			:: 	scalemodes,chg,factor
	REAL(kind=8), DIMENSION(:), ALLOCATABLE		:: x, y, z, time,mychg,allchg
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: u, v, w,&
							ux, uy, uz,&
							vx, vy, vz,&
							wx, wy, wz,&
							uold, uxold, uyold, uzold,&
							vold, vxold, vyold, vzold,&
							wold, wxold, wyold, wzold,&
							utemp, vtemp, wtemp, temp_r, &
							bx, by, bz,&
							bxx, bxy, bxz,&
							byx, byy, byz,&
							bzx, bzy, bzz,&
							bxold, bxxold, bxyold, bxzold,&
							byold, byxold, byyold, byzold,&
							bzold, bzxold, bzyold, bzzold,&
 							bxtemp, bytemp, bztemp, &
              theta, thetaold, thetatemp, &
              thetax, thetay, thetaz, &
              thetaxold, thetayold, thetazold
																	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	:: kx, ky, kz						
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: uhat, vhat, what, &
							bxhat, byhat, bzhat, &
							rhsuhatfix, rhsvhatfix, &
							rhswhatfix, rhsbxhatfix, &
							rhsbyhatfix, rhsbzhatfix, &
							nonlinuhat, nonlinvhat, &
							nonlinwhat, phat, temp_c, &
              thetahat, nonlinthetahat, rhsthetahatfix
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE 	::  realtemp
	! MPI and 2DECOMP variables
	TYPE(DECOMP_INFO)				::  decomp
	INTEGER(kind=MPI_OFFSET_KIND) 			::  filesize, disp
	INTEGER(kind=4)					::  p_row=0, p_col=0, numprocs, myid, ierr	
	
	! variables used for saving data and timing
	INTEGER(kind=4)					:: count, iol 
	INTEGER(kind=4)					:: i,j,k,n,t,allocatestatus
	INTEGER(kind=4)					:: ind, numberfile
	CHARACTER*100			 		:: name_config
	INTEGER(kind=4)					:: start, finish, count_rate 
	
	! initialisation of 2DECOMP&FFT
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	! do automatic domain decomposition
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	! get information about domain decomposition choosen
	CALL decomp_info_init(Nx,Ny,Nz,decomp)
	! initialise FFT library
	CALL decomp_2d_fft_init
	IF (myid.eq.0) THEN
	    PRINT *,'Grid:',Nx,'X',Ny,'Y',Nz,'Z'
		PRINT *,'dt:',dt
	END IF	
	ALLOCATE(x(1:Nx),y(1:Ny),z(1:Nz),time(1:Nt+1),mychg(1:6),allchg(1:6),&
				u(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),& 
 				v(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				w(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				ux(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				uzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				vzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				wzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				utemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				vtemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				wtemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				temp_r(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),& 
 				by(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				bz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				bxx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzx(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzy(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				theta(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetax(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetay(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetaz(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				byzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzyold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bzzold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				bxtemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				bytemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
 				bztemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetaold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetaxold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetayold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetazold(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				thetatemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)),&
				kx(1:Nx),ky(1:Ny),kz(1:Nz),&
				uhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
				vhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
	 			what(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
				bxhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
				byhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
	 			bzhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
	 			rhsuhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				rhsvhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				rhswhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
	 			rhsbxhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				rhsbyhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				rhsbzhatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				rhsthetahatfix(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				nonlinuhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				nonlinvhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				nonlinwhat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				nonlinthetahat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				phat(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
 				temp_c(decomp%zst(1):decomp%zen(1),&
   					decomp%zst(2):decomp%zen(2),&
   					decomp%zst(3):decomp%zen(3)),&
				realtemp(decomp%xst(1):decomp%xen(1),&
   					decomp%xst(2):decomp%xen(2),&
   					decomp%xst(3):decomp%xen(3)), stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP
	IF (myid.eq.0) THEN
	 	PRINT *,'allocated space'
	END IF	

	! setup fourier frequencies in x-direction
	DO i=1,Nx/2+1
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO	
	ind=1
	DO i=-Nx/2,Nx/2-1
		x(ind)=2.0d0*pi*REAL(i,kind(0d0))*Lx/REAL(Nx,kind(0d0))
		ind=ind+1
	END DO
	! setup fourier frequencies in y-direction
	DO j=1,Ny/2+1
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	ky(1+Ny/2)=0.0d0
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO	
	ind=1
	DO j=-Ny/2,Ny/2-1
		y(ind)=2.0d0*pi*REAL(j,kind(0d0))*Ly/REAL(Ny,kind(0d0))
		ind=ind+1
	END DO
	! setup fourier frequencies in z-direction
	DO k=1,Nz/2+1
		kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz  			
	END DO
	kz(1+Nz/2)=0.0d0
	DO k = 1,Nz/2 -1
		kz(k+1+Nz/2)=-kz(1-k+Nz/2)
	END DO	
	ind=1
	DO k=-Nz/2,Nz/2-1
		z(ind)=2.0d0*pi*REAL(k,kind(0d0))*Lz/REAL(Nz,kind(0d0))
		ind=ind+1
	END DO
	scalemodes=1.0d0/REAL(Nx*Ny*Nz,kind(0d0))
	IF (myid.eq.0) THEN
		PRINT *,'Setup grid and fourier frequencies'
	END IF	

	! Initial conditions for fluid flow field
	!initial conditions for Taylor-Green vortex
!	factor=2.0d0/sqrt(3.0d0)
!	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
!		u(i,j,k)=factor*sin(theta+2.0d0*pi/3.0d0)*sin(x(i))*cos(y(j))*cos(z(k))
!	END DO; END DO; END DO
!	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
!		v(i,j,k)=factor*sin(theta-2.0d0*pi/3.0d0)*cos(x(i))*sin(y(j))*cos(z(k))
!	END DO ; END DO ; END DO
!	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
!		w(i,j,k)=factor*sin(theta)*cos(x(i))*cos(y(j))*sin(z(k))
!	END DO ; END DO ; END DO

	time(1)=0.0d0
	factor=sqrt(3.0d0)
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
						+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
	END DO; END DO; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		v(i,j,k)=0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
						-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
	END DO ; END DO ; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
	END DO ; END DO ; END DO

	CALL decomp_2d_fft_3d(u,uhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(v,vhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(w,what,DECOMP_2D_FFT_FORWARD)
	
	! derivative of u with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD)	

	! derivative of v with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD)		

	! derivative of w with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD)		
	! save initial data
	n=0
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegax'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
	!omegay
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegay'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
	!omegaz
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegaz'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
  

	! Initial conditions for magnetic field
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		bx(i,j,k)=-0.5*( sin(y(j))*sin(z(k)))
	END DO; END DO; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		by(i,j,k)=0.5*(  sin(x(i))*sin(z(k)))
	END DO ; END DO ; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		bz(i,j,k)=cos(x(i))*cos(y(j))
	END DO ; END DO ; END DO

	CALL decomp_2d_fft_3d(bx,bxhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(by,byhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(bz,bzhat,DECOMP_2D_FFT_FORWARD)
	
	! derivative of u with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bxhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,bxx,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bxhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,bxy,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bxhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,bxz,DECOMP_2D_FFT_BACKWARD)	

	! derivative of v with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=byhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,byx,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=byhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,byy,DECOMP_2D_FFT_BACKWARD)	
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=byhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,byz,DECOMP_2D_FFT_BACKWARD)		

	! derivative of w with respect to x, y, and z 
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bzhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bzhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)		
	DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
		temp_c(i,j,k)=bzhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,bzz,DECOMP_2D_FFT_BACKWARD)		
	! save initial data
	n=0
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(bzy(i,j,k)-byz(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegabx'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
	!omegay
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(bxz(i,j,k)-bzx(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegaby'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
	!omegaz
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		realtemp(i,j,k)=REAL(byx(i,j,k)-bxy(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegabz'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)

      
        !start timer
        CALL system_clock(start,count_rate)
	DO n=1,Nt
		!fixed point iteration
		! store old fluid field terms
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			uold(i,j,k)=u(i,j,k)
			uxold(i,j,k)=ux(i,j,k)
			uyold(i,j,k)=uy(i,j,k)
			uzold(i,j,k)=uz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			vold(i,j,k)=v(i,j,k)
			vxold(i,j,k)=vx(i,j,k)
			vyold(i,j,k)=vy(i,j,k)
			vzold(i,j,k)=vz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			wold(i,j,k)=w(i,j,k)
			wxold(i,j,k)=wx(i,j,k)
			wyold(i,j,k)=wy(i,j,k)
			wzold(i,j,k)=wz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
                        rhsuhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*uhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
			rhsvhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*vhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
			rhswhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*what(i,j,k) 
		END DO ; END DO ; END DO
		! Store old Magnetic field terms
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			bxold(i,j,k)=bx(i,j,k)
			bxxold(i,j,k)=bxx(i,j,k)
			bxyold(i,j,k)=bxy(i,j,k)
			bxzold(i,j,k)=bxz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			byold(i,j,k)=by(i,j,k)
			byxold(i,j,k)=byx(i,j,k)
			byyold(i,j,k)=byy(i,j,k)
			byzold(i,j,k)=byz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			bzold(i,j,k)=bz(i,j,k)
			bzxold(i,j,k)=bzx(i,j,k)
			bzyold(i,j,k)=bzy(i,j,k)
			bzzold(i,j,k)=bzz(i,j,k)
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
                        rhsbxhatfix(i,j,k) = (dtInv+(0.5d0*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*bxhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
			rhsbyhatfix(i,j,k) = (dtInv+(0.5*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*byhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
			rhsbzhatfix(i,j,k) = (dtInv+(0.5*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*bzhat(i,j,k) 
		END DO ; END DO ; END DO

		
		chg=1
		DO WHILE (chg .gt. tol)
			! Fluid field
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(ux(i,j,k)+uxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(uy(i,j,k)+uyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(uz(i,j,k)+uzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(bxx(i,j,k)+bxxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(bxy(i,j,k)+bxyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(bxz(i,j,k)+bxzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinuhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(vx(i,j,k)+vxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(vy(i,j,k)+vyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(vz(i,j,k)+vzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(byx(i,j,k)+byxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(byy(i,j,k)+byyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(byz(i,j,k)+byzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinvhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(wx(i,j,k)+wxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(wy(i,j,k)+wyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(wz(i,j,k)+wzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(bzx(i,j,k)+bzxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(bzy(i,j,k)+bzyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(bzz(i,j,k)+bzzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinwhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				phat(i,j,k)=-1.0d0*( kx(i)*nonlinuhat(i,j,k)&
							+ky(j)*nonlinvhat(i,j,k)&
							+kz(k)*nonlinwhat(i,j,k))&
							/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+0.1d0**13)
			END DO ; END DO ; END DO

			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				uhat(i,j,k)=(rhsuhatfix(i,j,k)-nonlinuhat(i,j,k)-kx(i)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				vhat(i,j,k)=(rhsvhatfix(i,j,k)-nonlinvhat(i,j,k)-ky(j)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				what(i,j,k)=(rhswhatfix(i,j,k)-nonlinwhat(i,j,k)-kz(k)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO

			! derivative of u with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD)	

			! derivative of v with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD)	

			! derivative of w with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				utemp(i,j,k)=u(i,j,k)
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				vtemp(i,j,k)=v(i,j,k)
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				wtemp(i,j,k)=w(i,j,k)
			END DO ; END DO ; END DO

			CALL decomp_2d_fft_3d(uhat,u,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(vhat,v,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(what,w,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				u(i,j,k)=u(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				v(i,j,k)=v(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				w(i,j,k)=w(i,j,k)*scalemodes
			END DO ; END DO ; END DO

			! Magnetic field
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(bxx(i,j,k)+bxxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(bxy(i,j,k)+bxyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(bxz(i,j,k)+bxzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(ux(i,j,k)+uxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(uy(i,j,k)+uyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(uz(i,j,k)+uzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinuhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(byx(i,j,k)+byxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(byy(i,j,k)+byyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(byz(i,j,k)+byzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(vx(i,j,k)+vxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(vy(i,j,k)+vyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(vz(i,j,k)+vzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinvhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(bzx(i,j,k)+bzxold(i,j,k))&
							+(v(i,j,k)+vold(i,j,k))*(bzy(i,j,k)+bzyold(i,j,k))&
							+(w(i,j,k)+wold(i,j,k))*(bzz(i,j,k)+bzzold(i,j,k))&
							-(bx(i,j,k)+bxold(i,j,k))*(wx(i,j,k)+wxold(i,j,k))&
							-(by(i,j,k)+byold(i,j,k))*(wy(i,j,k)+wyold(i,j,k))&
							-(bz(i,j,k)+bzold(i,j,k))*(wz(i,j,k)+wzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_r,nonlinwhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				phat(i,j,k)=-1.0d0*( kx(i)*nonlinuhat(i,j,k)&
							+ky(j)*nonlinvhat(i,j,k)&
							+kz(k)*nonlinwhat(i,j,k))&
							/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+0.1d0**13)
			END DO ; END DO ; END DO

			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				bxhat(i,j,k)=(rhsbxhatfix(i,j,k)-nonlinuhat(i,j,k)-kx(i)*phat(i,j,k))/&
							(dtInv-(0.5d0*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				byhat(i,j,k)=(rhsbyhatfix(i,j,k)-nonlinvhat(i,j,k)-ky(j)*phat(i,j,k))/&
							(dtInv-(0.5d0*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				bzhat(i,j,k)=(rhsbzhatfix(i,j,k)-nonlinwhat(i,j,k)-kz(k)*phat(i,j,k))/&
							(dtInv-(0.5d0*RemInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO

			! derivative of u with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=bxhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bxx,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=bxhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bxy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bxz,DECOMP_2D_FFT_BACKWARD)	

			! derivative of v with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=byhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,byx,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=byhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,byy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=byhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,byz,DECOMP_2D_FFT_BACKWARD)	

			! derivative of w with respect to x, y, and z 
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=bzhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bzx,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=bzhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bzy,DECOMP_2D_FFT_BACKWARD)	
			DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
				temp_c(i,j,k)=bzhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL decomp_2d_fft_3d(temp_c,bzz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				bxtemp(i,j,k)=bx(i,j,k)
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				bytemp(i,j,k)=by(i,j,k)
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				bztemp(i,j,k)=bz(i,j,k)
			END DO ; END DO ; END DO

			CALL decomp_2d_fft_3d(bxhat,bx,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(byhat,by,DECOMP_2D_FFT_BACKWARD)	
			CALL decomp_2d_fft_3d(bzhat,bz,DECOMP_2D_FFT_BACKWARD)	

			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				bx(i,j,k)=bx(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				by(i,j,k)=by(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				bz(i,j,k)=bz(i,j,k)*scalemodes
			END DO ; END DO ; END DO
						
			mychg(1) =maxval(abs(utemp-u))
			mychg(2) =maxval(abs(vtemp-v))
			mychg(3) =maxval(abs(wtemp-w))
			mychg(4) =maxval(abs(bxtemp-bx))
			mychg(5) =maxval(abs(bytemp-by))
			mychg(6) =maxval(abs(bztemp-bz))
			CALL MPI_ALLREDUCE(mychg,allchg,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
			chg=allchg(1)+allchg(2)+allchg(3)+allchg(4)+allchg(5)+allchg(6)
			IF (myid.eq.0) THEN
				PRINT *,'chg:',chg
			END IF
		END DO
		time(n+1)=n*dt

                IF (myid.eq.0) THEN	
			PRINT *,'time',n*dt
		END IF 
		
                !save omegax, omegay, omegaz, omegabx, omegaby, and omegabz	
		!omegax
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegax'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
		!omegay
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegay'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
		!omegaz
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegaz'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
		!omegabx
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(bzy(i,j,k)-byz(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegabx'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
		!omegaby
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(bxz(i,j,k)-bzx(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegaby'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
		!omegabz
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			realtemp(i,j,k)=REAL(byx(i,j,k)-bxy(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegabz'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp,decomp)
                
	END DO
        
        CALL system_clock(finish,count_rate)

        IF (myid.eq.0) then
           PRINT *, 'Program took', REAL(finish-start)/REAL(count_rate), 'for main timestepping loop'
        END IF

	IF (myid.eq.0) THEN
		name_config = './data/tdata.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO n=1,1+Nt
			WRITE(11,*) time(n)
		END DO
		CLOSE(11)

		name_config = './data/xcoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO i=1,Nx
			WRITE(11,*) x(i)
		END DO
		CLOSE(11)	

		name_config = './data/ycoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO j=1,Ny
			WRITE(11,*) y(j)
		END DO
		CLOSE(11)
		
		name_config = './data/zcoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO k=1,Nz
			WRITE(11,*) z(k)
		END DO
		CLOSE(11)
		PRINT *,'Saved data'
	END IF
	
	
        ! clean up 
  	CALL decomp_2d_fft_finalize
  	CALL decomp_2d_finalize

	DEALLOCATE(x,y,z,time,mychg,allchg,&
		u,v,w,ux,uy,uz,vx,vy,vz,&
		wx,wy,wz,uold,uxold,&
		uyold,uzold,vold,vxold,&
		vyold,vzold,wold,wxold,&
		wyold,wzold,utemp,vtemp,&
 		wtemp,temp_r,&
		bx,by,bz,bxx,bxy,bxz,&
		byx,byy,byz,bzx,bzy,bzz,&
		bxold,bxxold,bxyold,bxzold,&
		byold,byxold,byyold,byzold,&
		bzold,bzxold,bzyold,bzzold,&
		bxtemp,bytemp,bztemp,kx,ky,kz,&
		uhat,vhat,what,bxhat,&
		byhat,bzhat,rhsuhatfix,rhsvhatfix,&
 		rhswhatfix,rhsbxhatfix,rhsbyhatfix,rhsbzhatfix,&
 		nonlinuhat,nonlinvhat,nonlinwhat,phat,&
 		temp_c,realtemp,stat=AllocateStatus)		
	IF (AllocateStatus .ne. 0) STOP
	IF (myid.eq.0) THEN
		PRINT *,'Program execution complete'
	END IF
	CALL MPI_FINALIZE(ierr)		

	END PROGRAM main
