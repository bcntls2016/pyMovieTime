PROGRAM Parms

IMPLICIT NONE
character (len = 15)	:: a(4)				! Crap we don't need
character (len = 15)	:: b(2)				! Crap we don't need
character (len = 15)	:: c(4)				! Crap we don't need
character (len = 40)	:: inputdata = "djogger.dat"! Input file containing the wavefunction
character (len = 40)	:: params = "params.py"	! Output file containing the parameters
character (len = 1) 	:: cchar = "#"				! Character to be skipped on routine titols
character (len = 2)		:: plane = "xy"
logical (kind = 4)	:: limp = .false.	! Treat impurity CM coord quantum/classical
logical (kind = 4)	:: tracking = .false.
integer (kind = 4)	:: nx, ny, nz		! Number of points of the grid in each direction
integer (kind = 4)	:: isalto, npi			! Number of lines to skip at reading
integer (kind=4)	:: nthreads = 1		! Number of cpu cores to use 
integer (kind=4)	:: denmode = 42		! Density reading mode 
integer (kind=4)	:: parammode = 42	! Density reading mode
Integer (Kind=4)	:: Icon = 13, Km1 = 4, npd = 13
real (kind = 4)		:: time,rcm(3)		! current time, rCM of droplet 
real (kind = 8)		:: xmax, ymax, zmax	! Maximum values of x, y and z
real (kind = 8)		:: hx, hy, hz		! x, y and z step for the grid
real (kind = 8)		:: xi, xf, yi, yf, zi, zf		! x, y and z step for the grid
real (kind = 8)		:: rimp(3)			! Vector storing the position of the impurity
real (kind = 8)		:: vimp(3)			! Vector storing the velocity of the impurity
real (kind = 8)		:: mimpur			! Impurity mass (au)
real (kind = 8)		:: Ekin	= 0			! Kinetic Energy
real (kind = 8)		:: nparticles = 1000! Number of He particles
real (kind = 8)		:: epsrho			! Kinetic Energy
real (kind = 8)		:: planeoffset = 0
real (kind = 8) , parameter	:: mp_u = 0.020614837554503673d0  	! proton mass in Angs**-2 * Kelvin**-1, go figure!
real (kind = 8)	, parameter	::	Ktops=7.638235070684233d0		! Convert K to ps

Namelist/Input/inputdata,params,nthreads,denmode,parammode,mimpur,plane,tracking,planeoffset, &
				nx,ny,nz,hx,hy,hz,npd,npi,Km1,icon,epsrho,xi,xf,yi,yf,zi,zf
read(5,nml=input)

mimpur = mimpur * mp_u

open (unit=1, file=inputdata)
open (unit=2, file=params)
select case (parammode)
	case (1) ! Statics
		call titols(1,cchar,isalto)
		read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,rimp
		write(2,*) "# params.py"
		write(2,*)
		select case (plane)
			case ("xy")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", ymax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(2)
			case ("xz")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(3)
			case ("yz")				
				write(2,6000) "xmax = ", ymax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(2)
				write(2,6000) "yimp = ", rimp(3)
		end select
	case (2) ! Dynamics
		call titols(1,cchar,isalto)
		read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp
		Ekin = 0.5 * mimpur * sum(vimp * vimp)
		write(2,*) "# params.py"
		write(2,*)
		select case (plane)
			case ("xy")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", ymax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(2)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(2)/Ktops
			case ("xz")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops
			case ("yz")				
				write(2,6000) "xmax = ", ymax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(2)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6002) "vximp = ", 100*vimp(2)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops				
		end select
		write(2,6000) "ekin = ", Ekin
	case (3) ! Dynamics
		read(1,*)
		read(1,*)
		read(1,*) a,time
		read(1,*) c,rcm
		rewind(1)
		call titols(1,cchar,isalto)
		read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp
		Ekin = 0.5 * mimpur * sum(vimp * vimp)		
		!write(2,7106) xmax, zmax, rimp(1), rimp(3), 100*sqrt(sum(vimp*vimp))/Ktops, Ekin, time, rcm(3)
		write(2,*) "# params.py"
		write(2,*)
		select case (plane)
			case ("xy")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", ymax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(2)
				write(2,6000) "xcom = ", rcm(1)
				write(2,6000) "ycom = ", rcm(2)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(2)/Ktops				
			case ("xz")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6000) "xcom = ", rcm(1)
				write(2,6000) "ycom = ", rcm(3)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops				
			case ("yz")		
				write(2,6000) "xmax = ", ymax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(2)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6000) "xcom = ", rcm(2)
				write(2,6000) "ycom = ", rcm(3)
				write(2,6002) "vximp = ", 100*vimp(2)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops				
		end select
		write(2,6000) "ekin = ", Ekin
		write(2,6000) "time = ", time
	case (4) ! Dynamics
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*) a,time
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*) a,rcm
		read(1,*) b,nparticles
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*) b,Ekin
		rewind(1)
		call titols(1,cchar,isalto)
		read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp
		write(2,*) "# params.py"
		write(2,*)
		select case (plane)
			case ("xy")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", ymax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(2)
				write(2,6000) "xcom = ", rcm(1)
				write(2,6000) "ycom = ", rcm(2)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(2)/Ktops				
			case ("xz")				
				write(2,6000) "xmax = ", xmax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(1)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6000) "xcom = ", rcm(1)
				write(2,6000) "ycom = ", rcm(3)
				write(2,6002) "vximp = ", 100*vimp(1)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops				
			case ("yz")		
				write(2,6000) "xmax = ", ymax
				write(2,6000) "ymax = ", zmax
				write(2,6000) "ximp = ", rimp(2)
				write(2,6000) "yimp = ", rimp(3)
				write(2,6000) "xcom = ", rcm(2)
				write(2,6000) "ycom = ", rcm(3)
				write(2,6002) "vximp = ", 100*vimp(2)/Ktops
				write(2,6002) "vyimp = ", 100*vimp(3)/Ktops				
		end select
		write(2,6000) "ekin = ", Ekin
		write(2,6000) "time = ", time
		write(2,6001) "nparticles = ", nparticles
	case default
		write(*,*)
		write(*,*) "You have chosen a 'parammode' unequal to {1-4}. Please modify 'density.settings' and choose one of:"
		write(*,*)
		write(*,*) "parammode = 1:	STATIC helium density/wave-function."
		write(*,*) "parammode = 2:	DYNAMIC helium wave-function. Time and COM info added later."
		write(*,*) "parammode = 3:	DYNAMIC helium wave-function. Time and COM info from WF-file"
		write(*,*) "parammode = 4:	DYNAMIC helium wave-function. Time and COM info from WF-file, and more"
		call EXIT(10)
end select
close(1)
close(2)

6000 format(A7, 1X, F7.2)
6001 format(A13, 1X, F7.2)
6002 format(A8, 1X, F7.2)
7100 format(2F15.5,1x,20E12.3)
7105 format(F5.2, 2x, F5.2, 2x, F6.2, 2x, F6.2, 2x, F6.1, 2x, F6.1)
7106 format(F5.2, 2x, F5.2, 2x, F6.2, 2x, F6.2, 2x, F6.1, 2x, F6.1, 2x, F5.1, 2x, F8.4)
7110 format()

END PROGRAM


!----------------------------------------------------------------------
!--                 Subroutine TITOLS                               ---
!----------------------------------------------------------------------
!
!   Esta rutina posiciona el puntero de lectura de un fichero
!   hasta la primera linea cuya primera columna NO comience por el
!   caracter  'cchar' y devuelve en la variable isalto el numero de
!   lineas que se ha saltado.
!
!   La variable ulog indica el numero de unidad logica a utilizar.
!
subroutine titols(ulog,cchar,isalto)
implicit none
character (LEN=1) :: pchar, pcolumn, cchar
integer :: nl,i,ulog,isalto

nl = 0
pcolumn = cchar
do while(pcolumn.eq.cchar)
	read(ulog,5000,end = 9000) pchar
	pcolumn = pchar
	nl = nl + 1
end do
rewind(ulog)
isalto = nl - 1
do i = 1,isalto
	read(ulog,*)
end do
return
9000 print *,' Ey Mister this file is too short ...'  
     error stop 'STOP due to severe errors in routine TITOLS'
5000 format(A1)
end
