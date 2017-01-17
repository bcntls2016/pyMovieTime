program ReadWavefunction

Use Para_DerivnD

Implicit Real*8(A-H,O-Z)
Complex   (Kind=8), Allocatable :: psi0(:,:,:), psi02D(:,:), psi2D(:,:), Dxpsi2D(:,:), Dypsi2D(:,:)
complex(kind = 8), ALLOCATABLE	:: invar(:)			! Lambda (internal variables)
Real      (Kind=8), Allocatable :: den0(:,:,:), x0(:), y0(:), z0(:), den2D(:,:),x(:), y(:), Wx(:), Wy(:)
real (kind=8) :: mimpur = 132.9054d0
real (kind=8) :: planeoffset = 0
integer (kind=4) 	:: Ox, Oy, Oz
integer (kind=4)	:: iplaneoffset = 0
integer (kind = 4)	:: ninvar			! Number of components of lambda
integer (kind=4)	:: denmode = 42
integer (kind=4)	:: parammode = 42	! Density reading mode 
integer (kind=4) :: nthreads = 4
Integer   (Kind=4) :: nn(2), Icon=13
Character (len =1) :: cchar="#"
Logical (kind = 4) :: limp
logical (kind = 4)	:: tracking = .false.
Character (len=2)  :: plane = "xy"
Character (len=80) :: results = "density.res"
Character (len=80) :: den2dout = "den.dat", currents = "current.dat", inputdata = "djogger.dat"
Character (len=80) :: den1d0x = "den0-x.dat", den1d0y = "den0-y.dat", den1d0z = "den0-z.dat"
Character (len=80) :: den1dx = "den-1.dat", den1dy = "den-2.dat"
Data nx/436/, ny/436/, nZ/436/, hx/0.4d0/, hy/0.4d0/, hy/0.4d0/, &
		npd/13/, Km1/4/, ndmax/1/, nthreads/4/, npi/4/
Data fac/158.66624d0/,epsrho/1.d-6/
Data xi/-40.d0/, xf/40.d0/, yi/-40.d0/, yf/40.d0/, zi/-40.d0/, zf/40.d0/

Namelist/Input/nthreads,denmode,parammode,mimpur,plane,tracking,planeoffset, &
				nx,ny,nz,hx,hy,hz,npd,npi,Km1,icon,epsrho,xi,xf,yi,yf,zi,zf
read(5,nml=input)

K = npd    ! Km1 will be the number of derivatives for the Taylor expansion

Open(Unit=6, File=results)
Write(6,Input)

open (10, File = inputdata)
call titols(10,cchar,isalto)
select case (denmode)
	case (1)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,limp,ximp,yimp,zimp
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) den0
		psi0 = sqrt(den0)
	case (2)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,limp,ximp,yimp,zimp
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
	case (3)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,ximp,yimp,zimp,vximp,vyimp,vzimp
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
	case (4)
		read(10,*) xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0,ximp,yimp,zimp,vximp,vyimp,vzimp,ninvar
		ALLOCATE(invar(ninvar))
		Allocate(den0(nx0,ny0,nz0))
		Allocate(psi0(nx0,ny0,nz0))
		read(10,*) invar
		read(10,*) psi0
		den0 = Conjg(psi0) * psi0
	case default
		write(*,*)
		write(*,*) "You have chosen a 'denmode' unequal to {1,2,3,4}. Please modify '2Dden.settings' and choose one of:"
		write(*,*)
		write(*,*) "denmode = 1:	STATIC helium REAL density"
		write(*,*) "denmode = 2:	STATIC helium COMPLEX density supporting vorticity"
		write(*,*) "denmode = 3:	DYNAMIC helium density"
		write(*,*) "denmode = 4:	DYNAMIC helium density + electronic state of impurity"
		call EXIT(10)
end select
close(10)


Write(6,'("xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0...:",/,1p6E15.6,0p,3I5)') &
           xmax0,ymax0,zmax0,hx0,hy0,hz0,nx0,ny0,nz0

den0 = Abs(psi0)**2
dxyz0 = hx0*hy0*hz0
Ox = nx0/2+1
Oy = ny0/2+1
Oz = nz0/2+1

Write(6,'("den(0,0,0)....:",1p,E15.6)')den0(Ox,Oy,Oz)
Write(6,'("Norma(den)....:",1p,E15.6)')sum(den0)*dxyz0
Write(6,'("Norma(psi)....:",1p,E15.6)')sum((Abs(psi0)**2))*dxyz0

select case (plane)
	case ("xy")
		if (tracking) then
			iplaneoffset = nint((zimp+planeoffset)/hz0)
		else
			iplaneoffset = nint(planeoffset/hz0)
		end if
		Allocate(psi02D(nx0,ny0))
		Psi02D = Psi0(:,:,Oz+iplaneoffset)
	case ("xz")
		if (tracking) then
			iplaneoffset = nint((yimp+planeoffset)/hy0)
		else
			iplaneoffset = nint(planeoffset/hy0)
		end if
		ny0 = nz0
		hy0 = hz0
		ny = nz
		hy = hz
		yi = zi
		yf = zf
		Allocate(psi02D(nx0,ny0))
		Psi02D = Psi0(:,Oy+iplaneoffset,:)
	case ("yz")
		if (tracking) then
			iplaneoffset = nint((ximp+planeoffset)/hx0)
		else
			iplaneoffset = nint(planeoffset/hx0)
		end if		
		nx0 = ny0
		hx0 = hy0
		nx = ny
		hx = hy
		xi = yi
		xf = yf
		ny0 = nz0
		hy0 = hz0
		ny = nz
		hy = hz
		yi = zi
		yf = zf
		Allocate(psi02D(nx0,ny0))
		Psi02D = Psi0(Ox+iplaneoffset,:,:)
	case default
		write(*,*) "No correct plane selected, now exiting !"
		call exit
end select

Allocate(x(nx))
Allocate(y(ny))
Allocate(psi2D(nx,ny))
Allocate(den2D(nx,ny))
Allocate(Dxpsi2D(nx,ny))
Allocate(Dypsi2D(nx,ny))

nn(1) = nx; nn(2) = ny
xmax = (nx/2)*hx
ymax = (ny/2)*hy

Do ix = 1,nx
	x(ix) = -xmax+(ix-1)*hx
EndDo

Do iy = 1,ny
	y(iy) = -ymax+(iy-1)*hy
EndDo

Call INTERxy(nx,ny,x,y,Psi2D,nx0,ny0,hx0,hy0,Psi02D,K,KM1)
den2D = Abs(psi2D)**2

Open(Unit=1,file=den1dx)
Do ix = 1,nx
	Write(1,'(1p,2E15.6)')x(ix),den2D(ix,ny/2+1)
EndDo
Close (1)

Open(Unit=1,file=den1dy)
Do iy = 1,ny
	Write(1,'(1p,2E15.6)')y(iy),den2D(nx/2+1,iy)
EndDo
Close (1)

Call Init_deriv_p(npd,ndmax,nthreads)
Call DerivnD(1,nn,hx,1,psi2D,Dxpsi2D,Icon)
Call DerivnD(1,nn,hy,2,psi2D,Dypsi2D,Icon)

Pi=4.0d0*Datan(1.0d0)
TwoPi=2.0d0*Pi

ixi = (xi+xmax)/hx + 1.5
ixf = (xf+xmax)/hx + 1.5
iyi = (yi+ymax)/hy + 1.5
iyf = (yf+ymax)/hy + 1.5

write(6,'("ixi, ixf, iyi,iyf....:",4I5)')ixi,ixf,iyi,iyf
write(6,'("xi, xf, yi,yf....:",1p,4E15.6)')x(ixi),x(ixf),y(iyi),y(iyf)

!
! Calculo de la circulacion
!
nyi = iyf - iyi + 1
nxi = ixf - ixi + 1
Allocate(Wy(nyi))
Allocate(Wx(nxi))

If(nyi.Lt.npi)Then
	npi = 2*(nyi/2)
	Write(6,'("I had changed npi..:",I4)')npi
Endif  

If(nxi.Lt.npi)Then
	npi = 2*(nxi/2)      
	Write(6,'("I had changed npi..:",I4)')npi
Endif  

Call Simps(npi,2,hx,x,Wx,nxi)
Call Simps(npi,2,hy,y,Wy,nyi)

j = 0
xnorm1 = 0.d0
xnorm2 = 0.d0
Do ix = ixi,ixf
	j = j+1
	aux1 = Dimag(Dxpsi2D(ix,iyi)/Psi2D(ix,iyi))
	aux2 = Dimag(Dxpsi2D(ix,iyf)/Psi2D(ix,iyf))
	xnorm1 = xnorm1 + aux1*Wx(j)
	xnorm2 = xnorm2 + aux2*Wx(j)
EndDo
xnorm1 = xnorm1*hx
xnorm2 = xnorm2*hx

j=0
ynorm1 = 0.d0
ynorm2 = 0.d0
Do iy = iyi,iyf
	j = j+1
	aux1 = Dimag(Dypsi2D(ixi,iy)/Psi2D(ixi,iy))
	aux2 = Dimag(Dypsi2D(ixf,iy)/Psi2D(ixf,iy))
	ynorm1 = ynorm1 + Aux1*Wy(j)
	ynorm2 = ynorm2 + Aux2*Wy(j)
EndDo
ynorm1 = ynorm1*hy
ynorm2 = ynorm2*hy

Write(6,'("xnorm1,xnorm2...:",1p,2E15.6)')xnorm1,xnorm2
Write(6,'("ynorm1,ynorm2...:",1p,2E15.6)')ynorm1,ynorm2

xnorm = -xnorm1+xnorm2
ynorm = ynorm1-ynorm2
circ = (xnorm + ynorm)/Twopi
Write(6,'("Circulacion.....:",1p,E15.6)')circ

Open(Unit=8,File=currents)
vxmax = -1.d10
vxmin =  1.d10
vymax = -1.d10
vymin =  1.d10
Do ix = 1,nx
	Do iy = 1,ny
		If(den2D(ix,iy).Gt.epsrho) Then
			Stox = fac*dimag(Dxpsi2D(ix,iy)/psi2D(ix,iy))
			Stoy = fac*dimag(Dypsi2D(ix,iy)/psi2D(ix,iy))
			If(Stox.Gt.vxmax) vxmax = Stox
			If(Stoy.Gt.vymax) vymax = Stoy
			If(Stox.Lt.vxmin) vxmin = Stox
			If(Stoy.Lt.vymin) vymin = Stoy
			Write(8,'(1p,5E18.10)')x(ix), y(iy), den2D(ix,iy), Stox, Stoy
		Else
			Write(8,'(1p,5E18.10)')x(ix),y(iy),den2D(ix,iy),0.0d0,0.0d0
		Endif
	EndDo
	Write(8,*)
EndDo
Close(8)
Write(6,'("VxMin, VxMax....:",1p,2E15.6)')vxmin,vxmax
Write(6,'("VyMin, VyMax....:",1p,2E15.6)')vymin,vymax
Close(6)
End program ReadWavefunction