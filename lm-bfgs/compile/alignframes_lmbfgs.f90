!-------------------------------------------------------
! alignframes_lbfgsb                       Dated 10/2/15
!-------------------------------------------------------
! JLR 5/14
! JLR 10/14 Allow specifying zeroframe to not be first frame
!      2/15 Allow paths to be specified; allow arbitrary length extension change
!-------------------------------------------------------

PROGRAM Alignframes_lbfgsb
implicit none
integer, parameter                :: dp=kind(0.d0) 
real(dp), parameter               :: pi=3.141592654
! Stack parameters
integer                           :: nxyz(3),mxyz(3),mode,nx,ny
real(dp)                          :: dmin,dmax,dmean!,cell(6)
! I/O
character(150)                    :: inlist,inmovie,tempname,vectorfile
character(150)                    :: inpath,outpath
character(20)                     :: vecext
integer                           :: io
! Counters
integer                           :: i,filenumb
! Coordinate infor
integer                           :: boxsize,nframes,framefirst,framelast,zeroframe
! Images
real,allocatable                  :: inimage(:,:)
real(dp), allocatable             :: holder(:,:)
real(dp), allocatable             :: weightmask(:,:)
double complex, allocatable       :: currentpartsc(:,:,:)
! Particle movement
real(dp)                          :: rmax1,rmax2,psize
real(dp), allocatable             :: xshift(:),yshift(:)
! Statistics
real(dp)                          :: mean,stddev,nsigma
! FFTW
real(dp), allocatable             :: image_fftw(:,:)
double complex, allocatable       :: imagec_fftw(:,:)
real(dp)                          :: plan2df, plan2dr
! optimization
real(dp)                          :: rmax1sq,rmax2sq
real(dp)                          :: boxsizesqinv,nframesinv,radtodeg
integer                           :: boxsizeby2plus1
! temp factor
real(dp)                          :: bfactor
! minimizer
real(dp)                          :: corrvsaverage,dcorrvsaveragedx,dcorrvsaveragedy
real(dp),allocatable              :: cx(:),cxtemp(:),cy(:),cytemp(:),coeff(:),coefftemp(:)
! l-bfgsb-b
integer                           :: nvar ! number of variables (n)
integer,  parameter               :: m = 15, iprint = 1!100
real(dp), parameter               :: pgtol  = 0
real(dp)                          :: factr ! factr=1.0d7 is moderate, 1.0d1 is extreme, 0 is not use
character(len=60)                 :: task, csave
logical                           :: lsave(4)
integer                           :: isave(44)
real(dp)                          :: funval
real(dp)                          :: dsave(29)
integer,  allocatable             :: nbd(:), iwa(:)
real(dp), allocatable             :: lowbound(:), upbound(:), grad(:), wa(:)
! temp
!real(dp)                          :: testderx,testdery
!real(dp), parameter               :: epsilon=0.000001

include '/usr/include/fftw3.f'

! TEMP

read (5,*) inlist
read (5,*) boxsize,psize,nsigma,rmax1,rmax2
read (5,*) bfactor
read (5,*) framefirst,framelast,zeroframe
read (5,*) factr
read (5,'(a)') inpath
read (5,'(a)') outpath
read (5,*) vecext

! File stream
! Stream 1: Input movie
! Stream 2: 
! Stream 3: Input list
! Stream 4: None
! Stream 7: Output vector file (shift file for whole frame) 
! Stream 8: 

open(unit=3,file=inlist,status="old")
boxsizesqinv = 1/dble(boxsize**2)
boxsizeby2plus1 = boxsize/2+1
! Restate resolution limits in terms of pixel radii
rmax1=(psize*boxsize)/rmax1
if (rmax1<1) rmax1=1.01 ! Never use the origin of the FT for alignment
rmax1sq=rmax1**2
rmax2=(psize*boxsize)/rmax2
rmax2sq=rmax2**2
radtodeg=180_dp/pi


! Allocate FT arrays and plan
! FFTW arrays
allocate (image_fftw(boxsize,boxsize))                           ; image_fftw=0
allocate (imagec_fftw(boxsizeby2plus1,boxsize))                  ; imagec_fftw=0
allocate (weightmask(boxsizeby2plus1,boxsize))                   ; weightmask=0
allocate (holder(boxsize,boxsize))                               ; holder=0
write (*,*) "Generating weightmaskc"
call temperaturemask(weightmask,boxsize,psize,bfactor)
! Plan FFTs
call Dfftw_plan_dft_r2c_2d(plan2df,boxsize,boxsize,image_fftw,imagec_fftw,fftw_measure)
call Dfftw_plan_dft_c2r_2d(plan2dr,boxsize,boxsize,imagec_fftw,image_fftw,fftw_measure)
write (*,*) "FTs Planned"

filenumb=0
Eachmovie: do
  ! Read a filename from list
  read (3,*,iostat=io) tempname
  ! Determine input image name
  Fileexists: if (io>0) then
    write (*,*) "problem with input list format"
    exit
  else if (io<0) THEN
    write (*,*) "End of file list reached"
    exit
  else
    call Changepath(tempname,inmovie,inpath)
    call Changepath(tempname,inmovie,inpath)
    write (*,*) inmovie
    filenumb=filenumb+1
    ! Open input movie
    call Imopen(1,inmovie,"RO")
    call Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
    nx         = nxyz(1)
    ny         = nxyz(2)
    if (framelast.eq.0) framelast=nxyz(3)
    nframes=framelast-framefirst+1
    nframesinv=1/dble(nframes)
    nxyz=(/boxsize,boxsize,100/)
    mxyz=nxyz
    dmin=1E6
    dmax=-1E6
    dmean=0
    ! Null motion model, where each frame is independent
    nvar=nframes*2
    ! Open output stack of averaged particles
    call Replaceextension_arblen(inmovie,tempname,vecext)
    call Changepath(tempname,vectorfile,outpath)
    ! Open vector file
    ! Open input coordinate file and output vectorfield file
    open(7,file=vectorfile,status='replace')
    ! Allocate arrays if necessary
    if (filenumb.eq.1) then
      write (*,*) "First time through loop. Allocating arrays"
      ! Allocate image arrays
      allocate (inimage(nx,ny))                                        ; inimage=0
      ! Complex Arrays
      allocate (currentpartsc(boxsizeby2plus1,boxsize,nframes))        ; currentpartsc=0
      ! Shift info
      allocate (xshift(nframes))                                       ; xshift=0
      allocate (yshift(nframes))                                       ; yshift=0
      ! minimizer
      allocate (lowbound(nvar))                                        ; lowbound=0
      allocate (upbound(nvar))                                         ; upbound=0
      allocate (nbd(nvar))                                             ; nbd=0
      allocate (grad(nvar))                                            ; grad=0
      allocate (iwa(3*nvar))
      allocate (wa(2*m*nvar+5*nvar+11*m*m+8*m) )
      ! Coefficients
      allocate (cx(nvar/2))                                            ; cx=0 
      allocate (cy(nvar/2))                                            ; cy=0 
      allocate (cxtemp(nvar/2))                                        ; cxtemp=0
      allocate (cytemp(nvar/2))                                        ; cytemp=0
      allocate (coeff(nvar))                                           ; coeff=0
      allocate (coefftemp(nvar))                                       ; coefftemp=0
    end if
    ! Prepare output stack header
    ! Set minimizer conditions
    coeff=0 ! Initial guess for minimizer
    nbd=0        ! 2=All variables have a lower and upper boundary
    !(if nbd=0, no bounds, 1 lower bound only, 2 upper and lower, 3 upper only
    Eachframe: do i=1,nframes
      write (*,*) "Reading frame",i+framefirst-1
      inimage=0
      call Imposn(1,i+framefirst-2,0) 
      call Irdpas(1,inimage(:,:),nx,ny,0,nx-1,0,ny-1)
      mean=(sum(inimage))*real(boxsizesqinv)
      inimage=inimage-real(mean)
      stddev=sqrt(sum(dble(inimage)**2)*boxsizesqinv)
      where (abs(dble(inimage)).gt.nsigma*stddev)
        inimage=0
      end where
      ! Put the particle frame in the real and complex stacks
      call Cram(dble(inimage(:,:)),holder,nx,ny,boxsize,boxsize)
      call Floatpart(holder,image_fftw,boxsize,boxsize) ! Float each frame
      call Dfftw_execute(plan2df)
      currentpartsc(:,:,i)=imagec_fftw
    end do Eachframe
    ! Set the gradients at initial conditions
    initialgrad: do i=1,nvar,2
      grad(i)  =dcorrvsaveragedx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,coeff,(i+1)/2,nvar)
      grad(i+1)=dcorrvsaveragedy(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,coeff,(i+1)/2,nvar)
    end do initialgrad
    do i=1,nvar,2
       write (*,*) "xshift", coeff(i)
       write (*,*) "yshift", coeff(i+1)
     end do
    ! Run the minizer
    task = 'START'
    runminimizer: do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
      ! This is the call to the L-BFGS-B code. 
      call setulb(nvar,m,coeff,lowbound,upbound,nbd,funval,grad,factr,pgtol,&
                  wa,iwa,task,iprint,&
                  csave,lsave,isave,dsave)
      if (task(1:2).eq.'FG') then
        ! calculate function value at this point
        funval=corrvsaverage(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,coeff,nvar)
        do i=1,nvar,2
          grad(i)  =dcorrvsaveragedx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,coeff,(i+1)/2,nvar)
          grad(i+1)=dcorrvsaveragedy(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,coeff,(i+1)/2,nvar)
        end do
    !!$    ! check that the derivatives work
    !!$    do i=1,nvar,2
    !!$      coefftemp=coeff
    !!$      coefftemp(i)=coefftemp(i)+epsilon
    !!$      testderx=(corrvsaverage_motionmodel(currentpartsc,boxsize,nframes,rmax1sq,rmax2sq,coefftemp,nvar)-&
    !!$                corrvsaverage_motionmodel(currentpartsc,boxsize,nframes,rmax1sq,rmax2sq,coeff,nvar))/epsilon
    !!$      coefftemp=coeff
    !!$      coefftemp(i+1)=coefftemp(i+1)+epsilon
    !!$      testdery=(corrvsaverage_motionmodel(currentpartsc,boxsize,nframes,rmax1sq,rmax2sq,coefftemp,nvar)-&
    !!$               corrvsaverage_motionmodel(currentpartsc,boxsize,nframes,rmax1sq,rmax2sq,coeff,nvar))/epsilon
    !!$      write (*,*) "derivtest",(i+1)/2,grad(i),testderx,grad(i)/testderx,grad(i+1)/testdery
    !!$    end do
      end if
    end do runminimizer
    ! Determine the shifts from the coefficients
    do i=1,nvar,2
      xshift((i+1)/2)=coeff(i)
      yshift((i+1)/2)=coeff(i+1)
    end do
    ! Output the shifted particles
    ! write out shifts
    write (*,'(A1,A11,2A12)') " "," Frame","x-shift","y-shift"
    write (7,'(A1,A11,2A12)') "#"," Frame","x-shift","y-shift"
    do i=1,nframes
      write (7,'(I12,2F12.6)') i+framefirst-1,xshift(i)-xshift(zeroframe),yshift(i)-yshift(zeroframe)
      write (*,*)              i+framefirst-1,xshift(i)-xshift(zeroframe),yshift(i)-yshift(zeroframe)
    end do
    !-------------------------------------------------
    ! Close output files
    !-------------------------------------------------
    ! Close files
    call Imclose(1)
    close(7)
  end if Fileexists
end do Eachmovie
close(3)

end program Alignframes_lbfgsb
!-------------------------------------------------
subroutine Shiftparticle(particlein,particleout,nx,ny,xshift,yshift)
!-------------------------------------------------
implicit none
integer, parameter            :: dp=kind(0.d0)   
integer, intent (in)          :: nx,ny
real(dp), intent(in)          :: xshift,yshift
real(dp)                      :: xsh,ysh
double complex, intent (in)   :: particlein(nx/2+1,ny)
double complex, intent (out)  :: particleout(nx/2+1,ny)
real(dp)                      :: phase
double complex                :: pshft
integer                       :: kxarray,kxwrap,kyarray,kywrap ! Needed for the DO construct
real(dp), parameter           :: twopi=6.2831853071796

xsh=twopi*xshift/dble(nx)
ysh=twopi*yshift/dble(ny)
particleout=0
do kyarray = 1,ny
  if (kyarray.le.ny/2+1) then
    kywrap=kyarray-1
  else
    kywrap=kyarray-ny-1
  end if
  do kxarray = 1,nx/2+1
    kxwrap=kxarray-1
    phase=dble(xsh)*dble(kxwrap)+dble(ysh)*dble(kywrap)
    pshft=dcmplx(cos(phase),sin(phase))
    particleout(kxarray,kyarray)=particlein(kxarray,kyarray)*pshft
  end do 
end do
end subroutine Shiftparticle
!-------------------------------------------------------
subroutine Floatpart(inimage,outimage,nx,ny)
!-------------------------------------------------------
implicit none
integer, parameter                :: dp=kind(0.d0)         
integer, intent(in)               :: nx, ny
real(dp), intent(in)              ::inimage(nx,ny)
real(dp), intent(out)             ::outimage(nx,ny)
real(dp)                          :: npixels,periave
! Calculate perimeter values
npixels = 2*dble(nx)+2*(dble(ny)-2)
periave=SUM(inimage(:,1))+SUM(inimage(:,ny))+SUM(inimage(1,2:ny-1))+SUM(inimage(nx,2:ny-1))
periave=periave/npixels
! Subtract the perimeter average from the input image
outimage=inimage-periave
return
end
!-------------------------------------------------------
subroutine temperaturemask(resmask,ftsize,psize,bfactor)
!-------------------------------------------------------
! Prepare mask for considering specific resolution
! and temperature factor
! temp=EXP(-B/4d**2) ; d=psize*nx/rad so 1/d**2=rad**2/(psize**2)*(nx**2)
implicit none
integer, parameter         :: dp=kind(0.d0)  
integer                    :: ftsize
real(dp)                   :: resmask(ftsize/2+1,ftsize)
real(dp)                   :: bfactor,radiussq,temperature,psize,oneoverdsquared
integer                    :: j,i,vox,lin,kywrap

do j=-ftsize/2+1,ftsize/2,1 
  do i=0,ftsize/2,1
    radiussq=dble(i)**2+dble(j)**2
    oneoverdsquared=radiussq/((psize**2)*(dble(ftsize**2)))
    temperature=exp(-0.25*oneoverdsquared*bfactor)
    vox = i+1
    lin = j+ftsize/2                  
    if (lin<ftsize/2)  kywrap=lin+ftsize/2+1
    if (lin>=ftsize/2) kywrap=lin-ftsize/2+1
    resmask(vox,kywrap)=temperature
  end do
end do

end subroutine temperaturemask

!------------------------------------------------------------------------------------------
subroutine Cram(image1,image2,nx1,ny1,nx2,ny2)
! Place DP image1 into the array for DP image2, padding or cropping as necessary
! Image1: double precision array
! Image2: double precision array
! nx1,ny1: extent of image1 array
! nx2,ny2: extent of image2 array
!------------------------------------------------------------------------------------------
implicit none
integer, parameter         :: dp=kind(0.d0)
integer                    :: nx1,ny1,nx2,ny2
integer                    :: x1start,x1stop,x2start,x2stop
integer                    :: y1start,y1stop,y2start,y2stop
real(dp)                   :: Image1(nx1,ny1)
real(dp)                   :: Image2(nx2,ny2)


image2=0
if (nx1.eq.nx2) then
  x1start=1 ; x1stop=nx1
  x2start=1 ; x2stop=nx2
else if (nx1.gt.nx2) then
  x1start=1+(nx1-nx2)/2  ; x1stop=x1start+nx2-1
  x2start=1 ; x2stop=nx2
else ! (nx2 must be gt nx1)
  x1start=1 ; x1stop=nx1
  x2start=1+(nx2-nx1)/2  ; x2stop=x2start+nx1-1
end if
if (ny1.eq.ny2) then
  y1start=1 ; y1stop=ny1
  y2start=1 ; y2stop=ny2
else if (ny1.gt.ny2) then
  y1start=1+(ny1-ny2)/2  ; y1stop=y1start+ny2-1
  y2start=1 ; y2stop=ny2
else ! (nx2 must be gt nx1)
  y1start=1 ; y1stop=ny1
  y2start=1+(ny2-ny1)/2  ; y2stop=y2start+ny1-1
end if

image2(x2start:x2stop,y2start:y2stop)=image1(x1start:x1stop,y1start:y1stop)
end subroutine Cram

!-------------------------------------------------------------
subroutine Replaceextension_arblen(filenamein,filenameout,extension)
!-------------------------------------------------------------
! replace the extension of a file with a new extension of arbitrary length
implicit none
character(20), intent(IN)     :: extension
character(80), intent(IN)     :: filenamein
character(80), intent(OUT)    :: filenameout
integer                       :: lenroot,lenext,posext
filenameout=""
lenext=len_trim(extension)
posext=index(filenamein,".",back=.true.)
lenroot=posext-1 
filenameout(1:lenroot)=filenamein(1:lenroot) ! Transfer root of name
filenameout(lenroot+1:lenroot+1)="."
filenameout(lenroot+2:lenroot+lenext+1)=trim(extension)
filenameout=trim(filenameout)
return
end subroutine

!-------------------------------------------------------
subroutine Changepath(filenamein,filenameout,path)
!-------------------------------------------------------

implicit none
character(150), intent(IN)    :: path
character(150), intent(IN)    :: filenamein
character(150), intent(OUT)   :: filenameout
integer                       :: lenname,lenpath,pospath

filenameout=""
lenpath=len_trim(path)
pospath=index(filenamein,"/",back=.true.)
lenname=len_trim(filenamein)
filenameout(1:lenpath)=trim(path)
filenameout(lenpath+1:lenpath+1+lenname-pospath)=trim(filenamein(pospath+1:lenname))
filenameout=trim(filenameout)
return
end subroutine


!-------------------------------------------------------
double precision function corrvsaverage(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,nvar)
!-------------------------------------------------------
implicit none
integer, parameter         :: dp=kind(0.d0)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize, nframes,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq
real(dp), intent(in)       :: shifts(nvar)
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp)                   :: xsh(nframes),ysh(nframes),cphase(nframes),sphase(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos,i
real(dp)                   :: twopibyboxsize,radiussq,phase,nterms,nframesinv
real(dp)                   :: rerefcomp,imrefcomp,retestcomp,imtestcomp
double complex             :: refcomp,testcomp

nterms=0_dp
nframesinv=1/dble(nframes)
twopibyboxsize=twopi/dble(boxsize)
! Convert from shifts to x-/y-shifts
xsh=0;ysh=0
do i=1,nvar,2
  xsh((i+1)/2)=shifts(i)*twopibyboxsize
  ysh((i+1)/2)=shifts(i+1)*twopibyboxsize
end do

corrvsaverage=0
eachky: do kyarray=1,boxsize
  if (kyarray.le.boxsize/2+1) then
    kywrap=kyarray-1
  else
    kywrap=kyarray-boxsize-1
  end if
  eachkx: do kxarray=1,boxsize/2+1
    kxwrap=kxarray-1
    radiussq=dble(kxwrap)**2+dble(kywrap)**2
    if (radiussq.ge.rmax1sq.and.radiussq.le.rmax2sq) then
      rerefcomp=0
      imrefcomp=0
      do zpos=1,nframes
        phase=dble(xsh(zpos))*dble(kxwrap)+dble(ysh(zpos))*dble(kywrap)
        cphase(zpos)=cos(phase)
        sphase(zpos)=sin(phase)
        refcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)*nframesinv
        rerefcomp=rerefcomp+dble(refcomp)*cphase(zpos)-aimag(refcomp)*sphase(zpos)
        imrefcomp=imrefcomp+dble(refcomp)*sphase(zpos)+aimag(refcomp)*cphase(zpos)
      end do  
      rerefcomp=rerefcomp/dble(nframes)
      imrefcomp=imrefcomp/dble(nframes)
      do zpos=1,nframes
        testcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        retestcomp=dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos)
        imtestcomp=dble(testcomp)*sphase(zpos)+aimag(testcomp)*cphase(zpos)
        corrvsaverage=corrvsaverage+rerefcomp*retestcomp+imrefcomp*imtestcomp
        nterms=nterms+1_dp
      end do
    end if
  end do eachkx
end do eachky
corrvsaverage=-corrvsaverage/nterms

end function corrvsaverage

!-------------------------------------------------------
double precision function dcorrvsaveragedx(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,j,nvar)
!-------------------------------------------------------
implicit none
integer, parameter         :: dp=kind(0.d0)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize, nframes,j,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp), intent(in)       :: shifts(nvar)
real(dp)                   :: cphase(nframes),sphase(nframes),xsh(nframes),ysh(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos,i
real(dp)                   :: radiussq,phase
real(dp)                   :: rerefcomp,imrefcomp,retestcomp,imtestcomp
real(dp)                   :: drerefcomp,dimrefcomp,dretestcomp,dimtestcomp
double complex             :: refcomp,testcomp
real(dp)                   :: twopibyboxsize,nterms,nframesinv

nframesinv=1/dble(nframes)
nterms=0_dp
twopibyboxsize=twopi/dble(boxsize)
rerefcomp=0;imrefcomp=0;drerefcomp=0;dimrefcomp=0
xsh=0;ysh=0
do i=1,nvar,2
  xsh((i+1)/2)=shifts(i)*twopibyboxsize
  ysh((i+1)/2)=shifts(i+1)*twopibyboxsize
end do

dcorrvsaveragedx=0
eachky: do kyarray=1,boxsize
  if (kyarray.le.boxsize/2+1) then
    kywrap=kyarray-1
  else 
    kywrap=kyarray-boxsize-1
  end if
  eachkx: do kxarray=1,boxsize/2+1
    kxwrap=kxarray-1
    radiussq=dble(kxwrap)**2+dble(kywrap)**2
    if (radiussq.ge.rmax1sq.and.radiussq.le.rmax2sq) then
      rerefcomp=0
      imrefcomp=0
      dimrefcomp=0
      drerefcomp=0
      do zpos=1,nframes ! Get rerefcomp, imrefcomp, drerefcomp, and dimrefcomp
        phase=dble(xsh(zpos))*dble(kxwrap)+dble(ysh(zpos))*dble(kywrap)
        cphase(zpos)=cos(phase) ! define sin and cos of phase for later use
        sphase(zpos)=sin(phase)
        refcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)*nframesinv
        rerefcomp=rerefcomp+dble(refcomp)*cphase(zpos)-aimag(refcomp)*sphase(zpos)
        imrefcomp=imrefcomp+dble(refcomp)*sphase(zpos)+aimag(refcomp)*cphase(zpos)
        if (zpos.eq.j) then ! get derivative of rerefcomp and imrefcomp when zpos=j
          drerefcomp=-(dble(refcomp)*sphase(zpos)+aimag(refcomp)*cphase(zpos))&
                     *twopibyboxsize*kxwrap
          dimrefcomp=(dble(refcomp)*cphase(zpos)-aimag(refcomp)*sphase(zpos))&
                      *twopibyboxsize*kxwrap
        end if
      end do
      rerefcomp=rerefcomp/dble(nframes)
      imrefcomp=imrefcomp/dble(nframes)
      retestcomp=0
      imtestcomp=0
      dimtestcomp=0
      dretestcomp=0
      do zpos=1,nframes ! Get retestcomp, imtestcomp, dretestcomp, dimtestcomp, and dF/dx
        testcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        retestcomp=dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos)
        imtestcomp=dble(testcomp)*sphase(zpos)+aimag(testcomp)*cphase(zpos)
        dcorrvsaveragedx=dcorrvsaveragedx+retestcomp*drerefcomp+imtestcomp*dimrefcomp
        if (zpos.eq.j) then ! only situation where dretestcomp and dimtestcomp ne 0
          dretestcomp=(-dble(testcomp)*sphase(zpos)-aimag(testcomp)*cphase(zpos))*twopibyboxsize*kxwrap
          dimtestcomp=(dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos))*twopibyboxsize*kxwrap
          dcorrvsaveragedx=dcorrvsaveragedx+rerefcomp*dretestcomp+imrefcomp*dimtestcomp
        end if
        nterms=nterms+1_dp
      end do
    end if
  end do eachkx
end do eachky
dcorrvsaveragedx=-dcorrvsaveragedx/nterms

end function dcorrvsaveragedx
!-------------------------------------------------------
double precision function dcorrvsaveragedy(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,j,nvar)
!-------------------------------------------------------
implicit none
integer, parameter         :: dp=kind(0.d0)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize, nframes,j,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp), intent(in)       :: shifts(nvar)
real(dp)                   :: cphase(nframes),sphase(nframes),xsh(nframes),ysh(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos,i
real(dp)                   :: radiussq,phase
real(dp)                   :: rerefcomp,imrefcomp,retestcomp,imtestcomp
real(dp)                   :: drerefcomp,dimrefcomp,dretestcomp,dimtestcomp
double complex             :: refcomp,testcomp
real(dp)                   :: twopibyboxsize,nterms,nframesinv

nframesinv=1/dble(nframes)
nterms=0_dp
twopibyboxsize=twopi/dble(boxsize)
rerefcomp=0;imrefcomp=0;drerefcomp=0;dimrefcomp=0
! Convert from shifts to -x/y-shifts
xsh=0;ysh=0
do i=1,nvar,2
  xsh((i+1)/2)=shifts(i)*twopibyboxsize
  ysh((i+1)/2)=shifts(i+1)*twopibyboxsize
end do

dcorrvsaveragedy=0
eachky: do kyarray=1,boxsize
  if (kyarray.le.boxsize/2+1) then
    kywrap=kyarray-1
  else 
    kywrap=kyarray-boxsize-1
  end if
  eachkx: do kxarray=1,boxsize/2+1
    kxwrap=kxarray-1
    radiussq=dble(kxwrap)**2+dble(kywrap)**2
    if (radiussq.ge.rmax1sq.and.radiussq.le.rmax2sq) then
      rerefcomp=0
      imrefcomp=0
      dimrefcomp=0
      drerefcomp=0
      do zpos=1,nframes ! Get rerefcomp, imrefcomp, drerefcomp, and dimrefcomp
        phase=dble(xsh(zpos))*dble(kxwrap)+dble(ysh(zpos))*dble(kywrap)
        cphase(zpos)=cos(phase)
        sphase(zpos)=sin(phase)
        refcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)*nframesinv
        rerefcomp=rerefcomp+dble(refcomp)*cphase(zpos)-aimag(refcomp)*sphase(zpos)
        imrefcomp=imrefcomp+dble(refcomp)*sphase(zpos)+aimag(refcomp)*cphase(zpos)
        if (zpos.eq.j) then ! get derivative of rerefcomp and imrefcomp when zpos=j
          drerefcomp=-(dble(refcomp)*sphase(zpos)+aimag(refcomp)*cphase(zpos))&
                              *twopibyboxsize*kywrap
          dimrefcomp=(dble(refcomp)*cphase(zpos)-aimag(refcomp)*sphase(zpos))&
                              *twopibyboxsize*kywrap
        end if
      end do
      rerefcomp=rerefcomp/dble(nframes)
      imrefcomp=imrefcomp/dble(nframes)
      retestcomp=0
      imtestcomp=0
      dimtestcomp=0
      dretestcomp=0
      do zpos=1,nframes ! Get retestcomp, imtestcomp, dretestcomp, dimtestcomp, and dF/dy
        testcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        retestcomp=dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos)
        imtestcomp=dble(testcomp)*sphase(zpos)+aimag(testcomp)*cphase(zpos)
        dcorrvsaveragedy=dcorrvsaveragedy+retestcomp*drerefcomp+imtestcomp*dimrefcomp
        if (zpos.eq.j) then ! only situation where dretestcomp and dimtestcomp ne 0
          dretestcomp=(-dble(testcomp)*sphase(zpos)-aimag(testcomp)*cphase(zpos))*twopibyboxsize*kywrap
          dimtestcomp=(dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos))*twopibyboxsize*kywrap
          dcorrvsaveragedy=dcorrvsaveragedy+rerefcomp*dretestcomp+imrefcomp*dimtestcomp
        end if
        nterms=nterms+1_dp
      end do
    end if
  end do eachkx
end do eachky
dcorrvsaveragedy=-dcorrvsaveragedy/nterms

end function dcorrvsaveragedy
