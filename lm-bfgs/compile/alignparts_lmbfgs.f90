!———————————————————————————
! alignparts_lmbfgs                         Dated 1/6/15
!-------------------------------------------------------
! JLR 5/14
!     5/14 Allowed for output of multiple stacks (for Relion)
!     6/14 Elliminate single stack output, parfile output
!     6/14 Always float with a ring and normalize stddev
!     6/14 Allow aligned box and output box to differ in size
!     6/14 Neigbour smoothing
!     8/14 Replaced partial derivatives with compact form to match manuscript
!     9/14 Added 1/nframes and 1/nterms to function and derivatives to keep user set constants small
!     9/14 Added stop statement when program runs into a partially off frame particle
!     9/14 Added stop statement when input file and output file are the same
!    11/14 Fixed bug in Smoothtrajectories when zeroframe is not frame 1
!    12/14 Allow root of output stack name to be modifed
!     2/15 Allow path of input and output files to be specified
!     3/15 Fixed bug that caused normalization mask to be .true. everywhere
!     4/15 Incorporate ML estimate from exposure damage and movement
!     6/15 Added motion flag
!
!-------------------------------------------------------
PROGRAM Alignparts_lbfgsb
implicit none
integer, parameter                :: dp=kind(0.d0) 
real(dp), parameter               :: pi=3.141592654
! Stack parameters
integer                           :: nxyz(3),mxyz(3),nxyzst(3),mode,nx,ny
character(80)                     :: title
real(dp)                          :: dmin,dmax,dmean,cell(6)
! I/O
character(150)                    :: movielist,coordlist,inmovie,incoord,vectorfile,outstk,tempname
character(150)                    :: inpath,outpath
character(20)                     :: vecext,stkext,stkmod
integer                           :: io
! Counters
integer                           :: i,j,filenumb,x,y
integer                           :: boxsize,subboxsize,nframes,xs,xf,ys,yf
integer                           :: framefirst_ali,framelast_ali,framefirst_avg,framelast_avg
! integer                         :: part
! Images
real,allocatable                  :: inbox(:,:)
real(dp), allocatable             :: normimage(:,:)
real(dp), allocatable             :: floatimage(:,:)
real(dp), allocatable             :: subimage(:,:)
real(dp), allocatable             :: weightmask(:,:)
real(dp), allocatable             :: expweightsc(:,:,:)
real(dp), allocatable             :: motweightsc(:,:,:)
double complex, allocatable       :: currentpartsc(:,:,:),averagec(:,:),imagec(:,:)
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
real(dp)                          :: boxsizesqinv,subboxsizesqinv,nframesinv,radtodeg
integer                           :: boxsizeby2plus1
! temp factor
real(dp)                          :: bfactor,smooth,exaggerate
integer                           :: zeroframe
! exposure weighting
integer                           :: expflag,motflag
real(dp)                          :: expframe,akv
! Neighbour averaging
integer                           :: neiflag, maxparts
real(dp)                          :: nwsig
real(dp),allocatable              :: traj(:,:,:),straj(:,:,:)
real(dp),allocatable              :: coord(:,:)
! minimizer
real(dp)                          :: corrvsaverage_cmplx,dcorrvsaveragedx_cmplx,dcorrvsaveragedy_cmplx
real(dp),allocatable              :: shifts(:),shiftstemp(:)
! l-bfgsb-b
integer                           :: nvar ! number of variables (n)
integer,  parameter               :: m = 15, iprint = -1!0!100 (iprint=-1 suppress completely, iprint=1 more etc.)
real(dp), parameter               :: pgtol  = 0
real(dp)                          :: factr ! factr=1.0d7 is moderate, 1.0d1 is extreme, 0 is not use
character(len=60)                 :: task, csave
logical                           :: lsave(4)
integer                           :: isave(44)
real(dp)                          :: funval
real(dp)                          :: dsave(29)
integer,  allocatable             :: nbd(:), iwa(:)
real(dp), allocatable             :: lowbound(:), upbound(:), grad(:), wa(:)
! masking and normalizing
logical, allocatable              :: mask(:,:)
integer                           :: noutmask
real(dp)                          :: ri,risq,radiussq
integer                           :: invflag
! coordinate file variables
integer                           :: xpos,ypos
real(dp)                          :: density
integer                           :: xposs,xposf,yposs,yposf
integer                           :: nparts,nskipped
! Derivative test (leave in code in order to check derivatives work when making changes)
!!$real(dp)                          :: testderx,testdery
!!$real(dp), parameter               :: epsilon=0.000001

include '/usr/include/fftw3.f'

read (5,'(a)') movielist
read (5,'(a)') coordlist
read (5,*) boxsize,subboxsize,ri,psize,nsigma,rmax1,rmax2
read (5,*) motflag,expflag,akv,expframe
read (5,*) bfactor,smooth,exaggerate,zeroframe,invflag
read (5,*) neiflag,maxparts,nwsig
read (5,*) framefirst_ali,framelast_ali,framefirst_avg,framelast_avg
read (5,*) factr
read (5,'(a)') inpath
read (5,'(a)') outpath
read (5,*) stkmod
read (5,*) vecext
read (5,*) stkext

! File stream
! Stream 1 : Input movie
! Stream 2 : Input coordinatefile
! Stream 3 : Input list of movies (movielist)
! Stream 4 : Output stack of averages (outstk)
! Stream 7 : Output vector file (shift file for whole frame) 
! Stream 8 : Input list of coordinate files (coordlist)

! Subbox extraction
if (subboxsize.eq.0) subboxsize=boxsize
if (subboxsize.gt.boxsize) subboxsize=boxsize
xs=(boxsize-subboxsize)/2 + 1
ys=xs
xf=xs+subboxsize-1
yf=xf

if (framefirst_avg.lt.framefirst_ali.or.framelast_avg.gt.framelast_ali) then
  write (*,*) "STOP: The frame range used for alignment must include the frame range used in the average"
  write (*,*) framefirst_avg,framelast_avg,framefirst_ali,framelast_ali
  stop
end if
nframes=framelast_ali-framefirst_ali+1
nframesinv=1/dble(nframes)
nvar=nframes*2

open(unit=3,file=movielist,status="old")
open(unit=8,file=coordlist,status="old")

boxsizesqinv = 1/dble(boxsize**2)
subboxsizesqinv = 1/dble(subboxsize**2)
boxsizeby2plus1 = boxsize/2+1
! Restate resolution limits in terms of pixel radii
rmax1=(psize*boxsize)/rmax1
if (rmax1<1) rmax1=1.01 ! Never use the origin of the FT for alignment
rmax1sq=rmax1**2
rmax2=(psize*boxsize)/rmax2
rmax2sq=rmax2**2
radtodeg=180_dp/pi
nxyzst=(/0,0,0/)
nx=0;ny=0 !avoids uninitialized warning
! Allocate FT arrays and plan
! FFTW arrays
allocate (image_fftw(boxsize,boxsize))                           ; image_fftw=0
allocate (imagec_fftw(boxsizeby2plus1,boxsize))                  ; imagec_fftw=0
allocate (imagec(boxsizeby2plus1,boxsize))                       ; imagec=0
allocate (averagec(boxsizeby2plus1,boxsize))                     ; averagec=0
allocate (weightmask(boxsizeby2plus1,boxsize))                   ; weightmask=0
allocate (floatimage(subboxsize,subboxsize))                     ; floatimage=0
allocate (normimage(subboxsize,subboxsize))                      ; normimage=0
allocate (subimage(subboxsize,subboxsize))                       ; subimage=0
allocate (inbox(boxsize,boxsize))                                ; inbox=0
allocate (mask(subboxsize,subboxsize))                           ; mask=.false.

! Prepare inverted circular mask of radius ri
risq=ri**2
noutmask=0 ! counter for number of pixels outside of mask
write (*,*) "Preparing masks"
do j=1,subboxsize
  do i=1,subboxsize
    radiussq=dble(j-subboxsize/2)**2 + dble(i-subboxsize/2)**2
    if (radiussq.le.risq) then
      mask(i,j)=.false.
    else
      mask(i,j)=.true.
      noutmask=noutmask+1
    end if
  end do
end do

write (*,*) "Generating weightmask"
call temperaturemask(weightmask,boxsize,psize,bfactor)
! Plan FFTs
call Dfftw_plan_dft_r2c_2d(plan2df,boxsize,boxsize,image_fftw,imagec_fftw,fftw_measure)
call Dfftw_plan_dft_c2r_2d(plan2dr,boxsize,boxsize,imagec_fftw,image_fftw,fftw_measure)
write (*,*) "FTs Planned"
nparts=0
filenumb=0
Eachmovie: do 
  ! Read a filename from list
  read (3,'(A80)',iostat=io) tempname ! read and input inmovie name
  Movieexists: if (io>0) then
    write (*,*) "problem with input movie list format"
    exit
  else if (io<0) THEN
    write (*,*) "End of movie list reached"
    exit
  end if Movieexists
  call Changepath(tempname,inmovie,inpath) ! Adjust path of input inmovie name
  ! Read a filename from list
  read (8,'(A80)',iostat=io) tempname ! read an input coordinate file name
  Coordexists: if (io>0) then
    write (*,*) "problem with input coordinate list format"
    exit
  else if (io<0) THEN
    write (*,*) "End of coordinate file list reached"
    exit
  end if Coordexists
  call Changepath(tempname,incoord,inpath) ! Adjust path of input inmovie name
  write (*,*) "Working on these files:"
  write (*,*) inmovie
  write (*,*) incoord
  filenumb=filenumb+1
  ! Open input movie
  call Imopen(1,inmovie,"RO")
  call Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
  nx         = nxyz(1)
  ny         = nxyz(2)
  call Replaceextension_modname(incoord,tempname,vecext,stkmod)  ! Generate output vector file name
  call Changepath(tempname,vectorfile,outpath) ! Adjust path of output vector file name
  call Replaceextension_modname(incoord,tempname,stkext,stkmod) ! Generate output stack name
  call Changepath(tempname,outstk,outpath) ! Adjust path of output stack file
  if (trim(outstk).eq.trim(inmovie).or.trim(vectorfile).eq.trim(inmovie).or.trim(vectorfile).eq.trim(outstk)) then
    write (*,*) "STOP: An output file has the same name as input file or another output file."
    write (*,*) "You probably don't want to do this."
    stop
  end if
  ! Open input coordinate file and output vectorfield file
  open(7,file=vectorfile,status='replace')
  nxyz=(/subboxsize,subboxsize,100/)
  mxyz=nxyz
  dmin=1E6
  dmax=-1E6
  dmean=0
  ! Open output stack
  call Imopen(4,outstk,"NEW")
  call Itrhdr(4,1)
  call Ialmod(4,2) ! Change to map mode 2 in case input was map mode 0
  call Ialsiz(4,nxyz,nxyzst)
  ! allocate arrays
  if (filenumb.eq.1) then
    write (*,*) "First time through loop. Allocating arrays"
    allocate (currentpartsc(boxsizeby2plus1,boxsize,nframes))        ; currentpartsc=0
    allocate (expweightsc(boxsizeby2plus1,boxsize,nframes))          ; expweightsc=1 ! keep at 1 if not used
    allocate (motweightsc(boxsizeby2plus1,boxsize,nframes))          ; motweightsc=1 ! keep at 1 if not used
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
    ! shifts 
    allocate (shifts(nvar))                                           ; shifts=0
    allocate (shiftstemp(nvar))                                       ; shiftstemp=0
    ! Neighbour averages
    allocate (coord(maxparts,2))                                      ; coord=0
    allocate(traj(maxparts,nframes,2))                                ; traj=0
    allocate(straj(maxparts,nframes,2))                               ; straj=0  
    if(expflag.eq.1) call Genexpweights(expweightsc,boxsize,nframes,psize,expframe,akv)
  end if
  ! Prepare output stack header
  ! Set minimizer conditions
  shifts=0 ! Initial guess for minimizer
  nbd=0        ! 2=All variables have a lower and upper boundary
  ! if nbd=0 "no bounds", if 1 "lower bound only", if 2 "upper and lower", if 3 "upper only"
  ! Open input coordinate file
  open(unit=2,file=incoord,status="old")
  ! Set output stack parameters to starting values
  dmin=1d6; dmax=-1d6; dmean=0
  !---------------------------------------------------------------------------------------------------
  ! Main Loop begins
  !---------------------------------------------------------------------------------------------------
  Eachcoord: do
    read(2,*,IOstat=io) xpos,ypos,density
    if (io>0) THEN
      write (*,*) "Skipping header line"
    else if (io<0) THEN
      write (*,*) "Last coordinate reached"
      exit
    else
      xposs=xpos-boxsize/2
      xposf=xpos+boxsize/2-1
      yposs=ypos-boxsize/2
      yposf=ypos+boxsize/2-1
      nparts=nparts+1
      write (*,*) "particle at coordinate",nparts,xpos,ypos
      if (xposs.lt.0.or.xposf.gt.nx.or.yposs.lt.0.or.yposf.gt.ny) then
        write (*,*) "STOP: Coordinate out of bounds for specified boxsize at record:",nparts
        write (*,*) "This problem can cause a major mismatch between the particles in a stack and star file"
        write (*,*) "Use the program removeoutlierscoordfile.f90 from"
        write (*,*) "https://sites.google.com/site/rubinsteingroup/micrograph-utilities" 
        write (*,*) "to remove offending particle positions"
        stop 
        nskipped=nskipped+1
        nparts=nparts-1
      else
        coord(nparts,:)=(/dble(xpos),dble(ypos)/) ! record the coordinate position for later
        Eachframe: do i=1,nframes ! Read the movie of individual particle into memory
          inbox=0
          call Imposn(1,i+framefirst_ali-2,0) 
          call Irdpas(1,inbox,boxsize,boxsize,xposs,xposf,yposs,yposf)
          ! remove outlier pixels from frame
          mean=(sum(inbox))*real(boxsizesqinv)
          stddev=sqrt(sum(dble(inbox)**2-dble(mean))*boxsizesqinv)
          where (abs(dble(inbox)).gt.nsigma*stddev)
            inbox=real(mean)
          end where
          image_fftw=dble(inbox)
          call Dfftw_execute(plan2df)
          currentpartsc(:,:,i)=imagec_fftw
        end do Eachframe
        initialgrad: do i=1,nvar,2
          grad(i)  =dcorrvsaveragedx_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,(i+1)/2,nvar,smooth)
          grad(i+1)=dcorrvsaveragedy_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,(i+1)/2,nvar,smooth)
        end do initialgrad
        ! Run the minizer
        task = 'START'
        runminimizer: do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
        ! This is the call to the L-BFGS-B code. 
          call setulb(nvar,m,shifts,lowbound,upbound,nbd,funval,grad,factr,pgtol,&
                    wa,iwa,task,iprint,&
                    csave,lsave,isave,dsave)
          if (task(1:2).eq.'FG') then
            ! calculate function value at this point
            funval=corrvsaverage_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,nvar,smooth)
            do i=1,nvar,2
              !write (*,*) coeff(i),coeff(i+1)
              grad(i)  =dcorrvsaveragedx_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,(i+1)/2,nvar,smooth)
              grad(i+1)=dcorrvsaveragedy_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,(i+1)/2,nvar,smooth)
            end do
            ! check that the derivatives work
!!$            Derivtest: do i=1,nvar,2
!!$              shiftstemp=shifts
!!$              shiftstemp(i)=shiftstemp(i)+epsilon
!!$              testderx=(corrvsaverage_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shiftstemp,nvar,smooth)-&
!!$                       corrvsaverage_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,nvar,smooth))/epsilon
!!$              shiftstemp=shifts
!!$              shiftstemp(i+1)=shiftstemp(i+1)+epsilon
!!$              testdery=(corrvsaverage_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shiftstemp,nvar,smooth)-&
!!$                        corrvsaverage_cmplx(currentpartsc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,nvar,smooth))/epsilon
!!$              write (*,*) "derivtest",(i+1)/2,grad(i),testderx,grad(i)/testderx,grad(i+1)/testdery
!!$            end do Derivtest
          end if
        end do runminimizer
        ! Determine the shifts from the coefficients
        do i=1,nvar,2
          xshift((i+1)/2)=shifts(i)
          yshift((i+1)/2)=shifts(i+1)
        end do
        xshift=xshift-xshift(zeroframe) ! Define one frame's shift as zero
        yshift=yshift-yshift(zeroframe)
        write (*,'(A1,A11,2A12)') " "," Frame","x-shift","y-shift"
        ! just store the trajectories for later processing
        do i=1,nframes 
          traj(nparts,i,1)=xshift(i)
          traj(nparts,i,2)=yshift(i)
          write (*,'(I12,2F12.6)') i+framefirst_ali-1,xshift(i),yshift(i)
        end do
      end if ! Coordinate is acceptable
    end if ! There is a coordinate
  end do Eachcoord
  !---------------------------------------------------------------------------------------------------
  ! Main Loop ends
  !---------------------------------------------------------------------------------------------------
  ! Smooth trajectories if requested begins
  !---------------------------------------------------------------------------------------------------
  if (neiflag.eq.1) then
    call Smoothtrajectories(traj,straj,coord,maxparts,nparts,nframes,nwsig)
  else
    straj=traj
  end if 
  do i=1,nparts
    if(motflag.eq.1) call Genmotweights(motweightsc,boxsize,nframes,straj(i,:,:))
    write (*,*) "particle",i,"applying smoothed shifts"
    averagec=0
    do j=1,nframes
      write (7,'(I12,2F12.6)') j+framefirst_ali-1,exaggerate*straj(i,j,1)+coord(i,1),exaggerate*straj(i,j,2)+coord(i,2)
      ! only include frames the are within the specified range
      if (j.ge.framefirst_avg-framefirst_ali+1.and.j.le.nframes-(framelast_ali-framelast_avg)) then
        xpos=INT(coord(i,1));ypos=INT(coord(i,2)) ! re-read the saved coordinate for this particle
        xposs=xpos-boxsize/2
        xposf=xpos+boxsize/2-1
        yposs=ypos-boxsize/2
        yposf=ypos+boxsize/2-1
        ! Re-read the particle from input movie
        call Imposn(1,j+framefirst_ali-2,0) 
        call Irdpas(1,inbox,boxsize,boxsize,xposs,xposf,yposs,yposf)
        mean=(sum(inbox))*real(boxsizesqinv)
        stddev=sqrt(sum(dble(inbox)**2-dble(mean))*boxsizesqinv)
        where (abs(dble(inbox)).gt.nsigma*stddev)
          inbox=real(mean)
        end where
        image_fftw=dble(inbox)
        call Dfftw_execute(plan2df)
        imagec(:,:)=imagec_fftw
        call Shiftparticle(imagec,currentpartsc(:,:,j),boxsize,boxsize,straj(i,j,1),straj(i,j,2))
        do y=1,boxsize
          do x=1,boxsizeby2plus1
            ! Full ML correction approach
            !averagec(x,y)=sum(currentpartsc(x,y,:)*motweightsc(x,y,:)*expweightsc(x,y,:))&
            !             /sum((expweightsc(x,y,:)*motweightsc(x,y,:))**2)
            !weighted average approach
            averagec(x,y)=sum(motweightsc(x,y,:)*expweightsc(x,y,:)*currentpartsc(x,y,:))/&
                          sum(abs(motweightsc(x,y,:)*expweightsc(x,y,:)))
          end do
        end do        
      end if
    end do
    imagec_fftw=averagec*nframesinv
    call Dfftw_execute(plan2dr) 
    image_fftw=image_fftw*boxsizesqinv
    ! Prepare particle Relion style
    if (invflag.eq.1) then
      image_fftw=-image_fftw
      write (*,'(A11)',advance='no') "Inverting;" 
    end if
    subimage(:,:)=image_fftw(xs:xf,ys:yf) ! Extract central part of the image (boxsize vs subboxsize)
    call Floatring(subimage,floatimage,mask,noutmask,subboxsize,subboxsize) ! Float the image
    call Normstddev(floatimage,normimage,mask,noutmask,subboxsize,subboxsize) ! Normalize the noise
    if (minval(normimage).lt.dmin) dmin=minval(normimage)
    if (maxval(normimage).gt.dmax) dmax=maxval(normimage)
    dmean=dmean+real(sum(normimage)*boxsizesqinv)
    call Imposn(4,i-1,0) ! write out particle (no neigbor smoothing)
    call Iwrsec(4,real(normimage))    
    write (7,*) "";write (7,*) "" ! Add two blank lines to end of vector file to help gnuplot
  end do
  !---------------------------------------------------------------------------------------------------
  ! Smooth trajectories if requested ends
  !---------------------------------------------------------------------------------------------------
  call Imclose(1) ! close input movie
  close(7) ! output vector file
  ! Set output header information 
  dmean=dmean/dble(nparts)
  nxyz=(/subboxsize,subboxsize,nparts/)
  cell=(/real(psize)*real(subboxsize),real(psize)*real(subboxsize),real(nparts),90.0,90.0,90.0/) 
  mxyz=nxyz
  call Ialsiz(4,nxyz,nxyzst)
  call Ialcel(4,cell)
  call Ialsam(4,mxyz)
  title="Aligned particles from alignparts_lbfgsb"
  call Iwrhdr(4,title,0,REAL(dmin),REAL(dmax),REAL(dmean))
  call Imclose(4)
  nparts=0
end do Eachmovie
close(3) ! List of movies
end program Alignparts_lbfgsb
!===========================================================
subroutine Shiftparticle(particlein,particleout,nx,ny,xshift,yshift)
!===========================================================
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
real(dp), parameter               :: twopi=6.2831853071796

xsh=twopi*xshift/REAL(nx)
ysh=twopi*yshift/REAL(ny)
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
!===========================================================
subroutine Floatpart(inimage,outimage,nx,ny)
!===========================================================
implicit none
integer, parameter      :: dp=kind(0.d0)         
integer, intent(in)     :: nx,ny
real(dp), intent(in)    :: inimage(nx,ny)
real(dp), intent(out)   :: outimage(nx,ny)
real(dp)                :: npixels,periave
! Calculate perimeter values
npixels = 2*real(nx)+2*(real(ny)-2)
periave=sum(inimage(:,1))+sum(inimage(:,ny))+sum(inimage(1,2:ny-1))+sum(inimage(nx,2:ny-1))
periave=periave/npixels
! Subtract the perimeter average from the input image
outimage=inimage-periave
return
end

!===========================================================
double precision function corrvsaverage_cmplx(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,nvar,smooth)
!===========================================================
implicit none
integer, parameter         :: dp=kind(0.d0)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize,nframes,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq,smooth
real(dp), intent(in)       :: shifts(nvar)
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp)                   :: xsh(nframes),ysh(nframes),cphase(nframes),sphase(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos,i
real(dp)                   :: twopibyboxsize,radiussq,phase,nframesinv,nterms
real(dp)                   :: rerefcomp,imrefcomp,retestcomp,imtestcomp
double complex             :: refcomp,testcomp

nterms=0
twopibyboxsize=twopi/dble(boxsize)
nframesinv=1/dble(nframes)
! Convert from shifts to x-/y-shifts
xsh=0;ysh=0
do i=1,nvar,2
  xsh((i+1)/2)=shifts(i)*twopibyboxsize
  ysh((i+1)/2)=shifts(i+1)*twopibyboxsize
end do

corrvsaverage_cmplx=0
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
      do zpos=1,nframes
        testcomp=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        retestcomp=dble(testcomp)*cphase(zpos)-aimag(testcomp)*sphase(zpos)
        imtestcomp=dble(testcomp)*sphase(zpos)+aimag(testcomp)*cphase(zpos)
        corrvsaverage_cmplx=corrvsaverage_cmplx+rerefcomp*retestcomp+imrefcomp*imtestcomp
        nterms=nterms+1_dp
      end do
    end if
  end do eachkx
end do eachky
corrvsaverage_cmplx=-corrvsaverage_cmplx*(1/dble(nterms))

! Smoothing constraint model 2nd order
do i=3,nframes
 corrvsaverage_cmplx=corrvsaverage_cmplx+smooth*(xsh(i)-2*xsh(i-1)+xsh(i-2))**2+&
               smooth*(ysh(i)-2*ysh(i-1)+ysh(i-2))**2
end do
end function corrvsaverage_cmplx

!===========================================================
double precision function dcorrvsaveragedx_cmplx(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,a,nvar,smooth)
!===========================================================
implicit none
integer, parameter         :: dp=kind(0.d0)
doublecomplex, parameter   :: i=(0,1)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize,nframes,a,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq,smooth
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp), intent(in)       :: shifts(nvar)
real(dp)                   :: xsh(nframes),ysh(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos
real(dp)                   :: radiussq,phase
double complex             :: Fjz,Fjz_conj,Sjz,Sjz_conj,FjaSja,Fja_conjSja_conj,sum1,sum2,funk
real(dp)                   :: twopibyboxsize,nframesinv,nterms

nterms=0
nframesinv=1/dble(nframes)
twopibyboxsize=twopi/dble(boxsize) 
xsh=0;ysh=0
do zpos=1,nvar,2
  xsh((zpos+1)/2)=shifts(zpos)*twopibyboxsize
  ysh((zpos+1)/2)=shifts(zpos+1)*twopibyboxsize
end do

dcorrvsaveragedx_cmplx=0
funk=0
FjaSja=0;Fja_conjSja_conj=0
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
      Fjz=0
      Fjz_conj=0
      Sjz=0
      Sjz_conj=0
      sum1=0
      sum2=0
      do zpos=1,nframes ! Get rerefcomp, imrefcomp, drerefcomp, and dimrefcomp
        phase=dble(xsh(zpos))*dble(kxwrap)+dble(ysh(zpos))*dble(kywrap)
        Sjz=dcmplx(cos(phase),sin(phase))
        Sjz_conj=conjg(Sjz)
        Fjz=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        Fjz_conj=conjg(Fjz)
        if (zpos.eq.a) then
          FjaSja=Fjz*Sjz
          Fja_conjSja_conj=Fjz_conj*Sjz_conj
        end if
        sum1=sum1+Fjz*Sjz
        sum2=sum2+Fjz_conj*Sjz_conj
      end do
      funk=funk+twopibyboxsize*i*kxwrap*(FjaSja*sum2-Fja_conjSja_conj*sum1)
      nterms=nterms+1_dp
    end if
  end do eachkx
end do eachky
dcorrvsaveragedx_cmplx=-dble(funk)*nframesinv*(1/dble(nterms))*nframesinv


! Second order smoothing (simplified)
if (a.eq.1) then
  dcorrvsaveragedx_cmplx=dcorrvsaveragedx_cmplx+&
  2*smooth*twopibyboxsize*(xsh(a+2)-2*xsh(a+1)+xsh(a))
else if (a.eq.2) then
  dcorrvsaveragedx_cmplx=dcorrvsaveragedx_cmplx+&
  2*smooth*twopibyboxsize*(-4*xsh(a+1)+5*xsh(a)-2*xsh(a-1)+xsh(a+2))
else if (a.ge.3.and.a.le.nframes-2) then
  dcorrvsaveragedx_cmplx=dcorrvsaveragedx_cmplx+&
  2*smooth*twopibyboxsize*(6*xsh(a)-4*xsh(a-1)+xsh(a-2)-4*xsh(a+1)+xsh(a+2))
else if (a.eq.nframes-1) then
  dcorrvsaveragedx_cmplx=dcorrvsaveragedx_cmplx+&
  2*smooth*twopibyboxsize*(5*xsh(a)-4*xsh(a-1)+xsh(a-2)-2*xsh(a+1))
else if (a.eq.nframes) then
  dcorrvsaveragedx_cmplx=dcorrvsaveragedx_cmplx+&
  2*smooth*twopibyboxsize*(xsh(a)-2*xsh(a-1)+xsh(a-2))
end if

end function dcorrvsaveragedx_cmplx

!===========================================================
double precision function dcorrvsaveragedy_cmplx(particlesc,weightmask,boxsize,nframes,rmax1sq,rmax2sq,shifts,a,nvar,smooth)
!===========================================================
implicit none
integer, parameter         :: dp=kind(0.d0)
doublecomplex, parameter   :: i=(0,1)  
real(dp), parameter        :: twopi=6.2831853071796
integer, intent(in)        :: boxsize,nframes,a,nvar
real(dp), intent(in)       :: rmax1sq,rmax2sq,smooth
double complex, intent(in) :: particlesc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)       :: weightmask(boxsize/2+1,boxsize)
real(dp), intent(in)       :: shifts(nvar)
real(dp)                   :: xsh(nframes),ysh(nframes)
integer                    :: kxarray,kyarray,kxwrap,kywrap,zpos
real(dp)                   :: radiussq,phase
double complex             :: Fjz,Fjz_conj,Sjz,Sjz_conj,FjaSja,Fja_conjSja_conj,sum1,sum2,funk
real(dp)                   :: twopibyboxsize,nframesinv,nterms

nterms=0
nframesinv=1/dble(nframes)
twopibyboxsize=twopi/dble(boxsize) 
xsh=0;ysh=0
do zpos=1,nvar,2
  xsh((zpos+1)/2)=shifts(zpos)*twopibyboxsize
  ysh((zpos+1)/2)=shifts(zpos+1)*twopibyboxsize
end do

dcorrvsaveragedy_cmplx=0
funk=0
FjaSja=0;Fja_conjSja_conj=0
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
      Fjz=0
      Fjz_conj=0
      Sjz=0
      Sjz_conj=0
      sum1=0
      sum2=0
      do zpos=1,nframes ! Get rerefcomp, imrefcomp, drerefcomp, and dimrefcomp
        phase=dble(xsh(zpos))*dble(kxwrap)+dble(ysh(zpos))*dble(kywrap)
        Sjz=dcmplx(cos(phase),sin(phase))
        Sjz_conj=conjg(Sjz)
        Fjz=particlesc(kxarray,kyarray,zpos)*weightmask(kxarray,kyarray)
        Fjz_conj=conjg(Fjz)
        if (zpos.eq.a) then
          FjaSja=Fjz*Sjz
          Fja_conjSja_conj=Fjz_conj*Sjz_conj
        end if
        sum1=sum1+Fjz*Sjz
        sum2=sum2+Fjz_conj*Sjz_conj
      end do
      funk=funk+twopibyboxsize*i*kywrap*(FjaSja*sum2-Fja_conjSja_conj*sum1)
      nterms=nterms+1_dp
    end if
  end do eachkx
end do eachky
dcorrvsaveragedy_cmplx=-dble(funk)*nframesinv*(1/dble(nterms))*nframesinv

! Second order smoothing (simplified)
if (a.eq.1) then
  dcorrvsaveragedy_cmplx=dcorrvsaveragedy_cmplx+&
  2*smooth*twopibyboxsize*(ysh(a+2)-2*ysh(a+1)+ysh(a))
else if (a.eq.2) then
  dcorrvsaveragedy_cmplx=dcorrvsaveragedy_cmplx+&
  2*smooth*twopibyboxsize*(-4*ysh(a+1)+5*ysh(a)-2*ysh(a-1)+ysh(a+2))
else if (a.ge.3.and.a.le.nframes-2) then
  dcorrvsaveragedy_cmplx=dcorrvsaveragedy_cmplx+&
  2*smooth*twopibyboxsize*(6*ysh(a)-4*ysh(a-1)+ysh(a-2)-4*ysh(a+1)+ysh(a+2))
else if (a.eq.nframes-1) then
  dcorrvsaveragedy_cmplx=dcorrvsaveragedy_cmplx+&
  2*smooth*twopibyboxsize*(5*ysh(a)-4*ysh(a-1)+ysh(a-2)-2*ysh(a+1))
else if (a.eq.nframes) then
  dcorrvsaveragedy_cmplx=dcorrvsaveragedy_cmplx+&
  2*smooth*twopibyboxsize*(ysh(a)-2*ysh(a-1)+ysh(a-2))
end if

end function dcorrvsaveragedy_cmplx

!===========================================================
subroutine temperaturemask(resmask,ftsize,psize,bfactor)
!===========================================================
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

!===========================================================
subroutine Changepath(filenamein,filenameout,path)
!===========================================================

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


!===========================================================
subroutine Replaceextension_modname(filenamein,filenameout,extension,flag)
!===========================================================
! replace the extension of a file with a new extension of arbitrary length
implicit none
character(20), intent(IN)     :: extension,flag
character(150), intent(IN)    :: filenamein
character(150), intent(OUT)   :: filenameout
integer                       :: lenroot,lenext,lenflag,posext

filenameout=""
lenext=len_trim(extension)
lenflag=len_trim(flag)
posext=index(filenamein,".",back=.true.)
lenroot=posext-1

if (lenflag.ne.0) then
  if (flag(lenflag:lenflag).eq.".") then
    lenflag=lenflag-1
  end if
end if

filenameout(1:lenroot)=filenamein(1:lenroot) ! Transfer root of name
filenameout(lenroot+1:lenroot+lenflag+1)=trim(flag)
filenameout(lenroot+lenflag+1:lenroot+lenflag+1)="."
filenameout(lenroot+lenflag+2:lenroot+lenext+lenflag+1)=trim(extension)
filenameout=trim(filenameout)
return
end subroutine

!===========================================================
subroutine Removeflag(filenamein,filenameout,flag)
!===========================================================
! remove the flag the extension of a file with a new extension of arbitrary length
implicit none
character(20), intent(IN)     :: flag
character(150), intent(IN)    :: filenamein
character(150), intent(OUT)   :: filenameout
integer                       :: lenin,lenflag,posflag

filenameout=""
lenflag=len_trim(flag)
lenin=len_trim(filenamein)
posflag=index(trim(filenamein),trim(flag),back=.true.)
write (*,*) trim(filenamein)
write (*,*) trim(flag)

if (posflag.eq.0) then
  write (*,'(2a)') "The input movie name does not contain the flag ",trim(flag)
  write (*,'(a)') "that you have requested this program remove to generate the output filename"
  stop
else
  filenameout(1:posflag-1)=filenamein(1:posflag-1)
  filenameout(posflag:lenin-lenflag)=filenamein(posflag+lenflag:lenin)
end if
return
end subroutine

!===========================================================
subroutine Normstddev(inimage,outimage,mask,noutmask,nx,ny)
!===========================================================
! Set the standard deviation in the area outside of the mask to 1
implicit none
integer, parameter     :: dp=kind(0.d0)         
integer                :: nx, ny
real(dp), intent(in)   :: inimage(nx,ny)
real(dp), intent(out)  :: outimage(nx,ny)
integer, intent(in)    :: noutmask
real(dp)               :: mean,stddev,noutmaskinv
logical, intent(in)    :: mask(nx,ny)

noutmaskinv=1/dble(noutmask)
mean=(sum(inimage,mask))*noutmaskinv
stddev=sqrt(sum((inimage-mean)**2,mask)*noutmaskinv)
write (*,'(A10,F7.3,A1)',advance='yes') "Stddev in",stddev,";"
outimage=inimage/stddev

return
end

!===========================================================
subroutine Floatring(inimage,outimage,mask,noutmask,nx,ny)
!===========================================================
! Float the image by the average pixel value outside of the mask
implicit none
integer, parameter    :: dp=kind(0.d0)     
integer               :: nx, ny
real(dp), intent(in)  :: inimage(nx,ny)
real(dp), intent(out) :: outimage(nx,ny)
integer, intent(in)   :: noutmask
real(dp)              :: meanoutmask,noutmaskinv
logical, intent(in)   :: mask(nx,ny)

noutmaskinv=1/dble(noutmask)
meanoutmask=sum(inimage,mask)*noutmaskinv
write (*,'(A9,F15.3,A1)',advance='no') "Float by",meanoutmask,";"
outimage=inimage-meanoutmask
return
end

!===========================================================
subroutine Smoothtrajectories(traj,straj,coord,maxparts,nparts,nframes,sigma)
!===========================================================
implicit none
integer, parameter    :: dp=kind(0.d0)  
integer, intent(in)   :: maxparts,nframes,nparts
real(dp),intent(in)   :: traj(maxparts,nframes,2)
real(dp),intent(in)   :: coord(maxparts,2)
real(dp),intent(out)  :: straj(maxparts,nframes,2)
integer               :: i,j,k
real(dp)              :: distsq,wght,twosigmasq,sigma
real(dp)              :: wghtsum

twosigmasq=2*sigma**2
eachpart: do i=1,nparts
  wghtsum=0
  eachcomp: do j=1,nparts
    distsq=(coord(i,1)-coord(j,1))**2+(coord(i,2)-coord(j,2))**2
    wght=exp(-distsq/(twosigmasq))
    !write (*,*) i,j,wght,distsq,-distsq/(twosigmasq),"weight"
    !straj(i,1,1)=0; straj(i,1,2)=0 
    eachframe: do k=1,nframes
      straj(i,k,1)= straj(i,k,1)+wght*traj(j,k,1)
      straj(i,k,2)= straj(i,k,2)+wght*traj(j,k,2)
    end do eachframe
    wghtsum=wghtsum+wght
  end do eachcomp
  straj(i,:,1)=straj(i,:,1)/wghtsum ! finish weigting and add back starting position
  straj(i,:,2)=straj(i,:,2)/wghtsum
end do eachpart
end subroutine

!===========================================================
subroutine Genexpweights(expweightsc,boxsize,nframes,psize,expframe,akv)
!===========================================================
implicit none
integer, parameter    :: dp=kind(0.d0)  
integer               :: boxsize, nframes
real(dp)              :: expweightsc(boxsize/2+1,boxsize,nframes)
real(dp)              :: psize, expframe, akv
integer               :: frame,i,j,vox,lin,kywrap
real(dp)              :: N,Ne,freq,radiussq,oneoverboxtimespix

write (*,*) "Generating Exposure weights"
write (*,*) "The curves used to perform exposure weighting were measured by"
write (*,*) "Drs. Tim Grant and Niko Grigorieff and the HHMI Janelia Research Campus"
write (*,*) "We gratefully acknowledge Drs. Grant and Grigorieff for providing this"
write (*,*) "unpublished data"

! Convert to equivalent exposure
if (akv.eq.100) then
  expframe=expframe/0.68
else if (akv.eq.120) then
  expframe=expframe/0.75 ! Just a guesss
else if (akv.eq.200) then
  expframe=expframe
else if (akv.eq.300) then
  expframe=expframe/1.22
else if (akv.eq.400) then
  expframe=expframe/1.33
else
  write (*,*) "STOP: Equivalent exposure for that accelerating voltage not known"
  stop
end if

oneoverboxtimespix=1/(dble(boxsize)*psize)

do frame=1,nframes
  N=expframe*dble(frame)
  ! Prepare mask for considering specific resolution
  do j=-boxsize/2+1,boxsize/2,1 
    do i=0,boxsize/2,1
      radiussq=dble(i)**2+dble(j)**2
      freq=sqrt(radiussq)*oneoverboxtimespix
      Ne=0.67586*freq**(-1.31489) ! Data from Tim Grant and Niko Grigorieff
      vox = i+1
      lin = j+boxsize/2                  
      if (lin<boxsize/2)  kywrap=lin+boxsize/2+1
      if (lin>=boxsize/2) kywrap=lin-boxsize/2+1
      expweightsc(vox,kywrap,frame)=exp(-N/(2*Ne))
    end do
  end do
end do
return
end subroutine


!===========================================================
subroutine Genmotweights(motweightsc,boxsize,nframes,traj)
!===========================================================
implicit none
integer, parameter        :: dp=kind(0.d0)
real(dp), parameter       :: pi = (3.1415926535897)
integer, intent(in)       :: boxsize, nframes
real(dp), intent(out)     :: motweightsc(boxsize/2+1,boxsize,nframes)
real(dp), intent(in)      :: traj(nframes,2)     
integer                   :: frame,i,j,vox,lin,kywrap
real(dp)                  :: fullrotbyboxsize,xsh,ysh,phi
real(dp), dimension(2000) :: alut=[(2*sin((dble(j)*pi/180)/2)/(dble(j)*pi/180),j=1,2000)]
! JLR 4/15 determines the weighting factor for Fourier components measured with a known drift during
!          the frame. Factor is (2*sin(phi/2))/phi, where phi is the phase change of the component
!          when phi=0 the factor is exactly 1. Position 1 of array alut corresponds to phi=0
!          Each position in alut corresponds to 1 degree (not radians)
!twopibyboxsize=2*pi/dble(boxsize)
fullrotbyboxsize=360.0/dble(boxsize) ! specify phi in degrees, rather than radians

do frame=1,nframes
  if (frame.lt.nframes) then
    xsh=(traj(frame+1,1)-traj(frame,1))*fullrotbyboxsize
    ysh=(traj(frame+1,2)-traj(frame,2))*fullrotbyboxsize
  else
    xsh=(traj(frame,1)-traj(frame-1,1))*fullrotbyboxsize
    ysh=(traj(frame,2)-traj(frame-1,2))*fullrotbyboxsize
  end if
  do j=-boxsize/2+1,boxsize/2,1 
    do i=0,boxsize/2,1     
      vox = i+1
      lin = j+boxsize/2                  
      if (lin<boxsize/2)  kywrap=lin+boxsize/2+1
      if (lin>=boxsize/2) kywrap=lin-boxsize/2+1
      phi=xsh*dble(i)+ysh*dble(j)
      if (abs(phi).lt.1) then
        motweightsc(vox,kywrap,frame)=1
      else if (abs(phi).gt.2000) then
        motweightsc(vox,kywrap,frame)=1d-6
      else
        motweightsc(vox,kywrap,frame)=alut(nint(abs(phi)))
      end if
    end do
  end do
end do

return
end subroutine

