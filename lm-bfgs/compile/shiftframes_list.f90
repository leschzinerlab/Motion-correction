!-------------------------------------------------------
! shiftframes_list                          Dated 9/2/15
!-------------------------------------------------------
! JLR 1/14 apply precalcuated shifts to movie frames
!     2/15 Allow input/output paths to be specified, rootname to be modified
!-------------------------------------------------------

program Shiftframes
implicit none
integer, parameter                :: dp=kind(0.d0) 
integer                           :: nx,ny,nz,ftsize
integer                           :: nxyz(3),nxyzst(3),mxyz(3),mode
real                              :: dmin,dmax,dmean
real                              :: dminali,dmaxali,dmeanali
real                              :: dminavg,dmaxavg,dmeanavg
! real images and movies
real,allocatable                  :: image(:,:)
real(dp),allocatable              :: holder(:,:)
integer                           :: sec,framefirst,framelast,framefirsti,framelasti
double complex, allocatable       :: shiftedimagec(:,:)
double complex, allocatable       :: avgimagec(:,:)
! input/output files
character(150)                    :: inlist,inmovie,shiftlist,shiftfile,outmovie,outavg,tempname
character(150)                    :: inpath,outpath
character(80)                     :: title
character(20)                     :: aliext,avgext,stkmod ! extentsion for aligned movies
integer                           :: aliflag,avgflag
integer                           :: io
! image statistics
real                              :: nsigma,stddev,mean
! optimization
real                              :: nxnyinv,nzinv
double precision                  :: ftsizesqinv
! shifts
integer                           :: shiftsec
real                              :: xshiftr,yshiftr
! fftw arrays
real(dp), allocatable             :: image_fftw(:,:)
double complex, allocatable       :: imagec_fftw(:,:)
real(dp)                          :: plan2dr,plan2df

include '/usr/include/fftw3.f'

read (5,*) ftsize,nsigma
read (5,'(a)') inlist ! list of movies to be aligned
read (5,'(a)') shiftlist ! corresponding list of shift files
read (5,'(a)') inpath
read (5,'(a)') outpath
read (5,*) stkmod
read (5,*) aliflag,aliext
read (5,*) avgflag,avgext
read (5,*) framefirsti,framelasti

write (*,*) "infile list",inlist

open(unit=4,file=inlist,status="old")
open(unit=7,file=shiftlist,status="old")

! data streams
! stream 1: input movie (inmovie)
! stream 2: output aligned movie (outmovie)
! stream 3: input shifts file (inshifts)
! stream 4: input list of movies (inlist)
! stream 5: reserved for executing script
! stream 6: reserved for screen output
! stream 7: input list of shifts
! stream 8, averaged image

! allocate fourier transform arrays
allocate (image_fftw(ftsize,ftsize))        ; image_fftw=0
allocate (imagec_fftw(ftsize/2+1,ftsize))   ; imagec_fftw=0
allocate (shiftedimagec(ftsize/2+1,ftsize)) ; shiftedimagec=0 
if (avgflag.eq.1) allocate (avgimagec(ftsize/2+1,ftsize))     ; avgimagec=0 

! plan ffts
write (*,*) "planning fourier transforms. size:",ftsize
call dfftw_plan_dft_r2c_2d(plan2df,ftsize,ftsize,image_fftw,imagec_fftw,fftw_estimate)
call dfftw_plan_dft_c2r_2d(plan2dr,ftsize,ftsize,imagec_fftw,image_fftw,fftw_estimate)

ftsizesqinv=1/(dble(ftsize)**2)
fileloop: do
  ! read a filename from list
  read (4,*,iostat=io) tempname
  ! determine input image name
  if (io>0) then
    write (*,*) "problem with input format"
    exit
  else if (io<0) then
    write (*,*) "end of file list reached"
    exit
  else
    write (*,*) "inmovie:",inmovie
    call Changepath(tempname,inmovie,inpath)  ! Adjust path of tempname to get input inmovie name
    call Changepath(inmovie,tempname,outpath) ! Adjust path of input movie name to get tempname for output files
    call Replaceextension_modname(tempname,outmovie,aliext,stkmod)
    call Replaceextension_modname(tempname,outavg,avgext,stkmod)
    write (*,*) "output aligned movie:",outmovie
    ! open input movie file
    call Imopen (1,inmovie,"old")
    call Irdhdr(1,nxyz,mxyz,mode,dmin,dmax,dmean)
    nx     = nxyz(1)
    ny     = nxyz(2)
    nz     = nxyz(3)
    nzinv=1/real(nz)
    nxnyinv=1/(real(nx)*real(ny))
    if (framelasti.eq.0) then
      framefirst=framefirsti; framelast=nz
    else
      framefirst=framefirsti;framelast=framelasti
    end if
    ! allocate images unique to this run
    allocate (image(nx,ny))                   ; image=0
    allocate (holder(nx,ny))                  ; holder=0
    ! open input shifts file and skip header
    read (7,*,iostat=io) shiftfile
    open(unit=3,file=shiftfile,status="old")
    do 
      read(3,'(i12,2f12.6)',iostat=io) shiftsec,xshiftr,yshiftr
      if (io.lt.0) then
        write (*,*) "reached eof while skipping shiftfile header. check parfile format"
        stop
      end if
      if (io.eq.0) then
        backspace(3)
        exit
      end if
    end do
    ! open output movie if requested
    if (aliflag.eq.1) then
      call Imopen(2,outmovie,"unknown") 
      call Itrhdr(2,1)
      dminali=0;dmaxali=0;dmeanali=0
    end if
    if (avgflag.eq.1) then
      call imopen(8,outavg,"unknown") 
      call itrhdr(8,1)
      dminavg=0;dmaxavg=0;dmeanavg=0
      avgimagec=0
    end if
    eachsection: do sec=framefirst,framelast
      ! read the shift values from the shift file until correct section is reached
      do      
        read (3,'(i12,2f12.6)') shiftsec,xshiftr,yshiftr
        if (shiftsec.eq.sec) exit
      end do
      call Imposn(1,sec-1,0)
      call Irdpas(1,image(:,:),nx,ny,0,nx-1,0,ny-1)
      mean=sum(image)*nxnyinv
      image=image-mean
      stddev=sqrt(sum(image**2)*nxnyinv) ! mean has been set to 0
      where (abs(image).gt.nsigma*stddev)
        image=0
      end where
      write (*,*) "working on frame:",sec,xshiftr,yshiftr
      call Cram(dble(image(:,:)),image_fftw,nx,ny,ftsize,ftsize)
      call Dfftw_execute(plan2df)
      call Shiftframe(imagec_fftw(:,:),shiftedimagec(:,:),ftsize,xshiftr,yshiftr)
      if (avgflag.eq.1) avgimagec=avgimagec+shiftedimagec ! accumulate average if requested
      ! write out shifted frame if requested
      if (aliflag.eq.1) then
        imagec_fftw=shiftedimagec
        call Dfftw_execute(plan2dr)
        call Cram(image_fftw*ftsizesqinv,holder,ftsize,ftsize,nx,ny)
        call Imposn(2,sec-1,0)
        call Iwrsec(2,real(holder))
        if (real(minval((holder))).lt.dminali) dminali=real(minval((holder)))
        if (real(maxval((holder))).gt.dmaxali) dmaxali=real(maxval((holder)))
        dmeanali=real(sum(holder))*nxnyinv
      end if
    end do eachsection
    if (aliflag.eq.1) then
      dmeanali=dmeanali*nzinv
      nxyz=(/nx,ny,framelast-framefirst+1/);nxyzst=(/0,0,0/)
      mxyz=nxyz
      call Ialsiz(2,nxyz,nxyzst)
      call Ialsam(2,mxyz)
      title="frames shifted by shiframes_list.f90"
      call Iwrhdr(2,title,0,dminali,dmaxali,dmeanali)
      call Imclose(2) ! close output movie
    end if
    if (avgflag.eq.1) then
      imagec_fftw=avgimagec
      call Dfftw_execute(plan2dr)
      call Cram(image_fftw*ftsizesqinv,holder,ftsize,ftsize,nx,ny)
      call Imposn(8,0,0)
      call Iwrsec(8,real(holder))
      dminavg=real(minval((holder)))
      dmaxavg=real(maxval((holder)))
      dmeanavg=real(sum(holder))*nxnyinv
      nxyz=(/nx,ny,1/);nxyzst=(/0,0,0/)
      mxyz=nxyz
      call Ialsiz(8,nxyz,nxyzst)
      call Ialsam(8,mxyz)
      title="average of frames shifted by shiframes_list.f90"
      call Iwrhdr(8,title,0,dminavg,dmaxavg,dmeanavg)
      call Imclose(8) ! close output movie
    end if
    call Imclose(1) ! close input movie
    close(3)
    deallocate (image)
    deallocate (holder)
  end if
end do fileloop

stop
end program Shiftframes
!-------------------------------------------------
subroutine Cram(image1,image2,nx1,ny1,nx2,ny2)
! place dp image1 into the array for dp image2, padding or cropping as necessary
! image1: double precision array
! image2: double precision array
! nx1,ny1: extent of image1 array
! nx2,ny2: extent of image2 array
!-------------------------------------------------
implicit none
integer           :: nx1,ny1,nx2,ny2
integer           :: x1start,x1stop,x2start,x2stop
integer           :: y1start,y1stop,y2start,y2stop
double precision  :: image1(nx1,ny1)
double precision  :: image2(nx2,ny2)


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
image2(x2start:x2stop,y2start:y2stop)=dble(image1(x1start:x1stop,y1start:y1stop))

end subroutine Cram
!-------------------------------------------------
subroutine Shiftframe(framein,frameout,ftsize,xshift,yshift)
!-------------------------------------------------
! JLR 5/14
implicit none 
integer, parameter            :: dp=kind(0.d0) 
integer, intent (IN)          :: ftsize
real, intent (IN)             :: xshift,yshift
double complex, intent (IN)   :: framein(ftsize/2+1,ftsize)
double complex, intent (OUT)  :: frameout(ftsize/2+1,ftsize)
real(dp)                      :: xsh,ysh,phase
double complex                :: pshft
integer                       :: kxarray,kxwrap,kyarray,kywrap ! needed for the do construct
real(dp), parameter           :: twopi=6.2831853071796

xsh=twopi*dble(xshift)/dble(ftsize)
ysh=twopi*dble(yshift)/dble(ftsize)

do kyarray = 1,ftsize
  if (kyarray.le.ftsize/2+1) kywrap=kyarray-1
  if (kyarray.ge.ftsize/2+2) kywrap=kyarray-ftsize-1
  do kxarray = 1,ftsize/2+1
    kxwrap=kxarray-1
    phase=dble(xsh)*dble(kxwrap)+dble(ysh)*dble(kywrap)
    pshft=dcmplx(cos(phase),sin(phase))
    frameout(kxarray,kyarray)=framein(kxarray,kyarray)*pshft
  end do 
end do
end subroutine Shiftframe

!---------------------------------------------------
subroutine Replaceextension_modname(filenamein,filenameout,extension,flag)
!---------------------------------------------------
! JLR 2/15
! replace the extension of a file with a new extension of arbitrary length
! add a 'flag' modifier to end of rootname of file
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

filenameout(1:lenroot)=filenamein(1:lenroot) ! transfer root of name
filenameout(lenroot+1:lenroot+lenflag+1)=trim(flag)
filenameout(lenroot+lenflag+1:lenroot+lenflag+1)="."
filenameout(lenroot+lenflag+2:lenroot+lenext+lenflag+1)=trim(extension)
filenameout=trim(filenameout)
return
end subroutine
!---------------------------------------------------
subroutine Changepath(filenamein,filenameout,path)
!---------------------------------------------------
! JLR 2/15
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

