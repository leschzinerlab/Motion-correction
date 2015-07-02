!———————————————————————————
! alignparts_starfilehandler                         Dated 31/3/15
!-------------------------------------------------------
! JLR 3/14


PROGRAM Alignparts_starfilehandler
implicit none
integer, parameter                :: dp=kind(0.d0) 
! I/O
character(150)                    :: instarfile,outstarfile,inmovie,inmic,outstk,outcoord,tempname,movielist,coordlist
character(150)                    :: moviepath,particlepath
character(150)                    :: prevmic
character(20)                     :: stkext,stkmod,movext,movmod,coordext,coordmod
! Starfile variables
character(50), dimension (50)     :: star_colcontent ! First position is column number, second is position in name
integer                           :: star_coltot,star_colpos,star_col
integer                           :: star_io
character(500)                    :: star_record,star_recordmod
character(5000)                   :: star_line
integer                           :: star_coordx,star_coordy,star_micname,star_imagename
! Logical
logical                           :: firstmic
! Image info
real(dp)                          :: xposreal,yposreal,density
integer                           :: nx,ny,boxsize,xpos,ypos,xposs,xposf,yposs,yposf
integer                           :: nparts,npartstot,nskipped,nmovies



read (5,'(a)') instarfile
read (5,'(a)') outstarfile
read (5,'(a)') movielist
read (5,'(a)') coordlist
read (5,*) boxsize,nx,ny
read (5,'(a)') moviepath
read (5,'(a)') particlepath
read (5,*) movmod
read (5,*) movext
read (5,*) stkmod
read (5,*) stkext

!-----------------------------------
! Open input and output .star files
!-----------------------------------
open(unit=1,file=instarfile,status="old")
open(unit=2,file=outstarfile,status="replace")
open(unit=3,file=movielist,status="replace")
open(unit=4,file=coordlist,status="replace")

!--------------------------------------------------------------
! File streams
! Stream 1 : In starfile
! Stream 2 : Out starfile
! Stream 3 : Out movielist  
! Stream 4 : Out coordlist
! Stream 7 ; Out coordfile (outcoord)
!--------------------------------------------------------------

!----------------------------------------------------------------------------------------
! Learn the columns of the input .star file and write the header of the output .star file
!----------------------------------------------------------------------------------------
write (*,'(a)') "Reading header from input .star file and writing header to output .star file"
call Star_learncolumns(1,star_colcontent,star_coltot)
call Star_writeheader(2,star_colcontent,star_coltot)
! Define which column contains micrographname, x-position, and y-position of particles
star_coordx=Star_colpos(star_colcontent,'_rlnCoordinateX')
star_coordy=Star_colpos(star_colcontent,'_rlnCoordinateY')
star_micname=Star_colpos(star_colcontent,'_rlnMicrographName')
star_imagename=Star_colpos(star_colcontent,'_rlnImageName')

!----------------------------------------------------------------------------------------
! Loop over .star file lines
!----------------------------------------------------------------------------------------
coordext="coord"
coordmod=""
density=0
nparts=0    ! Number of particles processed per film
npartstot=0 ! Number of particles transferred to new .star file
nskipped=0  ! Number of particles not transferred to new .star file
nmovies=0   ! Number of different movies observed in .star file
prevmic=""
firstmic=.true.
Eachstarline: do 
  read (1,'(a)',iostat=star_io) star_line ! Read a line from the starfile
  Starlineexists: if (star_io>0) then 
    write (*,*) "problem with starfile"
    exit
  else if (star_io<0) THEN
    write (*,*) "End of starfile reached"
    exit
  end if Starlineexists
  call Star_readrecord(star_line,star_coltot,star_micname,star_record) ! Read name of micrograph
  inmic=trim(star_record)
  if (trim(inmic).eq."") exit 
  ! Get x-coord and y-coord from star_line
  call Star_readrecord(star_line,star_coltot,star_coordx,star_record)
  read(star_record,*) xposreal
  xpos=int(xposreal) ! for some reason Relion writes down coordinates as real number
  call Star_readrecord(star_line,star_coltot,star_coordy,star_record)
  read(star_record,*) yposreal
  ypos=int(yposreal)
  xposs=xpos-boxsize/2
  xposf=xpos+boxsize/2-1
  yposs=ypos-boxsize/2
  yposf=ypos+boxsize/2-1
  Goodparticle: if (xposs.lt.0.or.xposf.gt.nx.or.yposs.lt.0.or.yposf.gt.ny) then
    write (*,*) "Coordinate out of bounds for specified boxsize at record:",nparts
    write (*,*) "This problem can cause a major mismatch between the particles in a stack and star file"
    write (*,*) "skipping the offending particle position and excluding it from the output star file"
    nskipped=nskipped+1
  else
    Newmicrograph: if (trim(inmic).ne.trim(prevmic)) then
      if (.not.firstmic) then
        close(7) ! Close old coordinate file
      end if
      nparts=0
      nmovies=nmovies+1
      firstmic=.false.
      prevmic=trim(inmic)
      ! Generate necessary filenames
      call Replaceextension_modname(inmic,tempname,movext,movmod) ! Generate name of input movie
      call Changepath(tempname,inmovie,moviepath)                       ! Adjust path of input movie
      call Replaceextension_modname(inmic,tempname,stkext,stkmod) ! Generate output stack name
      call Changepath(tempname,outstk,particlepath)                     ! Adjust path of output stack file
      call Replaceextension_modname(inmic,tempname,coordext,coordmod)    ! Generate output coordfile name
      call Changepath(tempname,outcoord,moviepath)                      ! Adjust path of output cordial     
      open(unit=7,file=outcoord,status="replace")
      write (7,'(A32)') "        x         y      density"! Write header for coordinate file
      write (3,*) trim(inmovie)
      write (4,*) trim(outcoord)
    end if Newmicrograph
    nparts=nparts+1;npartstot=npartstot+1
    ! write a line in the output coord file
    write (7,'(2I10,F12.1)') xpos,ypos,density ! write a line in the outcoord
    write (star_recordmod,'(i6.6,a1,a)') nparts,'@',outstk
    do star_col=1,star_coltot
       call Star_readrecord(star_line,star_coltot,star_col,star_record)
       if (star_col.ne.star_imagename) then
         write (2,'(a,a)',advance='no') trim(star_record),'  '
       else if  (star_col.eq.star_imagename) then
         write (2,'(a,a)',advance='no') trim(star_recordmod),'  '
       end if
    end do
    write (2,'(a)',advance='yes') ' '
  end if Goodparticle
end do Eachstarline
!

write (*,'(a)')    "**************************************"
write (*,'(a,i8)') "Number of .star file lines transferred:",npartstot
write (*,'(a,i8)') "Number of outlying particles removed:",nskipped
write (*,'(a,i8)') "Number of movies processed that contained usable particles:",nmovies
write (*,'(a)')    "**************************************"
! Close off :instarfile,outstarfile,movielist,coordlist,coordfile?
close(1) ! input starfile
close(2) ! output starfile
close(3) ! output movielist
close(4) ! output coordlist
close(7) ! output outcoord

END PROGRAM Alignparts_starfilehandler

!====================================================
integer function Star_colpos(star_colcontent,query)
!====================================================
! JLR 2/15
! Returns the position in columncontents of the query string
! e.g. Star_colpos(star_colcontent,'_rlnMagnification')
implicit none 
character(50), dimension (50), intent(IN) :: star_colcontent
character(*),intent(IN)                   :: query
integer,dimension(1)                      :: temp
temp=maxloc(index(star_colcontent,query))
star_colpos=temp(1)
return
end function

!====================================================
subroutine Star_learncolumns(stream,star_colcontent,star_coltot)
!====================================================
! JLR 2/15
! Returns the array columncontents and the integer star_coltot
implicit none
integer, intent(IN)                        :: stream
integer                                    :: i!,star_colpos
integer, intent(OUT)                       :: star_coltot ! total number of columns
character(200)                             :: headerline
character(50), dimension (50), intent(OUT) :: star_colcontent ! First position is column number, second is position in name

! Skip header
do
  read(stream,*) headerline
  if (index(headerline,'loop_').ne.0) then
    exit
  end if
end do
! Learn contents from header and learn total number of columns per line
i=0
do
  read(stream,*) headerline  ! read a line from the header
  if (index(headerline,'_').ne.1) then ! Determine if it has a column label
    backspace(stream)
    exit
  else
    i=i+1
    star_colcontent(i) = trim(headerline)
  end if
end do
star_coltot=i
return
end subroutine Star_learncolumns

!====================================================
subroutine Star_readrecord(star_line,star_coltot,star_col,star_record)
!====================================================
! JLR 2/15
! Given a line from a star file, retrieve the information at the specified column
implicit none
!character(50), dimension (50), intent(IN) :: star_colcontent
character(5000), intent(IN)               :: star_line
integer, intent(IN)                       :: star_coltot,star_col
character(500), intent(out)               :: star_record
integer                                   :: i,lowest,highest

highest=0;lowest=0
if (star_col.gt.star_coltot) then
  write (*,*) "Requested record is in a column that does not exist"
  stop
else if (star_col.eq.star_coltot) then ! if this is the last record, set its position
  highest=len_trim(star_line)
  lowest=index(star_line(:highest),' ',back=.true.)+1
else ! otherwise find its position
  highest=0
  lowest=1
  do i=1,star_col
    lowest=highest+1
    ! set lowest and highest to the col values of the blank spaces that surround text in the line
    do 
      highest = index(star_line(lowest:),' ')+lowest-1 ! Find the next nonblank
      if (lowest.ge.len(star_line)) then
        exit
      else if (highest-lowest.le.1) then ! two columns were adjacent, try again
        lowest=highest+1
      else ! two columns were not adjacent, found it
        exit
      end if
    end do
  end do
end if
read(star_line(lowest:highest),'(a)') star_record ! read the relevant record

return
end subroutine Star_readrecord

!====================================================
subroutine Star_writeheader(stream,star_colcontent,star_coltot)
!====================================================
! JLR 2/15
implicit none
integer, intent(IN)                        :: stream
integer                                    :: i!,star_colpos
integer, intent(IN)                        :: star_coltot ! total number of columns
character(50), dimension (50), intent(IN)  :: star_colcontent 

write (stream,'(a)') ''
write (stream,'(a5)') 'data_'
write (stream,'(a)') ''
write (stream,'(a5)') 'loop_'

do i=1,star_coltot
  write (stream,'(a,a2,i0)') trim(star_colcontent(i)),' #',i
end do
return
end subroutine Star_writeheader

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
