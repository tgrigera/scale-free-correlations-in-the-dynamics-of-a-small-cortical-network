module variables
!Run this code for differentvalues of Window,startx,starty
integer,parameter:: SystemW=300 !Size Over which I record activity (We used 500 on the article)
integer,parameter:: Window=100 !The actual window used, should change this value 

!Parameters For read File:
integer:: Run=1
real(8),parameter:: Threshold=0.3180
integer,parameter:: startx=1
integer,parameter:: starty=1

!Filenames
character(len=10) :: file_name
character(len=10) :: file_name2
character(len=10) :: file_name3
character(len=10) :: file_name4 
character(len=10) :: file_name5

integer,parameter:: maxdist=180 !Avoid Computing over distances larger than maxdist

integer(1),dimension(SystemW,SystemW):: element_state
real(8),dimension(0:2*maxdist):: Correl
integer(8),dimension(0:2*maxdist):: Events

integer,parameter:: t_start=5000 !We used 5000 in the article
integer,parameter:: t_end=35000  !We used 100000 in the article
integer,parameter:: t_interval=1  !We used 100000 in the article
end module


program main
use variables
implicit none
integer:: i
integer:: photo_counter 
real(8):: normalization !,mean2


write(file_name, '(i0)') Run
write(file_name2, '(i0)') int(Threshold*10000+0.5) !Add0.5 so that it is rounded properly
write(file_name3, '(i0)') Window
write(file_name4, '(i0)') startx
write(file_name5, '(i0)') starty

call system ('mkdir -p ResultsCorrelOriginal')

open(103,file='../../1-Dynamics_generator/Original/ResultsState/StateRun' &
& // trim(adjustl(file_name))// '-Tr0.'  &
& // trim(adjustl(file_name2))// '.dat' , access='DIRECT', recl=90000)





Events=0
Correl=0 
!Read The state
do photo_counter =1,int((t_end-t_start)/t_interval)
read(103,rec=photo_counter) element_state
call CalculateCorrel
enddo
close(103)


open(203,file='ResultsCorrelOriginal/CorrelRun' // trim(adjustl(file_name))// '-Tr0.'  &
& // trim(adjustl(file_name2))//'W'//trim(adjustl(file_name3))//'x'//trim(adjustl(file_name4))//'y'// &
& trim(adjustl(file_name5))// '.dat',status="unknown")

normalization=real(Correl(0))/real(Events(0))

do i=0, maxdist
if (Events(i)>0 ) then
write(203,*) i,real(Correl(i))/real(Events(i))/normalization,Correl(i), Events(i)
endif 
enddo
call flush(203)
close(203)

end program

!***************************************************
!SUBROUTINE CalculateCorrel
!***************************************************



subroutine CalculateCorrel
use variables; 
implicit none
integer:: x1,x2,y1,y2
integer dist

real(8):: activity_density
!integer:: Lmes=int(0.5*systemL)
real(8):: C1,C2

!Compute the total activity in this snapshot,setting state to 1 if neuron if active 
activity_density=0
do x1=startx,startx+Window-1 ; do y1=starty,starty+Window-1
if (element_state(x1,y1)==1) then
activity_density=activity_density+1.d0
else
element_state(x1,y1)=0
endif

enddo; enddo

if (activity_density<0.99) goto 333 !Skip Snapshots with no activity

activity_density=activity_density/real(Window*Window)
do x1=startx,startx+Window-1 ; do y1=starty,starty+Window-1
C1=real(element_state(x1,y1))-activity_density

do x2=startx,startx+Window-1 ; do  y2=starty,starty+Window-1 
dist=int(sqrt(0.5+(x2-x1)**2+(y2-y1)**2));
if (dist>maxdist) goto 222  
Events(dist)=Events(dist)+1
C2=real(element_state(x2,y2))-activity_density
Correl(dist)=Correl(dist)+C1*C2
222 continue
enddo; enddo;

enddo ;enddo





333 continue

end subroutine
