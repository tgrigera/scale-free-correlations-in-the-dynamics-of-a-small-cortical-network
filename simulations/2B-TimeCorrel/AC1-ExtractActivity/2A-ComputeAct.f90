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


integer(1),dimension(SystemW,SystemW):: element_state


integer,parameter:: t_start=5000 !We used 5000 in the article
integer,parameter:: t_end=35000  !We used 100000 in the article
integer,parameter:: t_interval=1  !We used 100000 in the article
real(8):: activity
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

call system ('mkdir -p ResultsActOriginal')

open(103,file='../../1-Dynamics_generator/Original/ResultsState/StateRun' &
& // trim(adjustl(file_name))// '-Tr0.'  &
& // trim(adjustl(file_name2))// '.dat' , access='DIRECT', recl=90000)



open(203,file='ResultsActOriginal/ActRun' // trim(adjustl(file_name))// '-Tr0.'  &
& // trim(adjustl(file_name2))//'W'//trim(adjustl(file_name3))//'x'//trim(adjustl(file_name4))//'y'// &
& trim(adjustl(file_name5))// '.dat',status="unknown")


 
!Read The state
do photo_counter =2,int((t_end-t_start)/t_interval)
read(103,rec=photo_counter) element_state
call CalculateAct
write(203,*) photo_counter,activitY
enddo
close(103)




call flush(203)
close(203)

end program

!***************************************************
!SUBROUTINE CalculateCorrel
!***************************************************



subroutine CalculateAct
use variables; 
implicit none
integer x1,y1
activity=0
do x1=startx,startx+Window-1 ; do y1=starty,starty+Window-1
if (element_state(x1,y1)==1) then
activity=activity+1
!element_state(x1,y1)=1
!else
!element_state(x1,y1)=0
endif

enddo; enddo


end subroutine
