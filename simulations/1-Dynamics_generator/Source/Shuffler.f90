		
module vars

integer,parameter:: SystemM=300


integer,parameter:: Ntimes=30000
integer(1),dimension(SystemM,SystemM):: element_state
integer(1),dimension(SystemM,SystemM):: element_stateSu
integer,dimension(SystemM,SystemM,2):: SuMatrix



 
end module


program main
use vars; 
implicit none
integer:: x1,y1,i,x2,y2,xsu,ysu
integer:: Nphoto 
integer::ts, tr

do x1=1,SystemM
do y1=1,SystemM
SuMatrix(x1,y1,1)=x1
SuMatrix(x1,y1,2)=y1
enddo
enddo


do x1=1,SystemM
do y1=1,SystemM
x2=1+int(rand()*SystemM)
y2=1+int(rand()*SystemM)
xsu=SuMatrix(x1,y1,1)
ysu=SuMatrix(x1,y1,2)

SuMatrix(x1,y1,1)=SuMatrix(x2,y2,1)
SuMatrix(x1,y1,2)=SuMatrix(x2,y2,2)


SuMatrix(x2,y2,1)=xsu
SuMatrix(x2,y2,2)=ysu

enddo
enddo






open(102,file="Input.dat", &
& form='unformatted', access='DIRECT', recl=90000)
open(103,file="Output.dat", access='DIRECT', recl=90000)


do Nphoto =1,Ntimes
if (nphoto==500*int(nphoto/500)) write(*,*) "Photo=",Nphoto

read(102,rec=Nphoto) element_state(:,:)


do x1=1,SystemM
do y1=1,SystemM
element_stateSu(x1,y1)=element_state(SuMatrix(x1,y1,1),SuMatrix(x1,y1,2))
enddo
enddo

write(103,rec=Nphoto)  element_stateSu(:,:)
enddo
call flush(103)
close(103)


end program
