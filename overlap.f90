!  overlap.f90 
!
!  FUNCTIONS:
!	overlap      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: overlap
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program over

	implicit none
	! Variables
		real(8)::h(8000),p532(8000),p607(8000),x(8000),overlap(100,8000),air(8000),pi,sa,sm,p(8000),pfernald(100,8000)
	real(8)::betaraman(100,8000),betafernald(100,8000),alpharaman(100,8000),alpha(100,8000),betan2(8000),betam(8000)
	real(8)::k1,k2,k3,k4,k5,k6,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,fun1,fun2,slope
	real(8),external::lesat1,lesat2
	real(8)::x1(4),x2(3),y1(4),y2(3)
	real::gardient
	integer::i,j,k,l,ma,mb,mc,md,sr,o,q,qr
	pi=3.1415926
	sm=8.0*pi/3.0
	sum7=0.0
	sum8=0.0
	fun1=0.0
	fun2=9999999.0
	open(unit=11,file="lanzhou_air.txt")
	open(unit=12,file="overlap.txt")
	open(unit=13,file="lanzhou.txt")
	open(unit=14,file="lanzhou_raman.txt")
	open(unit=15,file="lanzhou_fernald.txt")
	open(unit=16,file="xinhao.txt")
	OPEN(unit=17,file="diff.txt")
	open(unit=18,file="slope_k.txt")
	do i=1,1000
	read(11,*)air(i)
	air(i)=air(i)
    if(i>900)then
	overlap(:,i)=1.0
	end if
	end do
	betan2=0.78*air/sm
	betam=air/sm
	do j=1,1200
	read(13,*)h(j),p532(j),p607(j),p(j)
	h(j)=h(j)/1000.0
	end do
	x=p
	call qidian(x)
	p=x
	x=p
	call qidian(x)
	p=x
	x=p
	call qidian(x)
	p=x
	x=p532
	call qidian(x)
	p532=x
	x=p532
	call qidian(x)
	p532=x
	do k=1,100
	sa=0.0+real(k)
	betaraman(k,901)=betam(901)*sm/sa
	betafernald(k,901)=betam(901)*sm/sa
	do l=900,30,-1
	x2=h(l:l+2)
	y2=betam(l:l+2)
	slope=lesat2(x2,y2)
	sum1=slope/betam(l+1)
	sum6=slope
	sum2=(1.0-532.0**4/607.0**4)*betam(l+1)*sm
	sum3=(1.0-532.0/607.0)*betaraman(k,l+1)*sa
	x2=h(l:l+2)
	y2=p(l:l+2)
	slope=lesat2(x2,y2)
	sum4=slope/p(l+1)
	sum5=betam(l+1)+betaraman(k,l+1)
	k1=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.00375
	x1=h(l-1:l+2)
	y1=betam(l-1:l+2)
	slope=lesat1(x1,y1)
	sum1=slope/(betam(l+1)+betam(l))*2.0
	sum6=slope
	sum2=(1.0-532.0**4/607.0**4)*(betam(l+1)+betam(l))/2.0*sm
	sum3=(1.0-532.0/607.0)*(betaraman(k,l+1)-k1/2.0)*sa
	x1=h(l-1:l+2)
	y1=p(l-1:l+2)
	slope=lesat1(x1,y1)
	sum4=slope/(p(l+1)+p(l))*2.0
	sum5=(betam(l+1)+betam(l))/2.0+betaraman(k,l+1)-k1/2.0
	k2=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.00375
    x1=h(l-1:l+2)
	y1=betam(l-1:l+2)
	slope=lesat1(x1,y1)
	sum1=slope/(betam(l+1)+betam(l))*2.0
	sum6=slope
	sum2=(1.0-532.0**4/607.0**4)*(betam(l+1)+betam(l))/2.0*sm
	sum3=(1.0-532.0/607.0)*(betaraman(k,l+1)-k2/2.0)*sa
	x1=h(l-1:l+2)
	y1=p(l-1:l+2)
	slope=lesat1(x1,y1)
	sum4=slope/(p(l+1)+p(l))*2.0
	sum5=(betam(l+1)+betam(l))/2.0+betaraman(k,l+1)-k2/2.0
	k3=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.00375

	x2=h(l-1:l+1)
	y2=betam(l-1:l+1)
	slope=lesat2(x2,y2)
	sum1=slope/betam(l)
	sum6=slope
	sum2=(1.0-532.0**4/607.0**4)*betam(l)*sm
	sum3=(1.0-532.0/607.0)*(betaraman(k,l+1)-k3)*sa
	x2=h(l-1:l+1)
	y2=p(l-1:l+1)
	slope=lesat2(x2,y2)
	sum4=slope/p(l)
	sum5=betam(l)+betaraman(k,l+1)-k3
	k4=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.00375
    betaraman(k,l)=betaraman(k,l+1)-(k1+2.0*k2+2.0*k3+k4)/6.0
	end do
	do ma=30,900
do mb=ma+1,900
sum7=sum7+2.0*(1.0+532.0**4/607.0**4)*betam(mb)*sm
sum8=sum8+2.0*(1.0+532.0/607.0)*betaraman(k,mb)*sa
end do
sum7=(sum7+(1.0+532.0**4/607.0**4)*(betam(ma)+betam(901))*sm)*0.00375/2.0
sum8=(sum8+(1.0+532.0/607.0)*(betaraman(k,ma)+betaraman(k,901))*sa)*0.00375/2.0
overlap(k,ma)=p607(ma)*h(ma)**2*betam(901)/(p607(901)*h(901)**2*betam(ma))*exp(-1.0*sum7)*exp(-1.0*sum8)
sum7=0.0
sum8=0.0
end do
x=overlap(k,:)
call qidian(x)
overlap(k,:)=x
do qr=2,901
pfernald(k,qr)=h(qr)**2*p532(qr)/overlap(k,qr)
!pfernald(k,qr)=h(qr)**2*p532(qr)
end do
do mc=900,2,-1
sum9=(sa-sm)*(betam(mc)+betam(mc+1))*0.00375
sum10=pfernald(k,mc+1)/(betafernald(k,mc+1)+betam(mc+1))
sum11=sa*(pfernald(k,mc+1)+pfernald(k,mc)*exp(sum9))*0.00375
betafernald(k,mc)=pfernald(k,mc)*exp(sum9)/(sum10+sum11)-betam(mc)
end do
do md=201,800
fun1=fun1+abs(betaraman(k,md)-betafernald(k,md))**2
end do
write(17,*)k,fun1
if (fun1<=fun2)then
fun2=fun1
fun1=0.0
sr=sa
else
fun1=0.0
fun2=fun2
sr=sr
end if
end do
write(*,*)sr
do o=1,100
do q=2,900
write(12,*)o,h(q),overlap(o,q)
write(14,*)o,h(q),betaraman(o,q)
write(15,*)o,h(q),betafernald(o,q)
write(16,*)o,h(q),p(q)
end do
end do

	end program 


subroutine qidian(x)
implicit none
real(8)::ctm(8000),x(8000)
integer::i,j
real(8)::sum1,sum2,sum3
do i=1,7000
if(i==1)then
ctm(i)=(13.0*x(1)+10.0*x(2)+7.0*x(3)+4.0*x(4)+x(5)-2.0*x(6)-5.0*x(7))/28.0
else if(i==2)then
ctm(i)=(5.0*x(1)+4.0*x(2)+3.0*x(3)+2.0*x(4)+x(5)-x(7))/14.0
else if(i==3)then
ctm(i)=(7.0*x(1)+6.0*x(2)+5.0*x(3)+4.0*x(4)+3.0*x(5)+2.0*x(6)+x(7))/14.0
else if(i>3)then
ctm(i)=(x(i-3)+x(i-2)+x(i-1)+x(i)+x(i+1)+x(i+2)+x(i+3))/7.0
end if
end do
do j=1,7000
x(j)=ctm(j)
end do
return
end subroutine

function lesat2(x2,y2)
real(8)::x2(3),y2(3)
real(8)::sum1,sum2,sum3,sum4,b
real(8)::lesat2
integer::i,j
sum1=0.0
sum2=0.0
sum3=0.0
sum4=0.0
do i=1,3
sum1=sum1+x2(i)
sum2=sum2+y2(i)
sum3=sum3+x2(i)*y2(i)
sum4=sum4+x2(i)**2
end do
lesat2=(3.0*sum3-sum1*sum2)/(3.0*sum4-sum1**2)
b=sum2+lesat2*sum1
!do j=1,3
!write(*,*)b-lesat2*x2(j)
!end do
!write(*,*)b
sum1=0.0
sum2=0.0
sum3=0.0
sum4=0.0
return
end function

function lesat1(x1,y1)
real(8)::x1(4),y1(4)
real(8)::sum1,sum2,sum3,sum4
real(8)::lesat1
integer::i,j
sum1=0.0
sum2=0.0
sum3=0.0
sum4=0.0
do i=1,4
sum1=sum1+x1(i)
sum2=sum2+y1(i)
sum3=sum3+x1(i)*y1(i)
sum4=sum4+x1(i)**2
end do
sum1=sum1/4.0
sum2=sum2/4.0
sum3=sum3/4.0
sum4=sum4/4.0
lesat1=(sum3-sum1*sum2)/(sum4-sum1**2)
sum1=0.0
sum2=0.0
sum3=0.0
sum4=0.0
return
end function


	
