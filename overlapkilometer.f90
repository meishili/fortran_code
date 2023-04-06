!  overlapkilometer.f90 
!
!  FUNCTIONS:
!	overlapkilometer      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: overlapkilometer
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
module pinghua
implicit none
integer::N
contains
subroutine qidian(x)
implicit none
real(8)::ctm(4000),x(4000)
integer::i,j
real(8)::sum1,sum2,sum3
do i=1,3000
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
do j=1,3000
x(j)=ctm(j)
end do
return
end subroutine
end module


	program overlapkilometer
	use pinghua

	implicit none
	real(8)::h(4000),p532(4000),p607(4000),x(4000),overlap(100,4000),air(4000),pi,sa,sm,p(4000),pfernald(100,4000)
	real(8)::betaraman(100,4000),betafernald(100,4000),alpharaman(100,4000),alpha(100,4000),betan2(4000),betam(4000)
	real(8)::k1,k2,k3,k4,k5,k6,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,fun1,fun2
	integer::i,j,k,l,ma,mb,mc,md,sr,o,q,qr
	pi=3.1415926
	sm=8.0*pi/3.0
	sum7=0.0
	sum8=0.0
	fun1=0.0
	fun2=9999999.0
	open(unit=11,file="tazhongair.txt")
	open(unit=12,file="overlap.txt")
	open(unit=13,file="output.txt")
	open(unit=14,file="raman.txt")
	open(unit=15,file="fernald1.txt")
	open(unit=16,file="tazhong.txt")
	do i=1,4000
	read(11,*)air(i)
	air(i)=air(i)*1000.0
    if(i>600)then
	overlap(:,i)=1.0
	end if
	end do
	betan2=0.78*air/sm
	betam=air/sm
	do j=1,4000
	read(13,*)h(j),p532(j),p607(j)
	h(j)=h(j)/1000.0
	end do
	x=p607
	call qidian(x)
	p607=x
	p=p607/p532	
	x=p
	call qidian(x)
	p=x
	betaraman(:,601)=betam(601)
	betafernald(:,601)=betam(601)
	do k=1,100
	sa=0.0+real(k)
	do l=600,2,-1
	sum1=(betam(l+2)-betam(l))/0.015/betam(l+1)
	sum2=(1.0-607.0**4/532.0**4)*betam(l+1)*sm
	sum3=(1.0-607.0/532.0)*betaraman(k,l+1)*sa
	sum4=(p(l+2)-p(l))/0.015/p(l+1)
	sum5=betam(l+1)+betaraman(k,l+1)
	sum6=(betam(l+2)-betam(l))/0.015
	k1=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.0075
	sum1=(betam(l+1)-betam(l))/0.0075/(betam(l+1)+betam(l))*2.0
	sum2=(1.0-607.0**4/532.0**4)*(betam(l+1)+betam(l))/2.0*sm
	sum3=(1.0-607.0/532.0)*(betaraman(k,l+1)-k1/2.0)*sa
	sum4=(p(l+1)-p(l))/0.0075/(p(l+1)+p(l))*2.0
	sum5=(betam(l+1)+betam(l))/2.0+betaraman(k,l+1)-k1/2.0
	sum6=(betam(l+1)-betam(l))/0.0075
	k2=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.0075
	sum1=(betam(l+1)-betam(l))/0.0075/(betam(l+1)+betam(l))*2.0
	sum2=(1.0-607.0**4/532.0**4)*(betam(l+1)+betam(l))/2.0*sm
	sum3=(1.0-607.0/532.0)*(betaraman(k,l+1)-k2/2.0)*sa
	sum4=(p(l+1)-p(l))/0.0075/(p(l+1)+p(l))*2.0
	sum5=(betam(l+1)+betam(l))/2.0+betaraman(k,l+1)-k2/2.0
	sum6=(betam(l+1)-betam(l))/0.0075
	k3=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.0075
	sum1=(betam(l+1)-betam(l-1))/0.015/betam(l)
	sum2=(1.0-607.0**4/532.0**4)*betam(l)*sm
	sum3=(1.0-607.0/532.0)*(betaraman(k,l+1)-k3)*sa
	sum4=(p(l+1)-p(l-1))/0.015/p(l)
	sum5=betam(l)+betaraman(k,l+1)-k3
	sum6=(betam(l+1)-betam(l-1))/0.015
	k4=((sum1+sum2+sum3-sum4)*sum5-sum6)*0.0075
    betaraman(k,l)=betaraman(k,l+1)-(k1+2.0*k2+2.0*k3+k4)/6.0
	end do
	do ma=2,600
if(ma<600)then
do mb=ma+1,600
sum7=sum7+2.0*(1.0+607.0**4/532.0**4)*betam(mb)*sm
sum8=sum8+2.0*(1.0+607.0/532.0)*betaraman(k,mb)*sa
end do
sum7=(sum7+(1.0+607.0**4/532.0**4)*(betam(ma)+betam(601))*sm)*0.0075/2.0
sum8=(sum8+(1.0+607.0/532.0)*(betaraman(k,ma)+betaraman(k,601))*sa)*0.0075/2.0
else if(ma==600)then
sum7=(betam(600)+betam(601))/2.0*0.0075*(1.0+607.0**4/532.0**4)*sm
sum8=(betaraman(k,600)+betaraman(k,601))/2.0*0.0075*(1.0+607.0/532.0)*sa
end if
overlap(k,ma)=p607(ma)*h(ma)**2*betam(601)/(p607(600)*h(601)**2*betam(ma))*exp(-1.0*sum7)*exp(-1.0*sum8)
sum7=0.0
sum8=0.0
end do
x=overlap(k,:)
call qidian(x)
overlap(k,:)=x
do qr=2,4000
pfernald(k,qr)=h(qr)**2*p532(qr)/overlap(k,qr)
end do
do mc=600,2,-1
sum9=(sa-sm)*(betam(mc)+betam(mc+1))*0.0075
sum10=pfernald(k,mc+1)/(betafernald(k,mc+1)+betam(mc+1))
sum11=sa*(pfernald(k,mc+1)+pfernald(k,mc)*exp(sum9))*0.0075
betafernald(k,mc)=pfernald(k,mc)*exp(sum9)/(sum10+sum11)-betam(mc)
end do
do md=1,400
fun1=fun1+abs(betaraman(k,md)-betafernald(k,md))
end do
if (fun1<=fun2)then
fun2=fun1
fun1=0.0
sr=sa
else
fun1=0.0
fun2=fun2
end if
end do
write(*,*)sr
do o=1,100
do q=2,600
write(12,*)o,h(q),overlap(o,q)
write(14,*)o,h(q),betaraman(o,q)
write(15,*)o,h(q),betafernald(o,q)
end do
end do
	end program overlapkilometer

