
! CONVENTIONS AND CHANGES: 
! 1) All the periods have been changed by underscores
! 2) Instead using a derive type in fortran, fields of variable "vgm" has been split
!    into two new arrays: "vgm_ibin" and "vgm_ipair"
! 3) The integers n,p,n_out has been introduced to define the array sizes
!	 B= size of gamma_hat
!        n= size of vgm_pair
!        p= size of rho_grid
! 4) The output of the original R-function is included as the last argument in the subrutine
!------------------- hg function ---------------------------------------

subroutine hg(rho,m,value)

 implicit none
 integer,parameter:: dp = KIND(1.0d0)
 integer,intent(in)::m
 real(kind=dp),intent(in)::rho(m)
 real(kind=dp)::value(m),a,b,cc,fn(m),fnold(m),facn
 integer::n
 real(kind=dp),external::fgamma

 ! Note: the Gamma function was included as intrinsic in Gfortran 4.3

   a     = 3.0_dp/4.0_dp
   b     = 3.0_dp/4.0_dp
   cc    = 1.0_dp/2.0_dp
   fn    = fgamma(a)*fgamma(b)/fgamma(cc)
   facn  = 1.0_dp
   n     = 1
   fnold = 0.1_dp

   do while (maxval(abs(fn-fnold)/fnold) > 0.0001_dp)
      facn  = facn*n
      fnold = fn
      fn    = fn + fgamma(a + n)*fgamma(b + n)*rho**n/(fgamma(cc+n)*facn)
      n     = n + 1
    end do

   value=fn*fgamma(cc)/(fgamma(a)*fgamma(b))

end subroutine hg

!------------------- cor_sqrtabs function ------------------------------
subroutine cor_sqrtabs(rho,m,value) 

 implicit none

 interface 
   subroutine hg(rho,m,value)
      integer,parameter:: dp = KIND(1.0d0)
      integer,intent(in)::m
      real(kind=dp),intent(in)::rho(m)
      real(kind=dp)::value(m)
   end subroutine hg
 end interface

 integer,parameter:: dp = KIND(1.0d0)
 integer,intent(in)::m
 real(kind=dp),intent(in)::rho(m)
 real(kind=dp)::pi,value(m),val(m)
 real(kind=dp),external::fgamma
 pi=acos(-1.0_dp)


 ! Note: the Gamma function was included as intrinsi! in Gfortran 4.3
 call hg(rho**2,m,val)
 value=fgamma(0.75_dp)**2*((1.0_dp-rho**2)*val-1.0_dp)/(sqrt(pi)-fgamma(0.75_dp)**2)


end subroutine cor_sqrtabs

!------------------- approx_linear function ------------------------------
subroutine approx_linear(x,y,n,v,m,yleft,yright,value) 

 implicit none
 integer,parameter:: dp = KIND(1.0d0)
 integer,intent(in)::n,m
 integer::j,l
 real(kind=dp)::value(m)
 real(kind=dp),intent(in)::x(n),y(n),v(m),yleft,yright
 
 ! Loop in the components of v
 do j=1,m

    if(v(j)<=x(1)) then
       value(j)=yleft
    elseif(v(j)>=x(n)) then
       value(j)=yright
    else
       l=count(v(j)>x)
       value(j)=y(l)+(y(l+1)-y(l))*((v(j)-x(l))/(x(l+1)-x(l)))
    endif

 enddo
end subroutine approx_linear


subroutine cov_bin_fun(B,n,p,i,j,vgm_ibin,vgm_ipair,gamma_hat,mean_cv)

 implicit none

 interface 
   subroutine approx_linear(x,y,n,v,m,yleft,yright,value) 
      integer,parameter:: dp = KIND(1.0d0)
      integer,intent(in)::n,m
      real(kind=dp)::value(m)
      real(kind=dp),intent(in)::x(n),y(n),v(m),yleft,yright
   end subroutine approx_linear
 end interface

 interface 
   subroutine cor_sqrtabs(rho,m,value)
     integer,parameter:: dp = KIND(1.0d0)
     integer,intent(in)::m
     real(kind=dp),intent(in)::rho(m)
     real(kind=dp)::value(m)
   end subroutine cor_sqrtabs
 end interface

! Simple precision
!  integer,parameter:: dp = KIND(1.0)
! Double precision
  integer,parameter:: dp = KIND(1.0d0)

 integer,intent(in)::B,n,p,i,j
 integer,intent(in):: vgm_ipair(n,2),vgm_ibin(n)
 real(kind=dp),intent(in):: gamma_hat(B)
 real(kind=dp),intent(inout):: mean_cv

 integer::iter,k,l,m,m1,m2,m3,m4,n_bin_i,n_bin_j,nmax,ib,jb,n_out,unique_vgm_ibin(B)
 integer,dimension(:),allocatable::vrow1,vrow2,aux, &
                      vind1,vind2,vind3,vind4,ind,ind1,ind2,ind3,ind4
 integer,dimension(:,:),allocatable::aux1,aux2,mp,mp1,mp2, &
                      mpair1(:,:),mpair2(:,:),mpair3(:,:),mpair4(:,:)

 real(kind=dp):: rho_grid(p+1),f_grid(p+1)
 real(kind=dp),allocatable::cr(:),cv(:),tmp(:),mpair1_hat(:),mpair2_hat(:), &
                                    mpair3_hat(:),mpair4_hat(:),auxcv(:)

 logical,allocatable::ind_bin_i(:),ind_bin_j(:), &
                      mask(:),mask1(:),mask2(:),mask3(:),mask4(:)

 ! Memory allocation for ind_bin_i,ind_bin_j
 allocate(ind_bin_i(n),ind_bin_j(n))

! ------------ rho evaluation ------------------------------------------
 k=0; rho_grid(1:p) = (/((0.96_dp/(p-1))*(k-1),k=1,p)/)
 call cor_sqrtabs(rho_grid(1:p),p,f_grid(1:p))
 rho_grid(p+1) = 1.0_dp
 f_grid(p+1) = 1.0_dp

! ------------- Sort computations --------------------------------------
 ind_bin_i = (vgm_ibin == i)
 ind_bin_j = (vgm_ibin == j)
 n_bin_i = count(ind_bin_i)
 n_bin_j = count(ind_bin_j)

 ! Memory allocation for ind_bin_i,ind_bin_j
 n_out=n_bin_i*n_bin_j
 allocate(vrow1(n_out),vrow2(n_out))

 ! Computing the values for vrow1,vrow2
 do k=1,n_bin_j
    vrow1(1+(k-1)*n_bin_i:k*n_bin_i) = (/(l,l=1,n_bin_i)/)
 enddo
 do k=1,n_bin_j
    vrow2(1+(k-1)*n_bin_i:k*n_bin_i) = k
 enddo

 ! Memory allocation for mp1,mp2 and mp
 allocate(aux1(n_bin_i,2),aux2(n_bin_j,2), &
          mp1(n_out,2),mp2(n_out,2))

 ! Construction mp1
 aux1=vgm_ipair(pack((/(l,l=1,n)/),ind_bin_i),:)
 mp1=aux1(vrow1,:)

 ! Construction mp2
 aux2=vgm_ipair(pack((/(l,l=1,n)/),ind_bin_j),:)
 mp2=aux2(vrow2,:)

 ! Construction mp
 allocate(mp(n_out,4))
 mp(:,1:2)=mp1
 mp(:,3:4)=mp2

 ! Memory deallocation
 deallocate(aux1,aux2,vrow1,vrow2,mp1,mp2,ind_bin_i,ind_bin_j)

!---------------- mpair* matrices: --------------------------------------
! 1: pair(1)
! 2: pair(2)
! 3: pair(1)-pair(2)
! 4: index
! 5: gamma_hat, if different pairs --> go to mpair*_hat(:)

 ! Memory allocation
 allocate(mpair1(n_out,4),mpair2(n_out,4), &
          mpair3(n_out,4),mpair4(n_out,4))

 ! Initialize to zero
 mpair1 = 0; mpair2 = 0; mpair3 = 0; mpair4 = 0

 ! Compute the values for mpair*(:,1:3)
 mpair1(:,1)= merge(mp(:,2),mp(:,3),mp(:,2)<=mp(:,3))
 mpair1(:,2)= merge(mp(:,2),mp(:,3),mp(:,2)>=mp(:,3))
 mpair2(:,1)= merge(mp(:,1),mp(:,4),mp(:,1)<=mp(:,4))
 mpair2(:,2)= merge(mp(:,1),mp(:,4),mp(:,1)>=mp(:,4))
 mpair3(:,1)= merge(mp(:,1),mp(:,3),mp(:,1)<=mp(:,3))
 mpair3(:,2)= merge(mp(:,1),mp(:,3),mp(:,1)>=mp(:,3))
 mpair4(:,1)= merge(mp(:,2),mp(:,4),mp(:,2)<=mp(:,4))
 mpair4(:,2)= merge(mp(:,2),mp(:,4),mp(:,2)>=mp(:,4))

 mpair1(:,3)= mpair1(:,2)-mpair1(:,1)
 mpair2(:,3)= mpair2(:,2)-mpair2(:,1)
 mpair3(:,3)= mpair3(:,2)-mpair3(:,1)
 mpair4(:,3)= mpair4(:,2)-mpair4(:,1)

 ! Memory allocation for temporal arrays
 allocate(mask1(n_out),mask2(n_out),aux(n_out),mask3(n_out),mask4(n_out))
	
 ! Compute the mask for the mpair* (the pairs outside the diagonal)
 mask1=(mpair1(:,3) /= 0)
 mask2=(mpair2(:,3) /= 0)
 mask3=(mpair3(:,3) /= 0)
 mask4=(mpair4(:,3) /= 0)

 m1=count(mask1)
 m2=count(mask2)
 m3=count(mask3)
 m4=count(mask4)

 ! Memory allocation for vind*
 allocate(vind1(m1),vind2(m2),vind3(m3),vind4(m4))

 aux=(/(k,k=1,n_out)/)
 vind1= pack(aux,mask1)
 vind2= pack(aux,mask2)
 vind3= pack(aux,mask3)
 vind4= pack(aux,mask4)

 ! Memory deallocation
 deallocate(mask1,mask2,mask3,mask4)

 nmax= maxval(vgm_ipair)
 mpair1(vind1,4) = (mpair1(vind1,1)-1)*nmax- &
                   ((mpair1(vind1,1)-1)*mpair1(vind1,1))/2 + &
                    mpair1(vind1,2)-mpair1(vind1,1)

 mpair2(vind2,4) = (mpair2(vind2,1)-1)*nmax- &
                   ((mpair2(vind2,1)-1)*mpair2(vind2,1))/2 + &
                    mpair2(vind2,2)-mpair2(vind2,1)

 mpair3(vind3,4) = (mpair3(vind3,1)-1)*nmax- &
                   ((mpair3(vind3,1)-1)*mpair3(vind3,1))/2 + &
                    mpair3(vind3,2)-mpair3(vind3,1)

 mpair4(vind4,4) = (mpair4(vind4,1)-1)*nmax- &
                   ((mpair4(vind4,1)-1)*mpair4(vind4,1))/2 + &
                    mpair4(vind4,2)-mpair4(vind4,1)

 ! For gamma_hat at each bin
 allocate(ind1(m1),ind2(m2),ind3(m3),ind4(m4))

 ! Construction of the uniquenness for vgm_ibin
 unique_vgm_ibin=maxval(vgm_ibin)+1
 iter=1
 do k=1,int(maxval(vgm_ibin))
     if(any(vgm_ibin==k)) then
        unique_vgm_ibin(iter)=k
        iter=iter+1
     endif
 enddo

 do k=1,m1
   ind1(k)=count(unique_vgm_ibin<vgm_ibin(mpair1(vind1(k),4)))+1
 enddo
 do k=1,m2
   ind2(k)=count(unique_vgm_ibin<vgm_ibin(mpair2(vind2(k),4)))+1
 enddo
 do k=1,m3
   ind3(k)=count(unique_vgm_ibin<vgm_ibin(mpair3(vind3(k),4)))+1
 enddo
 do k=1,m4
   ind4(k)=count(unique_vgm_ibin<vgm_ibin(mpair4(vind4(k),4)))+1
 enddo

 ib = count(unique_vgm_ibin<i)+1
 jb = count(unique_vgm_ibin<j)+1

 ! Memory deallocation
 deallocate(mpair1,mpair2,mpair3,mpair4)

 ! Memory allocation
 allocate(cr(n_out),cv(n_out),tmp(n_out),mpair1_hat(n_out),mpair2_hat(n_out), &
                      mpair3_hat(n_out),mpair4_hat(n_out))

 ! Initialize to zero
 mpair1_hat = 0.0_dp; mpair2_hat = 0.0_dp
 mpair3_hat = 0.0_dp; mpair4_hat = 0.0_dp

 mpair1_hat(vind1) = gamma_hat(ind1)
 mpair2_hat(vind2) = gamma_hat(ind2)
 mpair3_hat(vind3) = gamma_hat(ind3)
 mpair4_hat(vind4) = gamma_hat(ind4)

 ! Memory deallocation
 deallocate(vind1,vind2,vind3,vind4,ind1,ind2,ind3,ind4)

 ! Compute cr
 cr = (-mpair1_hat-mpair2_hat+mpair3_hat+mpair4_hat)/ &
         (2.0_dp*sqrt(abs(gamma_hat(ib)*gamma_hat(jb))))!<--- (2.0_dp*sqrt((gamma_hat(ib)*gamma_hat(jb))))

 tmp=merge(abs(cr),(/(1.0_dp,k=1,n_out)/),abs(cr)<1.0_dp)
 cr=merge(tmp,-tmp,cr>=0.0_dp)
 
 ! Memory deallocation
 deallocate(mpair1_hat,mpair2_hat,mpair3_hat,mpair4_hat)

 ! Compute mult (using the function "approx" to compute a linear interpolation)
 allocate(mask(n_out))

 ! Memory allocation and compute the index array "ind"
 mask=(abs(cr)<=0.96_dp)
 m=count(mask)
 allocate(ind(m),auxcv(m))
 ind=pack(aux,mask)
 cv=1.0_dp
 call approx_linear(rho_grid,f_grid,p+1,abs(cr(ind)),m,f_grid(1),f_grid(p+1),auxcv)
 cv(ind)=auxcv
 !cv=cv*0.172402_dp*sqrt(sqrt((gamma_hat(ib)* gamma_hat(jb)))) !<-------------
 cv=cv*0.172402_dp*sqrt(sqrt(abs(gamma_hat(ib)* gamma_hat(jb))))
 mean_cv=sum(cv)/n_out

 ! Memory deallocation
 deallocate(mask,ind,cr,cv,aux,auxcv)

end subroutine cov_bin_fun

      
subroutine diag_cov_bin_fun(B,n,p,vgm_ibin,vgm_ipair,gamma_hat,mean_cv)

 implicit none

 interface 
   subroutine cov_bin_fun(B,n,p,i,j,vgm_ibin,vgm_ipair,gamma_hat,mean_cv)
     integer,parameter:: dp = KIND(1.0d0)
     integer,intent(in)::B,n,p,i,j
     integer,intent(in):: vgm_ipair(n,2),vgm_ibin(n)
     real(kind=dp),intent(in):: gamma_hat(B)
     real(kind=dp),intent(inout):: mean_cv
   end subroutine cov_bin_fun
 end interface

 integer,parameter:: dp = KIND(1.0d0)
 integer,intent(in)::B,n,p
 integer,intent(in):: vgm_ipair(n,2),vgm_ibin(n)
 real(kind=dp),intent(in):: gamma_hat(B)
 real(kind=dp),intent(inout):: mean_cv(B)
 integer:: i

 do i=1,B
   call cov_bin_fun(B,N,p,i,i,vgm_ibin,vgm_ipair,gamma_hat,mean_cv(i))
 enddo

end subroutine diag_cov_bin_fun


subroutine full_cov_bin_fun(B,n,p,vgm_ibin,vgm_ipair,gamma_hat,mean_cv)

 implicit none

 interface 
   subroutine cov_bin_fun(B,n,p,i,j,vgm_ibin,vgm_ipair,gamma_hat,mean_cv)
     integer,parameter:: dp = KIND(1.0d0)
     integer,intent(in)::B,n,p,i,j
     integer,intent(in):: vgm_ipair(n,2),vgm_ibin(n)
     real(kind=dp),intent(in):: gamma_hat(B)
     real(kind=dp),intent(inout):: mean_cv
   end subroutine cov_bin_fun
 end interface

 integer,parameter:: dp = KIND(1.0d0)
 integer,intent(in)::B,n,p
 integer,intent(in):: vgm_ipair(n,2),vgm_ibin(n)
 real(kind=dp),intent(in):: gamma_hat(B)
 real(kind=dp),intent(inout):: mean_cv(B,B)
 integer:: i,j

 do i=1,B
   do j=i,B
     call cov_bin_fun(B,N,p,i,j,vgm_ibin,vgm_ipair,gamma_hat,mean_cv(i,j))
     mean_cv(j,i)=mean_cv(i,j)
   enddo
 enddo

end subroutine full_cov_bin_fun



