!-------------------------------------------------------------------------------------------------
!		Subroutine for calculating the positive cell-interface flux using DPENO-AO9 scheme
!-------------------------------------------------------------------------------------------------
module DPENOAO9
	implicit none
	integer:: k, r
	integer:: klist(2:8)
	real(Kind=8):: delta(-3:32)
	real(Kind=8):: alpha(2,2:9)
	real(Kind=8):: dp(25)
	real(Kind=8):: Max_alpha
	integer, parameter :: p=9
	integer, parameter :: kp=(p-1)/2
	integer, parameter :: p_delta(8) = [0, 7, 13, 18, 22, 25, 27, 28]
	integer, parameter :: p_c(9) = [1, 2, 4, 7, 11, 15, 18, 20, 21]
	real(Kind=8), parameter:: kc=3.1415926535897932385/2.
	real(Kind=8), parameter:: sin_r_kc_half(1:p)=(/(sin(kc/2.)**(1-r),r=1,p,1)/)
	real(kind=8), parameter :: c(25) = sqrt([ &
	1./12.,                           &
	1./45.,                           1./45.,                   &
	59./6048.,                        17./6048.,                59./6048.,              &
	1213./226800.,                    163./226800.,             163./226800.,           1213./226800.,     &
	22789./6842880.,                  8807./34214400.,          4121./34214400.,        8807./34214400.,   22789./6842880., &
	71998./638512875.,                38891./1277025750.,       38891./1277025750.,     71998./638512875., &
	5187401./523069747200.,           14981203./2615348736000., 5187401./523069747200., &
	240391./166728481920.,            240391./166728481920., &
	44734915633./154821036883968000.])
	real(Kind=8), parameter:: a(-4:4, -4:0, 0:4)=reshape(&
	[12./60.,     -63./60.,     137./60.,     -163./60.,    137./60.,      0.,            0.,           0.,          0.,        &
	0.,          -3./12.,      13./12.,      -23./12.,     25./12.,       0.,            0.,           0.,          0.,        &
	0.,          0.,           2./6.,        -7./6.,       11./6.,        0.,            0.,           0.,          0.,        &
	0.,          0.,           0.,           -1./2.,       3./2.,         0.,            0.,           0.,          0.,        &
	0.,          0.,           0.,           0.,           1./1.,         0.,            0.,           0.,          0.,        &
	2./60.,      -13./60.,     37./60.,      -63./60.,     87./60.,       10./60.,       0.,           0.,          0.,        &
	0.,          -3./60.,      17./60.,      -43./60.,     77./60.,       12./60.,       0.,           0.,          0.,        &
	0.,          0.,           1./12.,       -5./12.,      13./12.,       3./12.,        0.,           0.,          0.,        &
	0.,          0.,           0.,           -1./6.,       5./6.,         2./6.,         0.,           0.,          0.,        &
	0.,          0.,           0.,           0.,           1./2.,         1./2.,         0.,           0.,          0.,        &
	4./420.,     -31./420.,    109./420.,    -241./420.,   459./420.,     130./420.,     -10./420.,    0.,          0.,        &
	0.,          -1./60.,      7./60.,       -23./60.,     57./60.,       22./60.,       -2./60.,      0.,          0.,        &
	0.,          0.,           2./60.,       -13./60.,     47./60.,       27./60.,       -3./60.,      0.,          0.,        &
	0.,          0.,           0.,           -1./12.,      7./12.,        7./12.,        -1./12.,      0.,          0.,        &
	0.,          0.,           0.,           0.,           2./6.,         5./6.,         -1./6.,       0.,          0.,        &
	3./840.,     -27./840.,    113./840.,    -307./840.,   743./840.,     365./840.,     -55./840.,    5./840.,     0.,        &
	0.,          -3./420.,     25./420.,     -101./420.,   319./420.,     214./420.,     -38./420.,    4./420.,     0.,        &
	0.,          0.,           1./60.,       -8./60.,      37./60.,       37./60.,       -8./60.,      1./60.,      0.,        &
	0.,          0.,           0.,           -3./60.,      27./60.,       47./60.,       -13./60.,     2./60.,      0.,        &
	0.,          0.,           0.,           0.,           3./12.,       13./12.,       -5./12.,       1./12.,      0.,        &
	4./2520.,    -41./2520.,   199./2520.,   -641./2520.,  1879./2520.,   1375./2520.,   -305./2520.,  55./2520.,   -5./2520., &
	0.,          -3./840.,    29./840.,      -139./840.,   533./840.,     533./840.,     -139./840.,   29./840.,    -3./840.,  &
	0.,          0.,          4./420.,       -38./420.,    214./420.,     319./420.,     -101./420.,   25./420.,    -3./420.,  &
	0.,          0.,          0.,            -2./60.,      22./60.,       57./60.,       -23./60.,     7./60.,      -1./60.,   &
	0.,          0.,          0.,            0.,           12./60.,       77./60.,       -43./60.,     17./60.,     -3./60.    &
	],[9,5,5])
end module DPENOAO9

subroutine flux_DPENOAO_9_P(f_j, numericalFlux)
	use DPENOAO9
	implicit none
	real(Kind=8), intent(in) :: f_j(-kp:kp)
	real(Kind=8), intent(out) :: numericalFlux
	integer :: Lp,Rp

!-------------------------------------------------------
	! Init delta
	do k = 1-kp, kp
		delta(k) = f_j(k) - f_j(k-1)
	end do
	
!-------------------------------------------------------
	! Forward pass procedure
	do r = 1, kp
		do k = r-kp, 0
			delta(p_delta(r+1)+k+1) = delta(p_delta(r)+k+1) - delta(p_delta(r)+k)
		end do
		dp(p_c(r+1)) = dp(p_c(r)) + c(p_c(r))*abs(delta(p_delta(r)))
		do k= 1, r-1
			delta(p_delta(r+1)+k+1) = delta(p_delta(r)+k+1) - delta(p_delta(r)+k)
			dp(p_c(r+1)+k) = min(dp(p_c(r)+k)   + c(p_c(r)+k)  *abs(delta(p_delta(r)+k)), &
								 dp(p_c(r)+k-1) + c(p_c(r)+k-1)*abs(delta(p_delta(r)+k)))
		enddo
		delta(p_delta(r+1)+r+1) = delta(p_delta(r)+r+1) - delta(p_delta(r)+r)
		dp(p_c(r+1)+r) = dp(p_c(r)+r-1) + c(p_c(r)+r-1)*abs(delta(p_delta(r)+r))
		do k = r+1, kp-1
			delta(p_delta(r+1)+k+1) = delta(p_delta(r)+k+1) - delta(p_delta(r)+k)
		end do
	end do
	
	do r = kp+1, 2*kp-1
		do k = r-kp, kp-1
			delta(p_delta(r+1)+k+1) = delta(p_delta(r)+k+1) - delta(p_delta(r)+k)
			dp(p_c(r+1)+k) = min(dp(p_c(r)+k)   + c(p_c(r)+k)  *abs(delta(p_delta(r)+k)), &
								 dp(p_c(r)+k-1) + c(p_c(r)+k-1)*abs(delta(p_delta(r)+k)))
		end do
		dp(p_c(r+1)+kp) = min(dp(p_c(r)+kp)   + c(p_c(r)+kp)  *abs(delta(p_delta(r)+kp)), &
							  dp(p_c(r)+kp-1) + c(p_c(r)+kp-1)*abs(delta(p_delta(r)+kp)))
	end do

!-------------------------------------------------------
	! Backward pass procedure
	if(dp(p_c(p-1)+kp) .lt. dp(p_c(p-1)+kp-1))then
		klist(p-1)=kp
		alpha(1,p-1)=c(p_c(p-1)+kp)*abs(delta(p_delta(p-1)+kp))*sin_r_kc_half(p-1)
	else
		klist(p-1)=kp-1
		alpha(1,p-1)=c(p_c(p-1)+kp)*abs(delta(p_delta(p-1)+kp))*sin_r_kc_half(p-1)
	endif
	
	if (dp(p_c(p-2)+klist(p-1)-1)+c(p_c(p-2)+klist(p-1)-1)*abs(delta(p_delta(p-2)+klist(p-1))) .gt.  &
		dp(p_c(p-1)+klist(p-1)))then
		klist(p-2)=klist(p-1)
		alpha(1,p-2)=(dp(p_c(p-1)+klist(p-1))-dp(p_c(p-2)+klist(p-2)))*sin_r_kc_half(p-2)
		alpha(2,p-2)=maxval(alpha(1,p-2:p-1))
	else
		klist(p-2)=klist(p-1)-1
		alpha(1,p-2)=(dp(p_c(p-1)+klist(p-1))-dp(p_c(p-2)+klist(p-2)))*sin_r_kc_half(p-2)
		alpha(2,p-2)=maxval(alpha(1,p-2:p-1))
	endif

	do r=p-2,kp+2,-1
		if (dp(p_c(r-1)+klist(r)-1)+c(p_c(r-1)+klist(r)-1)*abs(delta(p_delta(r-1)+klist(r))) .gt.  &
			dp(p_c(r)+klist(r)))then
            klist(r-1)=klist(r)
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		else
            klist(r-1)=klist(r)-1
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		endif
	enddo

	do r=kp+1,3,-1
		if (klist(r) == 0)then
            klist(r-1)=klist(r)
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		elseif (klist(r) == r-1)then
            klist(r-1)=klist(r)-1
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		elseif (dp(p_c(r-1)+klist(r)-1)+c(p_c(r-1)+klist(r)-1)*abs(delta(p_delta(r-1)+klist(r))) .gt.  &
				dp(p_c(r)+klist(r)))then
            klist(r-1)=klist(r)
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		else
            klist(r-1)=klist(r)-1
			alpha(1,r-1)=(dp(p_c(r)+klist(r))-dp(p_c(r-1)+klist(r-1)))*sin_r_kc_half(r-1)
			alpha(2,r-1)=min(maxval(alpha(1,r-1:r)),alpha(2,r))
		endif
	enddo
	
! -------------------------------------------------------
	! Compute the numerical flux with adaptive order
	Max_alpha=dp(p_c(2)+klist(2))
	do r=3,p-1
		Max_alpha=max(Max_alpha,alpha(1,r-1))
		if(alpha(2,r-1) .gt. Max_alpha) then
			Rp=klist(r); Lp=Rp+1-r
			numericalFlux=sum(a(Lp:Rp,Lp,Rp)*f_j(Lp:Rp))
			return
		endif
	enddo
	numericalFlux=sum(a(-kp:kp,-kp,kp)*f_j(-kp:kp))
end
