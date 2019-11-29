!id ising model without seeds
        
    program isingmodel
	implicit none
	
	double precision:: T,r1,ptrial, Estart,Eflip,Efinal,E_T,et,dT,dE,m,Et2,Cv,avg_E,avg2_E,mm,c2,ef,e1,e2,avg_mag
	integer ::i, row, right ,left,nmc,iT,c,c1,nn,nnn,k
	integer, parameter :: J = -1   
	integer, dimension(:),allocatable :: spin(:) 

	print*, 'give lattice sites'
	read*,nn
	print*,'give mc step'
	read*,nnn
	allocate(spin(nn))
	spin=0
	!av_energy=0
	
	
	
	do i=1,nn
		call random_number(r1)
	!	print*, random1
	  if (r1 .lt. 0.5) then
		spin(i) = +1
		else 
		spin(i) = -1
	  end if
	end do
	!print*,spin
!	print*,spin
	 T=2.010d0
	
	do iT=1,30
	t=t+dt
	
	if(iT.le.10) T=T-0.1
	if ((iT.gt.10).and.(iT.le.20)) T=T-0.09
	if ((iT.gt.20).and.(iT.le.30)) T=T-0.009
	
	
	 e1=0.0d0
	 e2=0.0d0
     m=0.0d0
	 c1=0.0d0
	 c2=0.0d0
		do nmc=1,nnn!no. of monte carlo steps
				
			E_T=0.0d0 
			
			c=0.0d0

			do i=1,nn
			
				right=i+1
				left=i-1     
			    if (i.eq.1) left=nn
				if  (i.eq.nn) right=1        
				Estart=J*(spin(i)*spin(right)+spin(i)*spin(left))                         
				Eflip =J*((-spin(i))*spin(right)+(-spin(i))*spin(left))
				dE=Eflip-Estart
               ! write(811,*)  Eflip,Estart        
		        if (dE.le.0) then 
				 spin(i)=-spin(i) 
				 Efinal=Eflip 
				 c1=c1+1
				
				 elseif(de.gt.0)then  
				 call random_number(r1)
				 ptrial = EXP(-dE/T)
				 
              !   write(500,*) de,ptrial
		         if(r1 .lt. ptrial)then
		  		 spin(i) = -spin(i)
		         Efinal=Eflip
		         c1=c1+1
		 !      print*, random1,ptrial
                else
                Efinal=Estart
                 c2=c2+1
                endif 
              E_T=E_T+Efinal 
              Et2=Et2+ Efinal**2
               endif 
          
             enddo
         !    write(101,*) e_t
             if(nmc.gt.(nnn/2))then
                
               ! m=0.0d0
                do i=1,nn
               
                m=m+spin(i) 
                enddo
              do i=1,nn
          	
                right=i+1
				left=i-1     
			    if (i.eq.1) left=nn
				if  (i.eq.nn) right=1  
                ef=J*(spin(i)*spin(right)+spin(i)*spin(left))  
	            e1=e1+(ef/2.0)
	            e2=e2+((ef/2.0)**2)
	           
	         enddo
           endif 
                
         enddo
		 write(400,*) c1,c2
            avg_mag=abs(m)/dble(nn*(nnn/2))
            avg_E=e1/dble(nn*(nnn/2)) 
            avg2_E=e2/dble(nn*(nnn/2)) 
            cv=(avg2_e-(avg_e)**2)/((t**2)*nn)
            write(111,*) t,avg_mag,avg_e ,cv	      
	   enddo
   end program isingmodel

