c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.

      subroutine celij(j1,j2,kind_p1,ini_kind_p2,lx2)
c
      include 'common.2D'  
      

      !- Changes looping to include BP-BP interactions
      if(kind_p1.eq.1)then
        kind_p2_start = 1
      else
        kind_p2_start = ini_kind_p2
      endif
      
      do kind_p2=kind_p2_start,2
c      do kind_p2=ini_kind_p2,2
        if(nc(j2,kind_p2).ne.0) then

        do ii=1,nc(j1,kind_p1)
          i = ibox(j1,kind_p1,ii)
         
          do jj=1,nc(j2,kind_p2)
           j = ibox(j2,kind_p2,jj)
            
            drx = xp(i) - xp(j)
            drz = zp(i) - zp(j)

            call periodicityCorrection(i,j,drx,drz,lx2)

            rr2 = drx*drx + drz*drz

            if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then
             dux = up(i) - up(j)
             duz = wp(i) - wp(j)

c            Calculating kernel & Normalized Kernel Gradient
             call kernel(drx,drz,i,j,j1,j2,rr2) 
             call kernel_correction(i,j)
   
             
c  ...  average density

            robar  = 0.5*( rhop(i) + rhop(j) )
            one_over_rhobar = 2.0/(rhop(i) + rhop(j))
	      cbar = 0.5*( cs(i)   + cs(j)   ) 

c  ...  inner product rv
c
            dot = drx*dux + drz*duz

c	  Used to calculate the time step due to viscosity

            visc_dt=max(dot/(rr2 + eta2),visc_dt) 


c  ...  pressure and viscous force (Monaghan 1992; Ann. Rev.
c			        Astron. Astrop. 30. Formula 3.3)
c         pm(j) is mass of particle j
c
            p_v = pr(i) + pr(j) !+ pi_visc


c  	Tensile correction (Monaghan , JCP.  159 (2000) 290- 311)
c		
c	Only to be activated with cubic spline kernel (Mon (8/21/02)
	
	     if (index_tensile.eq.1) then

c              ____ Tensile correction 
 
               fab=Wab*od_Wdeltap     !NOTE: We'll use a non-normalized
               fab=fab*fab            !kernel to calculate tensile correction
               fab=fab*fab            !It's the (Wab/Wdeltap)**4  of Monaghan's paper

             if (p(i).gt.0) then
                  Ra= 0.01 *pr(i)
             else
                  Ra= 0.2 *abs(pr(i))
             endif

             if (p(j).gt.0) then
                  Rb= 0.01 *pr(j)
             else
                  Rb= 0.2 *abs(pr(j))
             endif

             R=Ra+Rb
             p_v = p_v+ R*fab

	    endif

c            ax(i) = ax(i) - pm(j) * p_v * frxi
c            az(i) = az(i) - pm(j) * p_v * frzi
c
c            ax(j) = ax(j) + pm(i) * p_v * frxj
c            az(j) = az(j) + pm(i) * p_v * frzj
                 !- Debug code -
                i_particleCheck = -23197
                i_rankCheck = 0
                dpN = 0.01*vnorm_mass
599             format(a8,2(1x,i7),i3,2(1x,i2),2(1x,f12.6))            
                   !Place in separate ifs            
                  !Place in separate ifs            

          if(i.gt.nb)then
            ax(i) = ax(i) - pm(j) * p_v * frxi
            az(i) = az(i) - pm(j) * p_v * frzi
                 if(0.eq.i_rankCheck.and.i.eq.i_particleCheck)then
                    idrx = nint(drx/dpN)
                    idrz = nint(drz/dpN)
                    write(*,599)
     &              'pC:I,j ',i,j,id_call,idrx,idrz,
     &                       -pm(j)*p_v*frxi,-pm(j)*p_v*frzi
                  endif 
          elseif(i.gt.nbfm)then
            ax(i) = ax(i) - pm(j) * p_v * frxi * pm(i)
            az(i) = az(i) - pm(j) * p_v * frzi * pm(i)
          endif

          if(j.gt.nb)then
            ax(j) = ax(j) + pm(i) * p_v * frxj
            az(j) = az(j) + pm(i) * p_v * frzj
                 if(0.eq.i_rankCheck.and.j.eq.i_particleCheck)then
                    idrx = -nint((drx)/dpN)
                    idrz = -nint((drz)/dpN)
                    write(*,599)
     &              'pC:i,J ',i,j,id_call,idrx,idrz,
     &                        pm(i)*p_v*frxj, pm(i)*p_v*frzj
                  endif 
          elseif(j.gt.nbfm)then
            ax(j) = ax(j) + pm(i) * p_v * frxj * pm(j)
            az(j) = az(j) + pm(i) * p_v * frzj * pm(j)
          endif

            call gradients_calc(i,j,dux,duz)               

            call viscosity(dot,drx,drz,dux,duz,rr2,
     +             cbar,robar,one_over_rhobar,i,j,j1,j2,term2i,term2j)

c   ...         Thermal Energy
c		(Monaghan, JCP 110 (1994) 399- 406)
c		MON(7/25/01)

           term1i=0.5 * p_v *( frxi*dux+frzi*duz)
           term1j=0.5 * p_v *( frxj*dux+frzj*duz)

           aTE(i)=aTE(i)+pm(j) * (term1i+term2i)
           aTE(j)=aTE(j)+pm(i) * (term1j+term2j)


c  ...  density acceleration (Monaghan 1992; Ann. Rev. Astron. Astrop. 30. Formula 3.9)
c       using the derivative of the kernel, not the kernel itself

           dot2i = dux*frxi + duz*frzi
           dot2j = dux*frxj + duz*frzj
           ar(i) = ar(i) + pm(j)*dot2i
           ar(j) = ar(j) + pm(i)*dot2j

c
c  ...  XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
c

             pmj_Wab_over_rhobar = pm(j)*Wab*one_over_rhobar
             ux(i) = ux(i) - dux*pmj_Wab_over_rhobar  !pm(j) * dux * Wab / robar 
             wx(i) = wx(i) - duz*pmj_Wab_over_rhobar  !pm(j) * duz * Wab / robar


             pmi_Wab_over_rhobar = pm(i)*Wab*one_over_rhobar
             ux(j) = ux(j) + dux*pmi_Wab_over_rhobar   !pm(i) * dux * Wab / robar
             wx(j) = wx(j) + duz*pmi_Wab_over_rhobar   !pm(i) * duz * Wab / robar

c             ...  Vorticity calculation

             if(ipoute.eq.1.and.i_vort.eq.1.and.i.gt.nb.and.j.gt.nb)then
                  call vorticity_calc(i,j,dux,duz)
             endif 

c            if(i.gt.nbfm.and.i.lt.nbp1.and.j.gt.nb)then  !opposite forces on free-moving object
c              ax(i) = ax(i) + pm(j) * p_v * frxi *pm(i)
c              az(i) = az(i) + pm(j) * p_v * frzi *pm(i)
c            elseif(j.gt.nbfm.and.j.lt.nbp1.and.i.gt.nb)then
c              ax(j) = ax(j) + pm(i) * p_v * frxj *pm(j)
c              az(j) = az(j) + pm(i) * p_v * frzj *pm(j)
c            endif

	      endif ! if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then
         enddo  ! do jj=1,nc(j2,kind_p2)
        enddo  ! do ii=1,nc(j1,kind_p1)
	 endif  ! if(nc(j2,kind_p2).ne.0) then
	enddo  ! do kind_p2=kind_p2_start,2

      end

