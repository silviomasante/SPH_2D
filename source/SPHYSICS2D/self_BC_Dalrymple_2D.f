c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Dr Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
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

      subroutine self(j1,kind_p1,ini_kind_p2)
c
      include 'common.2D'
      

      do kind_p2=ini_kind_p2,2
       Sr002=0
       Sn00=0
       jj_start = 1
       if(kind_p1.gt.kind_p2)then
         jj_start = nplink_max + 1      
       endif

       do ii=1,nc(j1,kind_p1)
         i = ibox(j1,kind_p1,ii)

         
         if(kind_p1.eq.kind_p2)then
           jj_start = jj_start + 1
         end if
         
         do jj=jj_start,nc(j1,kind_p2)
            j = ibox(j1,kind_p2,jj)

            drx = xp(i) - xp(j)
            drz = zp(i) - zp(j)
            rr2 = drx*drx + drz*drz

	    if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then

            dux = up(i) - up(j)
            duz = wp(i) - wp(j)
             Sr002=Sr002+rr2
             Sn00=Sn00+1

c            Calculating kernel & Normalized Kernel Gradient
             call kernel(drx,drz,i,j,j1,j1,rr2) 
             call kernel_correction(i,j)

c
c  ...  average density
c

            robar  = 0.5*( rhop(i) + rhop(j) )
            one_over_rhobar = 2.0/(rhop(i) + rhop(j))
c

		cbar = 0.5*( cs(i)   + cs(j)   ) 
c
c  ...  inner product rv
c
            dot = drx*dux + drz*duz

c	  Used to calculate the time step due to viscosity


            visc_dt=max(dot/(rr2 + eta2),visc_dt) 


c
c  ...  pressure and viscous force (1992; 3.3)
c
c         pm(j) is mass of particle j
c
            p_v = pr(i) + pr(j)  !+  pi_visc



c  		 for tensile correction (Monaghan , JCP.  159 (2000) 290- 311)
c
c	Only to be activated with cubic spline kernel 
c

		if(index_tensile.eq.1) then   


c   _____ Tensile correction

             fab=Wab*od_Wdeltap              !NOTE: We'll use a non-normalized
             fab=fab*fab                     !kernel to calculate tensile correction
             fab=fab*fab                     !  It's the (Wab/Wdeltap)**4  of Monaghan's paper

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

             p_v = p_v+ (Ra+Rb)*fab

		endif

          ax(i) = ax(i) - pm(j) * p_v * frxi
          az(i) = az(i) - pm(j) * p_v * frzi

          ax(j) = ax(j) + pm(i) * p_v * frxj
          az(j) = az(j) + pm(i) * p_v * frzj

             call gradients_calc(i,j,dux,duz)               

             call viscosity(dot,drx,drz,dux,duz,rr2,
     +             cbar,robar,one_over_rhobar,i,j,j1,j1,term2i,term2j)



c   ...         Thermal Energy
c             (Monaghan, JCP 110 (1994) 399- 406)


              term1i=0.5 * p_v *( frxi*dux+frzi*duz)
              term1j=0.5 * p_v *( frxj*dux+frzj*duz)

              aTE(i)=aTE(i)+pm(j) * (term1i+term2i)
              aTE(j)=aTE(j)+pm(i) * (term1j+term2j)


c  ...  density acceleration (1992; 3.9)
c

          dot2i = dux*frxi + duz*frzi
          dot2j = dux*frxj + duz*frzj
          ar(i) = ar(i) + pm(j)*dot2i
          ar(j) = ar(j) + pm(i)*dot2j

c
c  ...  XSPH correction
c
          pmj_Wab_over_rhobar = pm(j)*Wab*one_over_rhobar
          ux(i) = ux(i) - dux*pmj_Wab_over_rhobar  !pm(j) * dux * Wab / robar
          wx(i) = wx(i) - duz*pmj_Wab_over_rhobar  !pm(j) * duz * Wab / robar


          pmi_Wab_over_rhobar = pm(i)*Wab*one_over_rhobar
          ux(j) = ux(j) + dux*pmi_Wab_over_rhobar   !pm(i) * dux * Wab / robar
          wx(j) = wx(j) + duz*pmi_Wab_over_rhobar   !pm(i) * duz * Wab / robar

c         SHIFT CORRECTION        

c          denom=rr2*sqrt(rr2)*Sn00*Sn00
c          num=BETA_sh_corr*Svmax*Sr002
c          coefficiente=num/denom
          
c          ux(i) = ux(i) - coefficiente*drx
          
c          wx(i) = wx(i) - coefficiente*drz
c          ux(j) = ux(j) + coefficiente*drx
          
c          wx(j) = wx(j) + coefficiente*drz        
 
c       absorbed fluid velocity
 
       
 
        VfX(i)=-coeffvel(i)*(-dPcdx(i))
        VfZ(i)=-coeffvel(i)*(-dPcdz(i)-rho0*grz)         
        VfX(j)=-coeffvel(j)*(-dPcdx(j))
        VfZ(j)=-coeffvel(j)*(-dPcdz(j)-rho0*grz) 
         
         versx=drx/sqrt(rr2)
         versz=drz/sqrt(rr2)
         
         
       diff(i)=((-VfX(j)*versx-VfZ(j)*versz)*(saturazione(j))**7)/rr2
c     + sqrt(rr2)

       diff(j)=((-VfX(i)*versx-VfZ(i)*versz)*(saturazione(i))**7)/rr2
c     + sqrt(rr2)  
         
         dmpItest=pm(i)+(dmpdt(i)+diff(i)*pVol(j)*(fmp(i)-fmp(j)))*dt2
         dmpJtest=pm(j)+(dmpdt(j)+diff(i)*pVol(i)*(fmp(j)-fmp(i)))*dt2
         

        if (((dmpItest.ge.0.d0).or.(dmpItest.le.125.d0))
     + .and.((dmpJtest.ge.0.d0).or.(dmpJtest.le.125.d0))) then
         dmpdt(i)= dmpdt(i)+diff(i)*pVol(j)*(fmp(i)-fmp(j))
         dmpdt(j)= dmpdt(j)+diff(i)*pVol(i)*(fmp(j)-fmp(i))
        endif
        
        dmpItest=pm(i)+(dmpdt(i)+diff(j)*pVol(j)*(fmp(i)-fmp(j)))*dt2
        dmpJtest=pm(j)+(dmpdt(j)+diff(j)*pVol(i)*(fmp(j)-fmp(i)))*dt2
         

        if (((dmpItest.ge.0.d0).or.(dmpItest.le.125.d0))
     + .and.((dmpJtest.ge.0.d0).or.(dmpJtest.le.125.d0))) then
         dmpdt(i)= dmpdt(i)+diff(j)*pVol(j)*(fmp(i)-fmp(j))
         dmpdt(j)= dmpdt(j)+diff(j)*pVol(i)*(fmp(j)-fmp(i))
        endif
        
c             ...  Vorticity calculation

        if(ipoute.eq.1.and.i_vort.eq.1.and.i.gt.nb.and.j.gt.nb)then
                  call vorticity_calc(i,j,dux,duz)
        endif 

	  endif
        enddo
      enddo
c      write(*,*) 'AHHHHHHHH  ',coefficiente*drx/ux(i)
	enddo

      end
c

