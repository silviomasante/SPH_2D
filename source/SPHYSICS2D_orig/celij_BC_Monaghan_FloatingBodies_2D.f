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

      subroutine celij(j1,j2,kind_p1,ini_kind_p2,lx2)
c
      include 'common.2D'  
      
      
      !do kind_p2=ini_kind_p2,2
      !- Changes looping to include BP-BP interactions
      if(kind_p1.eq.1)then
        kind_p2_start = 1
      else
        kind_p2_start = ini_kind_p2
      endif
      do kind_p2=kind_p2_start,2
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
               
c            ...  average density

             robar  = 0.5*( rhop(i) + rhop(j) )
             one_over_rhobar = 2.0/( rhop(i) + rhop(j) )
c
             cbar = 0.5*( cs(i) + cs(j) ) 

c            ...  inner product rv
c
             dot = drx*dux + drz*duz

c            Used to calculate the time step duu to viscosity

             visc_dt=max(dot/(rr2 + eta2),visc_dt) 
c
c            ...  pressure and viscous force (Monaghan 1992; Ann. Rev.
c                 Astron. Astrop. 30. Formula 3.3)
c                 pm(j) is mass of particle j
c
             p_v = pr(i) + pr(j)


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


             if(i.gt.nb.and.j.gt.nb) then   !both fluid particles

               ax(i) = ax(i) - pm(j) * p_v * frxi
               az(i) = az(i) - pm(j) * p_v * frzi


               ax(j) = ax(j) + pm(i) * p_v * frxj
               az(j) = az(j) + pm(i) * p_v * frzj


               call gradients_calc(i,j,dux,duz)                                       

               call viscosity(dot,drx,drz,dux,duz,rr2,
     +               cbar,robar,one_over_rhobar,i,j,j1,j2,term2i,term2j)


c              ...  density acceleration
c              using the derivative of the kernel, not the kernel itself

               dot2i = dux*frxi + duz*frzi
               dot2j = dux*frxj + duz*frzj
               ar(i) = ar(i) + pm(j)*dot2i
               ar(j) = ar(j) + pm(i)*dot2j

c   ...        Thermal Energy

               term1i=0.5 * p_v *( frxi*dux+frzi*duz)
               term1j=0.5 * p_v *( frxj*dux+frzj*duz)

               aTE(i)=aTE(i)+pm(j) * (term1i+term2i)
               aTE(j)=aTE(j)+pm(i) * (term1j+term2j)

c
c              ...  XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
c

               pmj_Wab_over_rhobar = pm(j)*Wab*one_over_rhobar
               ux(i) = ux(i) - dux*pmj_Wab_over_rhobar  !pm(j) * dux * Wab / robar
               wx(i) = wx(i) - duz*pmj_Wab_over_rhobar  !pm(j) * duz * Wab / robar


               pmi_Wab_over_rhobar = pm(j)*Wab*one_over_rhobar
               ux(j) = ux(j) + dux*pmi_Wab_over_rhobar   !pm(i) * dux * Wab / robar
               wx(j) = wx(j) + duz*pmi_Wab_over_rhobar   !pm(i) * duz * Wab / robar


c              ...  Vorticity calculation


             if(ipoute.eq.1.and.i_vort.eq.1.and.i.gt.nb.and.j.gt.nb)then
                  call vorticity_calc(i,j,dux,duz)
             endif  

             else if(i.gt.nb.and.j.le.nb)then			!j is boundary particle
               call MonaghanBC(i,j,drx,drz,dot,dux,duz,
     &                             fxbMON,fzbMON)   
               ax(i) = ax(i) + fxbMON
               az(i) = az(i) + fzbMON

               !- Viscous wall force -
               temp = pVol(i)*visc_wall*
     &           ((drx*frxi+drz*frzi)/(rr2 + eta2))
               ax(i) = ax(i) + temp*dux
               az(i) = az(i) + temp*duz

                          
               if(j.gt.nbfm)then  !Opposite Forces on free-moving object (Not per unit mass, see subroutine whole_body_motion)
                 ax(j) = ax(j) - fxbMON*pm(i)
                 az(j) = az(j) - fzbMON*pm(i)

                 !- Viscous wall force -
                 temp = pVol(i)*visc_wall*
     &             ((drx*frxj+drz*frzj)/(rr2 + eta2))
                 ax(j) = ax(j) - temp*dux*pm(i)
                 az(j) = az(j) - temp*duz*pm(i)

               end if

             else if(j.gt.nb.and.i.le.nb)then			!i is boundary particle
               call MonaghanBC(j,i,-drx,-drz,dot,-dux,-duz,
     &                             fxbMON,fzbMON)   
               ax(j) = ax(j) + fxbMON
               az(j) = az(j) + fzbMON

               !- Viscous wall force -
               temp = pVol(j)*visc_wall*
     &           ((drx*frxj+drz*frzj)/(rr2 + eta2))
               ax(j) = ax(j) - temp*dux
               az(j) = az(j) - temp*duz

               
               if(i.gt.nbfm)then  !Opposite Forces on free-moving object (Not per unit mass, see subroutine whole_body_motion)
                 ax(i) = ax(i) - fxbMON*pm(j)
                 az(i) = az(i) - fzbMON*pm(j)

                 !- Viscous wall force -
                 temp = pVol(i)*visc_wall*
     &             ((drx*frxi+drz*frzi)/(rr2 + eta2))
                 ax(i) = ax(i) + temp*dux*pm(j)
                 az(i) = az(i) + temp*duz*pm(j)

               end if

            else if(i.lt.nbp1.and.j.lt.nbp1)then !BP-BP
              if(i.gt.nbfm)then !i is free-moving BP
                if(j.gt.nbfm)then  !j is free-moving BP
                  if(num_FB.gt.1)then
                    !Loop to prevent BP within a Floating Body interacting with a BP within the same Floating Body
                    i_num_FB = 1
                    i_BP_BP_interaction = 0
                    id_num_FB_i = i_FB_Pointer_Info(i)
                    id_num_FB_j = i_FB_Pointer_Info(j)
                    do while (i_BP_BP_interaction.lt.1)
                      if(i_num_FB.lt.num_FB)then
                        i_num_FB = i_num_FB + 1
c                        if(i.gt.nb_FB(i_num_FB-1).and.
c     &                     i.le.nb_FB(i_num_FB)  .and.
c     &                    (j.le.nb_FB(i_num_FB-1).or.
c     &                     j.gt.nb_FB(i_num_FB)      ) )then
                         if(id_num_FB_i.ne.id_num_FB_j)then
                          call MonaghanBC
     &                         (i,j,drx,drz,dot,dux,duz,
     &                                     fxbMON,fzbMON)
                          ax(i) = ax(i) + fxbMON*pm(j)
                          az(i) = az(i) + fzbMON*pm(j)

                          ax(j) = ax(j) - fxbMON*pm(j)
                          az(j) = az(j) - fzbMON*pm(j)

                          i_BP_BP_interaction = 1
                        endif
                      else
                        i_BP_BP_interaction = 1
                      endif
                    enddo
                  endif
                else  !j is prescribed BP
                  call MonaghanBC(i,j,drx,drz,dot,dux,duz,
     &                             fxbMON,fzbMON)
                  ax(i) = ax(i) + fxbMON*pm(j)
                  az(i) = az(i) + fzbMON*pm(j)

                end if    !End of:  if(j.gt.nbfm)then
              else if(j.gt.nbfm)then   !j is free-moving BP, i is prescibed BP
                call MonaghanBC(j,i,-drx,-drz,dot,-dux,-duz,
     &                          fxbMON,fzbMON)
                ax(j) = ax(j) + fxbMON*pm(i)
                az(j) = az(j) + fzbMON*pm(i)

              endif   !End of:  if(i.gt.nbfm)then
            endif  ! Interaction fluid-fluid or fluid-boundary
           
           endif ! if q<2
         enddo ! Box jj
        enddo ! Box ii
       endif  ! Box jj is not empty
      enddo   ! Kind of particle

      end