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

      subroutine pre_self(j1,kind_p1,ini_kind_p2)
c      implicit double precision(a-h,o-z)
      include 'common.2D'

	do kind_p2=ini_kind_p2,2

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

c              if(i.ge.nstart.and.j.ge.nstart)then

               drx = xp(i) - xp(j)
               drz = zp(i) - zp(j)
               rr2 = drx*drx + drz*drz
            
               if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then
               
                  call kernel(drx,drz,i,j,j1,j1,rr2)             
                  
				V_j = pVol(j)  				 		                          			 
                  frxh = V_j*frx
                  frzh = V_j*frz

                  aM_a11(i) = aM_a11(i) - frxh*drx
                  aM_a12(i) = aM_a12(i) - frxh*drz
                  aM_a21(i) = aM_a21(i) - frzh*drx
                  aM_a22(i) = aM_a22(i) - frzh*drz
	
                  V_i = pVol(i)
                  frxh = V_i*frx
                  frzh = V_i*frz

                  aM_a11(j) = aM_a11(j) - frxh*drx
                  aM_a12(j) = aM_a12(j) - frxh*drz
                  aM_a21(j) = aM_a21(j) - frzh*drx
                  aM_a22(j) = aM_a22(j) - frzh*drz 
				   
c                  neigh(i)=neigh(i)+1
c                  neigh(j)=neigh(j)+1				 				     
                endif ! if q<2
c              endif ! if(i.ge.nstart.and.j.ge.nstart)then  
         enddo  ! Box jj
        enddo ! Box ii
      enddo ! Kind of particle
      return
      end
c

