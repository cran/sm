ccccccFORTRAN subroutine for 3d density estimation cccccc

c     This code was written by Stuart Young in connection with his Ph.D.
c     thesis, "Graphics and Inference in Nonparametric Modelling",
c     University of Glasgow, 1996.
c     For wider distribution, some compiler-specific aspects of the code 
c     were modified, initially by B.D. Ripley and subsequently by
c     A.W. Bowman, August 1997.


      subroutine npcont(x,y,z,xx,yy,zz,hx,hy,hz,lng,size,
     &     fmax,q,limit,dh,npoly,c,xcoord,ycoord,zcoord)

c     Takes a 3-dimensional data set and produces sets of points which
c     comprise 3-sided polygons, on the surface of a contour of a
c     3-dimensional nonparametric density estimate. 
      

c     Version 7 : Splus
c     (assumes density height <> contour)
c     (can handle >4 cut-points per cube)


c     Parameters to be supplied:
c     x,y,z 		- data for the three variables
c     xx,yy,zz	- 3 equally-spaced vectors spanning the range
c     of the 3 variables
c     hx,hy,hz	- smoothing parameters for each variable
c     lng		- the number of observations on each variable
c     size		- the length of xx,yy and zz (controls the
c     'roughness' of the plotted contour)
c     fmax		- the maximum value of the density, estimated
c     over all the observed points
c     q		- the value of the density at the contour level
c     chosen
c     limit		- the maximum number of polygons that can be
c     found (for setting array sizes)

c     The following parameters are calculated and returned:
c     dh		- a 3-dimensional matrix of density heights,
c     calculated on the (xx,yy,zz) grid.
c     npoly		- the number of polygons found
c     c		- a 3-dimensional matrix of polygon coordinates
c     xcoord		- vector of x-coordinates for the polygons
c     ycoord		- vector of y-coordinates for the polygons
c     zcoord		- vector of z-coordinates for the polygons


c     ***************************************************************** 
c     Initialization
c     ***************************************************************** 
      integer lng, size, limit, npoly, a, b, d, i, j, k, l
      integer count, numcub, tri(3)
      double precision sgn, fli

      double precision x(lng), y(lng), z(lng)
      double precision xx(size), yy(size), zz(size), hx, hy, hz
      double precision fmax, q, dh(size,size,size), c(limit,3,3)
      double precision xcoord(limit*3), ycoord(limit*3)
      double precision zcoord(limit*3) 
      double precision tmp, sum, dh1, dh2, temp(12,3)


c     ***************************************************************** 
c     Calculate grid of density heights, dh
c     ***************************************************************** 

      do 1 i=1,size
         do 2 j=1,size
            do 3 k=1,size
               tmp=0
               do 4 l=1,lng
                  sum = ((xx(i)-x(l))/hx)**2 + ((yy(j)-y(l))/hy)**2 + 
     &                 ((zz(k)-z(l))/hz)**2
                  tmp=tmp+dexp(-0.5*sum)
   4           continue
               dh(i,j,k)=tmp/(lng*hx*hy*hz*(6.283185**1.5))
c               dh(i,j,k)=tmp/(lng*hx*hy*hz*15.74961)
               dh(i,j,k)=int(dh(i,j,k)*90/fmax)
   3        continue
   2     continue
   1  continue


c     ***************************************************************** 
c     Main loop to check all cubes
c     i,j,k loop ensures all 12 edges of each cube are checked
c     ***************************************************************** 

      npoly=0
      numcub=size-1

      do 5 d=1,numcub
         do 6 b=1,numcub
            do 7 a=1,numcub
               count=0
               do 8 i=1,3
                  do 9 j=0,1
                     do 10 k=0,1
                        if (i.eq.1) then
                           dh1=dh(a,b+j,d+k)
                           dh2=dh(a+1,b+j,d+k)
                           if (sgn(dh1-q)*sgn(dh2-q).eq.-1) then
                              count=count+1
                              temp(count,1)=fli(a,dh1,dh2,q)
                              temp(count,2)=b+j
                              temp(count,3)=d+k
                           endif
                        endif
                        if (i.eq.2) then
                           dh1=dh(a+j,b,d+k)
                           dh2=dh(a+j,b+1,d+k)
                           if (sgn(dh1-q)*sgn(dh2-q).eq.-1) then
                              count=count+1
                              temp(count,1)=a+j
                              temp(count,2)=fli(b,dh1,dh2,q)
                              temp(count,3)=d+k
                           endif
                        endif
                        if (i.eq.3) then
                           dh1=dh(a+j,b+k,d)
                           dh2=dh(a+j,b+k,d+1)
                           if (sgn(dh1-q)*sgn(dh2-q).eq.-1) then
                              count=count+1
                              temp(count,1)=a+j
                              temp(count,2)=b+k
                              temp(count,3)=fli(d,dh1,dh2,q)
                           endif
                        endif            
   10                continue
    9             continue
    8          continue                             

c     ***************************************************************** 
c     Store the points as 3-sided polygons
c     ***************************************************************** 

  100          if (count.gt.2) then
                  call angle(count,temp,tri)
                  npoly=npoly+1
                  do 11 i=1,3
                     do 12 j=1,3
                        c(npoly,i,j)=temp(tri(i),j)
   12                continue
   11             continue
                  count=count-1
                  goto 100
               endif

    7       continue
    6    continue
    5 continue

c     ***************************************************************** 
c     End of main loop
c     ***************************************************************** 


c     ***************************************************************** 
c     Store coordinates
c     ***************************************************************** 

      count=0
      do 13 i=1,npoly
         do 14 j=1,3
            count=count+1
            xcoord(count)=c(i,j,1)
            ycoord(count)=c(i,j,2)
            zcoord(count)=c(i,j,3)
   14    continue
   13 continue
      return
      end


c     ***************************************************************** 
c     Functions and subroutines
c     ***************************************************************** 


c     ***************************************************************** 
c     Function to linearly interpolate between point r and point (r+1)  
c     with density heights s and t respectively, where q is the         
c     required density height.
c     ***************************************************************** 

      double precision function fli(r,s,t,q)
      integer r
      double precision s, t, q
      fli=(abs(t-q)*r+abs(s-q)*(r+1))/abs(s-t)
      return
      end

c     ***************************************************************** 
c     Sgn function
c     ***************************************************************** 

      double precision function sgn(w)
      double precision w
      if (w.lt.0) sgn=-1
      if (w.gt.0) sgn=1
      if (w.eq.0) sgn=0
      return
      end

c     ***************************************************************** 
c     Subroutine to find one 3-sided polygon in a 'count'-sided one.    
c     This is made from the last pt and the two points collinear with   
c     it.  These 2 points are found by considering all combinations of  
c     angles between the last point and two of the other points.  The   
c     largest of these angles involves the two points we require.
c     ***************************************************************** 

      subroutine angle(count,temp,tri)
      double precision amax, temp(12,3), vlen(2), v(2,3), dotp, ang
      integer count, numa(2), tri(3), i, j, k
      amax=0.

      do 1 i=1,(count-2)
         do 2 j=(i+1),(count-1)

c     ***************************************************************** 
c     Calculate the two vectors from the last point to points i and j   
c     and obtain their lengths.
c     ***************************************************************** 

            vlen(1)=0.
            vlen(2)=0.
            do 3 k=1,3
               v(1,k)=temp(i,k)-temp(count,k)
               v(2,k)=temp(j,k)-temp(count,k)
               vlen(1)=vlen(1)+v(1,k)**2
               vlen(2)=vlen(2)+v(2,k)**2
   3        continue

c     ***************************************************************** 
c     Calculate the dot product, and hence obtain the angle between     
c     the two vectors.
c     ***************************************************************** 

            dotp=v(1,1)*v(2,1)+v(1,2)*v(2,2)+v(1,3)*v(2,3)
            ang=acos(dotp/sqrt(vlen(1)*vlen(2)))

c     ***************************************************************** 
c     If this is the largest angle found so far, then store it.         
c     ***************************************************************** 

            if (ang.gt.amax) then
               amax=ang
               numa(1)=i
               numa(2)=j
            endif
    2    continue
    1 continue

c     ***************************************************************** 
c     End of loop: go back and look at next pair of vectors (if any).   
c     ***************************************************************** 


c     ***************************************************************** 
c     Store the points which define the 3-sided polygon.
c     ***************************************************************** 

      tri(1)=numa(1)
      tri(2)=numa(2)
      tri(3)=count
      return
      end
