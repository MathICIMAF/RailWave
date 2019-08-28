          subroutine EIGENSOLVERMXM (M,A,rows,cols,nev,ncv,w,v,
     &                               neigh,info,numIter,numOpS,
     &                               numOpBX,numReOR)

           integer           rows, nev, ncv, lworkl, info, ierr, i,
     &                       maxitr,ishfts,mode,rn,cn,cols,j,numIter,
     &                       numOpS,numOpBX,numReOR

           Double precision  tol

           character         bmat*1, which*2, trans*1

           logical           rvec

           integer           iparam(11),ipntr(15),
     &                       ipiv(rows),neigh(rows/3,cols/3-1)
           logical           select(ncv)

           Complex*16        M(rows,cols), A(rows,rows),
     &                       resid(rows),
     &                       v(rows,ncv), workd(3*rows),
     &                       workev(2*ncv), w(nev+1),
     &                       workl(3*ncv**2+5*ncv)

           Double precision  rwork(ncv)

           Complex*16        zero, one

           Double precision  dznrm2, dlapy2
           external          zcopy, mXx

           do 206 i = 1,11
               iparam(i) = 0.0D+0
206        continue
           do 207 i = 1,15
               ipntr(i) = 0.0D+0
207        continue
           do 209 i = 1,ncv
               rwork(i) = 0.0D+0
               select(i) = .FALSE.
209        continue
           do 210 i = 1,rows
               ipiv(i) = 0.0D+0
               resid(i) = dcmplx(rand(),rand())
               do 211 j = 1,ncv
                   v(i,j) = (0.0D+0, 0.0D+0)
211            continue
210        continue
           do 212 i = 1,3*rows
                workd(i) = (0.0D+0, 0.0D+0)
212        continue
           do 213 i = 1,2*ncv
                workev(i) = (0.0D+0, 0.0D+0)
213        continue
           do 214 i = 1,nev+1
                w(i) = (0.0D+0, 0.0D+0)
214        continue
           do 215 i = 1,3*ncv**2+5*ncv
                workl(i) = (0.0D+0, 0.0D+0)
215        continue

           zero = (0.0D+0, 0.0D+0)
           one = (1.0D+0, 0.0D+0)
           rn = rows/3
           cn = cols/3-1
           ldv = rows
           bmat  = 'G'
           which = 'LM'
           lworkl = 3*ncv**2+5*ncv
           tol    = 0.0
           ido    = 0
           info   = 1
           ishfts = 1
           maxitr = 300
           mode   = 3
           iparam(1) = ishfts
           iparam(3) = maxitr
           iparam(7) = mode
           trans='N'
           ierr = 0

           call zgetrf(rows,rows,A,rows,ipiv,ierr)

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

 10        continue
           call znaupd( ido, bmat, rows, which, nev, tol,
     &                    resid, ncv, v, ldv, iparam, ipntr,
     &                    workd, workl, lworkl, rwork, info)

           if (ido .eq. -1 ) then
               call mXx(M, rows, cols, workd(ipntr(1)), neigh,
     &                  workd(ipntr(2)))
               call zgetrs(trans, rows, 1, A, rows, ipiv,
     &                        workd(ipntr(2)), rows, ierr)
               go to 10
           else if ( ido .eq. 1) then
               call mXx(M, rows, cols, workd(ipntr(3)), neigh,
     &                  workd(ipntr(2)))
               call zgetrs(trans, rows, 1, A, rows, ipiv,
     &                        workd(ipntr(2)), rows, ierr)
              go to 10
           else if ( ido .eq. 2) then
                  call zcopy (rows, workd(ipntr(1)),1,
     &                      workd(ipntr(2)),1)
                go to 10
           end if
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
           if ( info .lt. 0 ) then
               print *, ' Error with _naupd, info = ',info
           else
               rvec = .true.               
               numIter = iparam(3)
               numOpS = iparam(9)
               numOpBX = iparam(10)
               numReOR = iparam(11)
               call zneupd(rvec, 'ALL', select, w, v, ldv, zero,
     &                     workev, bmat, rows, which, nev, tol,
     &                     resid, ncv, v, ldv, iparam, ipntr,
     &                     workd, workl, lworkl, rwork, ierr)
               if ( ierr .ne. 0 ) then
                  print *, ' Error with _neupd, info = ', ierr
               end if
           end if
 9000      continue
         end


!----------------------------------------------------------!
!         subroutine mXx (M, rows, cols, x, neighbors,w)
!             integer     rows, cols, rneigh, cneigh, i, nvi
!             integer     neighbors(rows/3, cols/3-1),
!     &                   indvi(cols/3-1)
!             Complex*16  M(rows, cols), w(rows), x(rows),
!     &                   sum1, sum2, sum3
!             integer     threei,threej
!
!             rneigh = rows/3
!             cneigh = cols/3-1
!             do 27 i = 1,rneigh
!                !Getting neighbors
!                nvi = 0
!                do 29 j = 1, cneigh
!                   if(neighbors(i,j) == -1) then
!                      go to 33
!                   end if
!                   indvi(j) = neighbors(i,j)
!                   nvi = nvi+1
!29              continue
!
!33              threei = 3*i
!                w(threei-2) = 0
!                w(threei-1) = 0
!                w(threei) = 0
!                sum1 = M(threei-2,1)*x(threei-2)
!                sum2 = M(threei-2,2)*x(threei-1)
!                sum3 = M(threei-2,3)*x(threei)
!                w(threei-2) = w(threei-2)+ sum1 + sum2 + sum3
!                sum1 = M(threei-1,1)*x(threei-2)
!                sum2 = M(threei-1,2)*x(threei-1)
!                sum3 = M(threei-1,3)*x(threei)
!                w(threei-1) = w(threei-1)+ sum1 + sum2 + sum3
!                sum1 = M(threei,1)*x(threei-2)
!                sum2 = M(threei,2)*x(threei-1)
!                sum3 = M(threei,3)*x(threei)
!               w(threei) = w(threei)+ sum1 + sum2 + sum3
!
!                do 28 k=1,nvi
!                   indvj = indvi(k)
!                   threej = 3*k
!                   sum1 = M(threei-2,threej+1)*x(3*indvj-2)
!                   sum2 = M(threei-2,threej+2)*x(3*indvj-1)
!                   sum3 = M(threei-2,threej+3)*x(3*indvj)
!                   w(threei-2) = w(threei-2)+ sum1 + sum2 + sum3
!                   sum1 = M(threei-1,threej+1)*x(3*indvj-2)
!                   sum2 = M(threei-1,threej+2)*x(3*indvj-1)
!                   sum3 = M(threei-1,threej+3)*x(3*indvj)
!                   w(threei-1) = w(threei-1)+ sum1 + sum2 + sum3
!                   sum1 = M(threei,threej+1)*x(3*indvj-2)
!                   sum2 = M(threei,threej+2)*x(3*indvj-1)
!                   sum3 = M(threei,threej+3)*x(3*indvj)
!                   w(threei) = w(threei)+ sum1 + sum2 + sum3
!28              continue
!27           continue
!             return
!       end
!----------------------------------------------------------!
       subroutine mXx (M, rows, cols, x, neighbors, w)

           integer     rows, cols, i, nvi, threei, threej

           integer     neighbors(rows/3, cols/3-1), indvi(cols/3-1)

           Complex*16  sum1, sum2, sum3

           Complex*16  M(rows, cols), w(rows), x(rows)


           do 27 i = 1, rows/3
               !Getting neighbors
               nvi = 0
               do 29 j = 1, cols/3-1
                   if(neighbors(i,j) == -1) then
                       go to 33
                   end if
                   indvi(j) = neighbors(i,j)
                   nvi = nvi+1
29             continue

33             threei = 3*i
               w(threei-2) = 0
               w(threei-1) = 0
               w(threei) = 0
               sum1 = M(threei-2,1)*x(threei-2)
               sum2 = M(threei-2,2)*x(threei-1)
               sum3 = M(threei-2,3)*x(threei)
               w(threei-2) = w(threei-2)+ sum1 + sum2 + sum3
               sum1 = M(threei-1,1)*x(threei-2)
               sum2 = M(threei-1,2)*x(threei-1)
               sum3 = M(threei-1,3)*x(threei)
               w(threei-1) = w(threei-1)+ sum1 + sum2 + sum3
               sum1 = M(threei,1)*x(threei-2)
               sum2 = M(threei,2)*x(threei-1)
               sum3 = M(threei,3)*x(threei)
               w(threei) = w(threei)+ sum1 + sum2 + sum3

               do 28 k=1,nvi
                   indvj = indvi(k)
                   threej = 3*k
                   sum1 = M(threei-2,threej+1)*x(3*indvj-2)
                   sum2 = M(threei-2,threej+2)*x(3*indvj-1)
                   sum3 = M(threei-2,threej+3)*x(3*indvj)
                   w(threei-2) = w(threei-2)+ sum1 + sum2 + sum3
                   sum1 = M(threei-1,threej+1)*x(3*indvj-2)
                   sum2 = M(threei-1,threej+2)*x(3*indvj-1)
                   sum3 = M(threei-1,threej+3)*x(3*indvj)
                   w(threei-1) = w(threei-1)+ sum1 + sum2 + sum3
                   sum1 = M(threei,threej+1)*x(3*indvj-2)
                   sum2 = M(threei,threej+2)*x(3*indvj-1)
                   sum3 = M(threei,threej+3)*x(3*indvj)
                   w(threei) = w(threei)+ sum1 + sum2 + sum3
28             continue
27         continue
       return
       end

