       subroutine EIGENSOLVERQ (M,A,n,nev,ncv,w,v,info)

            integer           n, nev, ncv, lworkl, info, ierr, j,
     &                        nconv, maxitr, ishfts, mode

            Double precision  tol

            character         bmat*1, which*2, trans*1

            logical           rvec

            integer           iparam(11),ipntr(15),ipiv(n)

            logical           select(ncv)

            Complex*16        M(n,n), A(n,n), resid(n),
     &                        v(n,ncv), workd(3*n),
     &                        workev(2*ncv), w(nev+1),
     &                        workl(3*ncv**2+5*ncv)

            Double precision  rwork(ncv)

            Complex*16        zero, one

           Double precision  dznrm2, dlapy2
           external          zaxpy, zcopy, dznrm2, dlapy2

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
           do 210 i = 1,n
               ipiv(i) = 0.0D+0
               resid(i) = dcmplx(rand(),rand())
               do 211 j = 1,ncv
                   v(i,j) = (0.0D+0, 0.0D+0)
211            continue
210        continue
           do 212 i = 1,3*n
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
           ldv = n
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

c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
          call zgetrf(n, n, A, n, ipiv, ierr)

 10        continue
           call znaupd( ido, bmat, n, which, nev, tol,
     &                    resid, ncv, v, ldv, iparam, ipntr,
     &                    workd, workl, lworkl, rwork, info)

           if (ido .eq. -1 ) then
               call zgemv('N', n, n, one, M, n,
     &                    workd(ipntr(1)), 1, zero, workd(ipntr(2)), 1)
               call zgetrs(trans, n, 1, A, n, ipiv,
     &                        workd(ipntr(2)), n, ierr)
               go to 10
           else if ( ido .eq. 1) then
               call zgemv('N', n, n, one, M, n,
     &                    workd(ipntr(3)), 1, zero, workd(ipntr(2)), 1)
               call zgetrs(trans, n, 1, A, n, ipiv,
     &                        workd(ipntr(2)), n, ierr)
               go to 10
           else if ( ido .eq. 2) then
                  call zcopy (n, workd(ipntr(1)),1,
     &                      workd(ipntr(2)),1)
                go to 10
           end if
c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%
           if ( info .lt. 0 ) then
               print *, ' Error with _naupd, info = ',info
           else
               rvec = .true.
               call zneupd(rvec, 'A', select, w, v, ldv, zero,
     &                     workev, bmat, n, which, nev, tol,
     &                     resid, ncv, v, ldv, iparam, ipntr,
     &                     workd, workl, lworkl, rwork, ierr)
               if ( ierr .ne. 0 ) then
                  print *, ' Error with _neupd, info = ', ierr
               end if
            end if
 9000       continue
         end
