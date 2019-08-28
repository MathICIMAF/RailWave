       subroutine EIGENSOLVERMXMAXB(M, A, rows, cols, nev, ncv,
     &                              w, v, neigh, info)

            integer           rows, cols, nev, ncv, info,
     &                        i, j, lworkl, ierr

            integer           iparam(11), ipntr(15),
     &                        neigh(rows/3, cols/3-1)

            Double precision  tol, rwork(ncv)

            character         bmat*1, which*2, trans*1

            logical           rvec, select(ncv)

            complex           C(rows, rows)

            Complex*16        zero, sigma

            Complex*16        M(rows, cols), A(rows, cols),
     &                        resid(rows), workd(3*rows),
     &                        workev(2*ncv), workl(3*ncv**2+5*ncv),
     &                        w(nev+1), v(rows, ncv)

           external          zcopy, mXx, aXb

           zero = (0.0D+0, 0.0D+0)
           ldv = rows
           bmat  = 'G'
           which = 'LM'
           trans= 'N'
           sigma = zero
           tol    = 0.0
           ido    = 0
           info   = 1
           lworkl = 3*ncv**2+5*ncv
           iparam(1) = 1
           iparam(3) = 300
           iparam(7) = 3


           do 100 i = 1,ncv
               rwork(i) = 0.0D+0
               select(i) = .FALSE.
100        continue

           do 101 i = 1,rows
               resid(i) = dcmplx(rand(),rand())
101        continue


!          Llamado a la factorizacion incompleta
!           call icholesky(A, rows, cols, neigh, C)

c          %-------------------------------------------%
c          | M A I N   L O O P (Reverse communication) |
c          %-------------------------------------------%
10         continue
           call znaupd(ido, bmat, rows, which, nev, tol,
     &                 resid, ncv, v, ldv, iparam, ipntr,
     &                 workd, workl, lworkl, rwork, info)
           if (ido .eq. -1 ) then
               call mXx(M, rows, cols, workd(ipntr(1)), neigh,
     &                  workd(ipntr(2)))
               call aXb(A, rows, cols, workd(ipntr(2)), neigh, C)
               go to 10
           else if ( ido .eq. 1) then
               call mXx(M, rows, cols, workd(ipntr(3)), neigh,
     &                  workd(ipntr(2)))
               call aXb(A, rows, cols, workd(ipntr(2)), neigh, C)
               go to 10
           else if ( ido .eq. 2) then
               call zcopy(rows, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               go to 10
           end if
c          %-----------------------------------------%
c          | Either we have convergence, or there is |
c          | an error.                               |
c          %-----------------------------------------%
           if ( info .lt. 0 ) then
           else
c              %-------------------------------------------%
c              | No fatal errors occurred.                 |
c              | Post-Process using ZNEUPD .               |
c              |                                           |
c              | Computed eigenvalues may be extracted.    |
c              |                                           |
c              | Eigenvectors may also be computed now if  |
c              | desired.  (indicated by rvec = .true.)    |
c              %-------------------------------------------%
               rvec = .true.
               call zneupd(rvec, 'A', select, w, v, ldv, sigma,
     &                     workev, bmat, rows, which, nev, tol,
     &                     resid, ncv, v, ldv, iparam, ipntr,
     &                     workd, workl, lworkl, rwork, ierr)
           end if
9000       continue
       end


c---------------------------------------------------------------c
       subroutine aXb (A, rows, cols, x, nei, C)

           integer     rows, cols, j, inf, lwork,
     &                 revcom, colx, coly, colz,
     &                 matvec, precondLeft,  dotProd
           parameter   (matvec=1, precondLeft=2, dotProd=3)

           integer     irc(7), icntl(7), info(3),
     &                 nei(rows/3, cols/3-1)

           complex     C(rows, rows), b(rows)

           complex*16  work(8*rows+5), A(rows, cols), x(rows)

           real*8      cntl(5), rinfo(3)

           complex*16  ZERO
           parameter   (ZERO=(0.0d0,0.0d0))

           complex*16  zdotc
           external    zdotc
           intrinsic   sqrt


           lwork = 8*rows+5
           inf = 0

           do j = 1,rows
               work(j) = ZERO
           enddo

           do j = 1, rows
               work(j+rows) = x(j)
           enddo
!          *******************************************************
!          ** Initialize the control parameters to default value
!          *******************************************************
!          * setup the monitoring CG variables to default values
           call init_zcg(icntl, cntl)!             *
           icntl(5) = 1
!          * Define the bound for the stopping criterion
           cntl(1)  = 1.0d-7
!          * Define the stream for the convergence history
           icntl(3) = 1
!          * Define the maximum number of iterations
           icntl(6) = 300 !rows+1
!          * Ask for the estimation of the condition number (if icntl(1) == 1)
           icntl(7) = 0
!          * Preconditioner is provided
           icntl(4) = 0

           icntl(1) = 1

!          *****************************************
!          ** Reverse communication implementation
!          *****************************************
10         call drive_zcg(rows, rows, lwork, work, irc,
     &                    icntl, cntl, info, rinfo)
           revcom = irc(1)
           colx   = irc(2)
           coly   = irc(3)
           colz   = irc(4)

           if (revcom.eq.matvec) then
!            * perform the matrix vector product for the CG iteration
!!           *        work(colz) <-- A * work(colx)
             call mXx(A, rows, cols, work(colx), nei, work(colz))
             goto 10
           else if (revcom.eq.precondleft) then
!            *      work(colz) <-- M * work(colx
             call zcopy(rows, work(colx), 1, work(colz), 1)
!             call mXx(C, rows, cols, work(colx), nei, work(colz))
!             do j = 1, rows
               !b(j) = work(colz+j-1)
!             enddo

!!            Resuelve el sistema con factoriazacion incompleta de cholesky
!             call cpotrs('L', rows, 1, C, rows, b, rows, inf)
!             do j = 1,rows
!               work(colz+j-1) = b(j)
!             enddo
             goto 10
           else if (revcom.eq.dotProd) then
!            *      perform the scalar product for the CG iteration
!            *      work(colz) <-- work(colx) work(coly)
             work(colz)= zdotc(rows, work(colx),1, work(coly),1)
             goto 10
           endif

           do j = 1, rows
             x(j) = work(j)
           enddo
       end


!----------------------------------------------------------!
       subroutine icholesky(A, rows, cols, neigh, L)

           integer     rows, cols, k, i, j, ix3, jx3, neix3

           integer     neigh(rows/3, cols/3-1)

           complex*16  A(rows, cols)
           complex     L(rows, rows)


           n = rows/3
           do 77 i = 1, n
             ix3 = i*3
             do 7777 j = 1, cols/3-1
               if (neigh(i,j) .ne. -1) then
                 neix3 = neigh(i,j)*3
                 jx3 = j*3
                 if (ix3 .lt. neix3) then
                   L(neix3-2,ix3-2) = A(ix3-2,jx3+1)
                   L(neix3-1,ix3-2) = A(ix3-2,jx3+2)
                   L(neix3,ix3-2) = A(ix3-2,jx3+3)

                   L(neix3-2,ix3-1) = A(ix3-1,jx3+1)
                   L(neix3-1,ix3-1) = A(ix3-1,jx3+2)
                   L(neix3,ix3-1) = A(ix3-1,jx3+3)

                   L(neix3-2,ix3) = A(ix3,jx3+1)
                   L(neix3-1,ix3) = A(ix3,jx3+2)
                   L(neix3,ix3) = A(ix3,jx3+3)
                 end if
               end if
7777         continue
             L(ix3-2,ix3-2) = A(ix3-2,1)

             L(ix3-1,ix3-2) = A(ix3-1,1)
             L(ix3-1,ix3-1) = A(ix3-1,2)

             L(ix3,ix3-2) = A(ix3,1)
             L(ix3,ix3-1) = A(ix3,2)
             L(ix3,ix3) = A(ix3,3)
77         continue

           do 444 k = 1, rows
             L(k,k) = sqrt(real(L(k,k)))
             do 4444 i = k+1, rows
               if (L(i,k) .ne. 0) then
                 L(i,k) = L(i,k)/L(k,k)
               end if
4444         continue
             do 555 j = k+1, rows
               do 5555 i = j, rows
                 if (L(i,j) .ne. 0) then
                   L(i,j) = L(i,j) - L(i,k)*L(j,k)
                 end if
5555           continue
555          continue
444        continue
       end
