C*
C*  Copyright (C) CERFACS 2000
C*
C*  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED TO
C*  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALLING
C*  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD AND
C*  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS.
C*
C*  1 - This software program provided in source code format ("the " Source
C*  Code ") and any associated documentation (the " Documentation ") are
C*  licensed, not sold, to you.
C*
C*  2 - CERFACS grants you a personal, non-exclusive, non-transferable and
C*  royalty-free right to use, copy or modify the Source Code and
C*  Documentation, provided that you agree to comply with the terms and
C*  restrictions of this agreement. You may modify the Source Code and
C*  Documentation to make source code derivative works, object code
C*  derivative works and/or documentation derivative Works (called "
C*  Derivative Works "). The Source Code, Documentation and Derivative
C*  Works (called " Licensed Software ") may be used by you for personal
C*  and non-commercial use only. " non-commercial use " means uses that are
C*  not or will not result in the sale, lease or rental of the Licensed
C*  Software and/or the use of the Licensed Software in any commercial
C*  product or service. CERFACS reserves all rights not expressly granted
C*  to you. No other licenses are granted or implied.
C*
C*  3 - The Source Code and Documentation are and will remain the sole
C*  property of CERFACS. The Source Code and Documentation are copyrighted
C*  works. You agree to treat any modification or derivative work of the
C*  Licensed Software as if it were part of the Licensed Software itself.
C*  In return for this license, you grant CERFACS a non-exclusive perpetual
C*  paid-up royalty-free license to make, sell, have made, copy, distribute
C*  and make derivative works of any modification or derivative work you
C*  make of the Licensed Software.
C*
C*  4- The licensee shall acknowledge the contribution of the Source Code
C*  in any publication of material dependent upon the use of the Source
C*  Code. The licensee shall use reasonable endeavours to send to CERFACS a
C*  copy of each such publication.
C*
C*  5- CERFACS has no obligation to support the Licensed Software it is
C*  providing under this license.
C*
C*  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO
C*  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE,
C*  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF
C*  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF
C*  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD
C*  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL NOT
C*  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR ANY
C*  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR
C*  USE OF THE LICENSED SOFTWARE.
C*
C*  6- For information regarding a commercial license for the Source Code
C*  and Documentation, please contact Mrs Campassens (campasse@cerfacs.fr)
C*
C*  7- This license is effective until terminated. You may terminate this
C*  license at any time by destroying the Licensed Software.
C*
C*    I agree all the terms and conditions of the above license agreement
C*
        subroutine zcg(n,b,x,r,z,p,q,D,E,aux,
     &                 irc,icntl,cntl,info,rinfo)
*
*
*  Purpose
*  =======
*  dcg solves the linear system Ax = b using the
*  Conjugate Gradient iterative method
*
* When preconditioning is used we solve :
*     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
*     x = M_2^{-1} y
*
*   Convergence test based on the normwise backward error for
*  the preconditioned system
*
* Written : June 2000
* Authors : Luc Giraud, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : August 2002
* Authors : Luc Giraud 
*             Parallel Algorithms - CERFACS
* Purpose : Remove the implicit type conversion from complex to real
*           in the complex module when computing alpha and rho.
*
* Updated : January 10, 2003
* Authors : Luc Giraud
*             Parallel Algorithms - CERFACS
* Purpose : Declare the variable bea as SAVE
*           (bug reported by an external user)
*
*  Arguments
*  =========
*
*  n       (input) INTEGER.
*           On entry, the dimension of the problem.
*           Unchanged on exit.
*
*
*  b        (input)  real*8/complex*16
*           Right hand side of the linear system.
*
*  x        (output) real*8/complex*16
*           Computed solution of the linear system.
*
*  r        (workspace)   real*8/complex*16
*           Vector used to store the residual
*
*  z        (workspace)  real*8/complex*16
*           Vector used to store the preconditioned residual
*
*  p        (workspace) real*8/complex*16
*           Vector used to store the descent direction
*
*  q        (workspace) real*8/complex*16
*           Vector used to store the gradient
*
*  D        (workspace) real*8/complex*16
*           Vector used to compute the approximation of the eigenvalues
*
*  E        (workspace) real*8/complex*16
*           Vector used to compute the approximation of the eigenvalues
*
*  aux      (workspace)  real*8/complex*16
*            Scalar used to store the ddot product
*
*  irc     (input/output) INTEGER array. length 3
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*
* icntl    (input) INTEGER array. length 7
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - error, default set in Init
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(6) : maximum number of iterations
*            icntl(7) : 0 - no estimation of the condition number
*                       1 - estimation of the condition number
*                           (actually estimations of the smallest and largest eigenvalues
*                            which ratio defines the condition number)
*
* cntl     (input) real*8 array, length 3
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*
* info     (output) INTEGER array, length 2
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - lwork too small
*                      -3 - convergence not achieved after icntl(7) iterations
*                      -4 - precondition type not set by user
*            info(2) : if info(1)=0 - number of iteration to converge
*                      if info(1)=-2 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
* rinfo    (output) real*8 array, length 3
*            if info(1)=0 
*              rinfo(1) : backward error for the unpreconditioned system
*            if icntl(7)=1 
*              rinfo(2) : largest estimated eigenvalue
*              rinfo(3) : smallest estimated eigenvalue
*              The ratio rinfo(2)/rinfo(3) is an estimation of the condition number
*
* Input variables
* ---------------
        integer  n, icntl(*)
        complex*16   b(*)
        real*8   cntl(*)
*
* Output variables
* ----------------
       integer  info(*)
       real*8   rinfo(*)
*
* Input/Output variables
* ----------------------
       integer  irc(*)
       complex*16 x(*), r(*), z(*), p(*), q(*), aux
       real*8 D(*), E(*)
*
* Local variables
* ---------------
       integer   j, iterMax, initGuess, iter, i
       integer   xptr, bptr, rptr, zptr, pptr, qptr, auxptr
       integer   Eptr, Dptr
       integer   typePrec, leftPrec, noPrec
       integer   iwarn, ihist, iest
       real*8    beta, bn, sA, sb, trueNormRes
       real*8    alpha, rho, rho1, dnormres, dnormx
       complex*16 zalpha
       real*8    alpha_old, be, bea
       parameter (noPrec = 0, leftPrec = 1)
*
*
       complex*16 ZERO, ONE
       parameter (ZERO = (0.0d0, 0.0d0), ONE = (1.0d0, 0.0d0))
       real*8 DZRO,DONE
       parameter (DZRO = 0.0d0, DONE = 1.0d0)
*
* Reverse communication variables
* -------------------------------
       integer retlbl
       DATA retlbl /0/
       integer matvec, precondLeft, prosca
       parameter(matvec=1, precondLeft=2, prosca=3)
*
* Saved variables
* ---------------
      save beta, bn, dnormres, retlbl, j
      save sA, sb, be, iter, bea
      save alpha, rho, rho1, dnormx
      save alpha_old
*
* Intrinsic function
* ------------------
        intrinsic dsqrt, dreal, dconjg
*
*       Executable statements
*
       iwarn      = icntl(2)
       ihist      = icntl(3)
       typePrec   = icntl(4)
       initGuess  = icntl(5)
       iterMax    = icntl(6)
       iest       = icntl(7)
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + n
       rptr     = bptr + n
       zptr     = rptr + n
       pptr     = zptr + n
       qptr     = pptr + n
       auxptr   = qptr + n
       Dptr     = auxptr + 1
       Eptr     = Dptr + itermax + 1
*
      if (retlbl.ne.0) then
          if (retlbl.eq.5) then
            goto 5
          else if (retlbl.eq.11) then
            goto 11
          else if (retlbl.eq.16) then
            goto 16
          else if (retlbl.eq.18) then
            goto 18
          else if (retlbl.eq.21) then
            goto 21
          else if (retlbl.eq.26) then
            goto 26
          else if (retlbl.eq.31) then
            goto 31
          else if (retlbl.eq.38) then
            goto 38
          else if (retlbl.eq.41) then
            goto 41
          else if (retlbl.eq.43) then
            goto 43
          endif
        endif
*
*
* intialization of various variables
*
        iter     = 0
*
        if (initGuess.eq.0) then
          do j=1,n
            x(j) = ZERO
          enddo
        endif
*
*        bn = snrm2(n,b,1)
*
        irc(1) = prosca
        irc(2) = bptr
        irc(3) = bptr
        irc(4) = auxptr
        retlbl = 5
        return
 5      continue
        bn = dsqrt(dreal(aux))
*
        if (bn.eq.DZRO) then
          do j=1,n
            x(j) = ZERO
          enddo  
          if (iwarn.ne.0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING : null right-hand side'
            write(iwarn,*) '           solution set to zero'
            write(iwarn,*)
          endif
          info(1)  = 0
          info(2)  = 0
          rinfo(1) = ZERO
          irc(1)   = 0
          return
        endif
*
* Compute the scaling factor for the backward error on the 
*  unpreconditioned sytem
*
       sA       = cntl(2)
       sb       = cntl(3)
       if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
         sb = bn
       endif
*
* Compute the first residual
*           Y = AX : r0 <-- A x
*
* The residual is computed only if the initial guess is not zero
*
       if (initGuess.ne.0) then
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = rptr
         retlbl = 11
         return
       endif
 11    continue
       if (initGuess.ne.0) then
         do j=1,n
           r(j) = b(j)-r(j)
         enddo
       else
         call zcopy(n,b,1,r,1)
       endif 
*
       iter = 0
*
*       loop : dcg iteration
*
*       REPEAT
 7      continue
        iter = iter + 1
*
*
*
* Compute the preconditioned residual if necessary
*      M_1Y = X : z <-- M_1^{-1} r
*
       if (typePrec.eq.noPrec) then
         call zcopy(n,r,1,z,1)
       else
         irc(1) = precondLeft
         irc(2) = rptr
         irc(4) = zptr
         retlbl = 16
         return
       endif
 16    continue
*
*
* rho = r^H * z
*
        irc(1) = prosca
        irc(2) = zptr
        irc(3) = rptr
        irc(4) = auxptr
        irc(5) = 1
        retlbl = 18
        return
 18     continue
        rho = dreal(aux)
*
        if ((rho.eq.DZRO).and.(iter.eq.1)) then
         if (iwarn.ne.0) then
              write(iwarn,*)
              write(iwarn,*) ' WARNING : You started from the solution '
              write(iwarn,*)
          endif
          irc(1)  = 0
          info(1) = 0
          info(2) = iter
          retlbl  = 0
          rinfo(1) = DZRO
          return
        endif
        if (iter.eq.1) then
          call zcopy(n,z,1,p,1)
        else
          beta = rho/rho1
          do i=1,n
            p(i) = z(i) + beta*p(i)
          enddo
        endif
*
*           Y = AX : q <-- A p
*
         irc(1) = matvec
         irc(2) = pptr
         irc(4) = qptr
         retlbl = 21
         return
 21      continue
*
* aux = p^H * q
*
        irc(1) = prosca
        irc(2) = pptr
        irc(3) = qptr
        irc(4) = auxptr
        irc(5) = 1
        retlbl = 26
        return
 26     continue
*
        if (iest.eq.1) alpha_old = alpha
*
        alpha = rho / dreal(aux)
*
        zalpha = alpha*DONE
        call zaxpy(n,zalpha,p,1,x,1)
        call zaxpy(n,-zalpha,q,1,r,1)
*
        if (iest.eq.1) then
         if (iter.gt.1) then
           D(iter-1) = D(iter-1) + 1.0d0/alpha_old
           D(iter) = beta/alpha_old
           E(iter-1)   = -sqrt(beta)/alpha_old
         else
           D(1) = DZRO
         endif
        endif
*
        rho1 = rho
*
* Check the convergence
*
* aux = r^H * r
*
        irc(1) = prosca
        irc(2) = rptr
        irc(3) = rptr
        irc(4) = auxptr
        irc(5) = 1
        retlbl = 31
        return
 31     continue
        dnormres = dsqrt(dreal(aux))
*
        if (sA.ne.DZRO) then
*
*         dnormx = snrm2(n,xCurrent,1)
*
          irc(1) = prosca
          irc(2) = xptr
          irc(3) = xptr
          irc(4) = auxptr
          irc(5) = 1
          retlbl = 38
          return
       else
          dnormx    = ONE
       endif
 38    continue
       if (sA.ne.DZRO) then
         dnormx = dsqrt(dreal(aux))
       else
         dnormx = ONE
       endif
*
       bea = dnormres/(sA*dnormx+sb)
*
* Check the convergence based on the Backward error for the
* unpreconditioned system
       if ((bea.le.cntl(1)).or.(iter.ge.iterMax)) then
*
* Use q vector to compute the residual to avoid "disturbing" the
* residual orthogonality if the iteration should continue
* q = A * x
*
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = qptr
         retlbl = 41
         return
       endif
 41    continue
*       
       if ((bea.le.cntl(1)).or.(iter.ge.iterMax)) then
* q = b - A * x
         do j=1,n
           q(j) = b(j) - q(j)
         enddo
*
* Compute the norm of the unpreconditioned residual
*
*        trueNormRes = snrm2(n,r,1)
*
         irc(1) = prosca
         irc(2) = qptr
         irc(3) = qptr
         irc(4) = auxptr
         retlbl = 43
         return
       endif
 43    continue
*
       if ((bea.le.cntl(1)).or.(iter.ge.iterMax)) then
         trueNormRes = dsqrt(dreal(aux))
         be          = trueNormRes/(sA*dnormx+sb)
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,'(I5,11x,E8.2,$)') iter,bea
           write(ihist,'(7x,E8.2)') be
         endif
         if ((be.le.cntl(1)).or.(iter.ge.iterMax)) then
* Return the backward errors
           rinfo(1) = be
*
* Compute the approximation of the condition number
           if ((iest.eq.1).and.(iter.gt.2)) then
             call dstev('N',iter-2,D,E,p,1,aux,i)
             rinfo(2) = D(1)
             rinfo(3) = D(iter-2)
           endif
*
           if ((iter.le.iterMax).and.(be.le.cntl(1))) then
             info(1) = 0
             if (ihist.ne.0) then
                write(ihist,*)
                write(ihist,'(A20)') 'Convergence achieved'
                write(ihist,'(A29,$)') 'B.E. on the unpreconditioned '
                write(ihist,'(A8,E8.2)') 'system: ', rinfo(1)
             endif
           else
             if (iwarn.ne.0) then
                write(iwarn,*)
                write(iwarn,*) ' WARNING : No convergence after '
                write(iwarn,*) iter,' iterations '
                write(iwarn,*)
             endif
             if (ihist.ne.0) then
                write(ihist,*)
                write(ihist,*) ' WARNING : No convergence after '
                write(ihist,*) iter,' iterations '
                write(ihist,*)
                write(ihist,'(A29,$)') 'B.E. on the unpreconditioned '
                write(ihist,'(A8,E8.2)') 'system: ', rinfo(1)
             endif
             info(1) = -3
            endif
            info(2) = iter
            if (ihist.ne.0) then
                 write(ihist,'(A10,I2)') 'info(1) = ',info(1)
                 write(ihist,'(A32,I5)')
     &                'Number of iterations (info(2)): ',info(2)
            endif
            irc(1)  = 0
            retlbl  = 0
            return
         endif
       else
* Save the backward error on a file if convergence history requested
          if (ihist.ne.0) then
            write(ihist,'(I5,11x,E8.2,$)') iter,bea
            write(ihist,'(9x,A2)') '--'
          endif
       endif
*
       goto 7
*
       end
        subroutine drive_zcg(n,nloc,lwork,work,
     &                         irc,icntl,cntl,info,rinfo)
*
*  Purpose
*  =======
*    DRIVE_DCG is the driver routine for solving the linear system
*  Ax = b where A is symmetric positive definite using the Conjugate 
*  Gradient iterative method with preconditioning.
*  This solver is implemented with a reverse communication scheme: control
*  is returned to the user for computing the (preconditioned) 
*  matrix-vector product or the dot products.
*  See the README file for an example of use.
*
*
*
*  Arguments
*  =========
*
*  n      (input) INTEGER.
*          On entry, the dimension of the problem.
*          Unchanged on exit.
*
*  nloc   (input) INTEGER.
*          In a parallel distributed envirionment, this corresponds
*          to the size of the subset of entries of the right hand side
*          and solution allocated to the calling process.
*          Unchanged on exit.
*
*
*  lwork   (input) INTEGER
*          size of the workspace
*          lwork >=  6*n+1, if icntl(7) == 0; i.e. No cond. number estimation
*          lwork >=  6*n+1+2*(itermax+1) else
*
*  work    (workspace) real*8 array, length lwork
*          work contains the required vector and matrices stored in the 
*          following order :
*            x  (n,1)       : computed solution.
*            b  (n,1)       : right hand side.
*            r  (n,1)       : vector workspace.
*            z  (n,1)       : vector workspace.
*            p  (n,1)       : vector workspace.
*            q  (n,1)       : vector workspace.
*            D(iterMax,1)   : diagonal of the tridiagonal symmetric Lanczos matrix.
*            E(iterMax,1)   : subdiagonal of the tridiagonal symmetric Lanczos matrix.
*            aux(1,1)       : dot product result.
*
*  irc     (input/output) INTEGER array. length 4
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*                     1 - perform the matrix vector product
*                     2 - perform the preconditioning
*                     3 - perform the scalar product
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*
* icntl    (input) INTEGER array. length 7
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - error, default set in Init
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(6) : maximum number of iterations
*            icntl(7) : 0 - no estimation of the condition number
*                       1 - estimation of the condition number
*                           (actually estimations of the smallest and largest eigenvalues
*                            which ratio defines the condition number)
*
*
* cntl     (input) real*8 array, length 3
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*
* info     (output) INTEGER array, length 3
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - lwork too small
*                      -3 - convergence not achieved after icntl(7) iterations
*                      -4 - precondition type not set by user
*            info(2) : if info(1)=0 - number of iteration to converge
*                      if info(1)=-2 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
* rinfo    (output) real*8 array, length 3
*            if info(1)=0 
*              rinfo(1) : backward error for the unpreconditioned system
*            if icntl(7)=1
*              rinfo(2) : largest estimated eigenvalue
*              rinfo(3) : smallest estimated eigenvalue
*              The ratio rinfo(2)/rinfo(3) is an estimation of the condition number
*
* Input variables
* ---------------
       integer n, nloc, lwork, icntl(*)
       real*8   cntl(*)
       real*8   sA, sb
* Output variables
* ----------------
       integer  info(*)
       real*8    rinfo(*)
* Input/Output variables
* ----------------------
       integer  irc(*)
       complex*16 work(*)
* Local variables
* ---------------
       integer xptr, bptr, rptr, zptr, pptr, qptr, Dptr, Eptr
       integer auxptr
       integer sizeWrk
       integer iwarn, ierr, ihist, itermax
       real*8 DZRO
       parameter (DZRO = 0.0d0)
*
       integer icheck
       DATA icheck /0/
       save icheck
*
*       Executable statements :
*
       ierr  = icntl(1)
       iwarn = icntl(2)
       ihist = icntl(3)
       itermax = icntl(6)
*
*
*
       if (icheck.eq.0) then
* Check the value of the arguments
         if (n.lt.1) then
            write(ierr,*)
            write(ierr,*)' ERROR : N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
         endif
         if ((icntl(4).ne.0).and.(icntl(4).ne.1)) then
            write(ierr,*)
            write(ierr,*)' ERROR : Undefined preconditioner '
            write(ierr,*)
            info(1) = -4
            irc(1)  = 0
            return
         endif
*
         if ((icntl(5).ne.0).and.(icntl(5).ne.1)) then
           icntl(5) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING : Undefined initial guess '
             write(iwarn,*) '           default x0 = 0 '
             write(iwarn,*)
           endif
         endif
         if (icntl(6).le.0) then
           icntl(6) = n
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,'(A31,$)') ' WARNING : Negative max number '
             write(iwarn,*) 'of iterations'
             write(iwarn,*) '           default ',n
             write(iwarn,*)
           endif
         endif
* Compute the size of the wokspace depending on the request to compute
* or not compute the estimated condition number.
       if (icntl(7).ne.0) then
         sizeWrk  = 6*nloc+1+2*(itermax+1)
       else
         sizeWrk  = 6*nloc+1
       endif
*
* Check if workspace is big enough for the problem.
*
         if (lwork.lt.sizeWrk) then
            write(ierr,*)
            write(ierr,*)' ERROR : Not enough space for the problem'
            write(ierr,*)
            info(1) = -2
            info(2) = sizeWrk
            irc(1)  = 0
            return
         endif
*
         info(3) = sizeWrk
         icheck = 1
*
* save the parameters in the history file
*
         if (ihist.ne.0) then
           write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR CG'
           write(ihist,*)
           write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr 
           if (iwarn.eq.0) then
             write(ihist,'(A27)') 'Warnings are not displayed:'
           else
             write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ',
     &                               iwarn
           endif 
           write(ihist,'(A13,I7)') 'Matrix size: ',n
           if (icntl(4).eq.0) then
             write(ihist,'(A18)') 'No preconditioning'
           elseif (icntl(4).eq.1) then
             write(ihist,'(A20)') 'Left preconditioning'
           endif
           if (icntl(5).eq.0) then
             write(ihist,'(A29)') 'Default initial guess x_0 = 0'
           else
             write(ihist,'(A27)') 'User supplied initial guess'
           endif
           write(ihist,'(A30,I5)') 'Maximum number of iterations: ',
     &                              icntl(6)
           write(ihist,'(A27,E8.2)') 'Tolerance for convergence: ', 
     &                                cntl(1) 
* 
           write(ihist,'(A53)') 
     &       'Backward error on the unpreconditioned system Ax = b:'
           sA       = cntl(2)
           sb       = cntl(3)
           if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
             write(ihist,'(A39)') 
     &       '    the residual is normalised by ||b||'
           else
             write(ihist,'(A34)') '    the residual is normalised by '
             write(ihist,'(A8,E8.2,$)') '        ', sA
             write(ihist,'(A11,E8.2)') ' * ||x|| + ', sb
           endif
*
           write(ihist,'(A31,I7)') 'Optimal size for the workspace:',
     &                              info(3)
           write(ihist,*) 
           write(ihist,'(A32,$)') 'Convergence history: b.e. on the'
           write(ihist,'(A24)') ' unpreconditioned system'
           write(ihist,'(A11,$)') ' Iteration '
           write(ihist,'(A27)') '  Approx. b.e.    True b.e.'
         endif
*
       endif
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + nloc
       rptr     = bptr + nloc
       zptr     = rptr + nloc
       pptr     = zptr + nloc
       qptr     = pptr + nloc
       auxptr   = qptr + nloc
       Dptr     = auxptr + 1
       Eptr     = Dptr + itermax + 1
*
*
       call zcg(nloc,work(bptr),work(xptr),
     &            work(rptr),work(zptr),work(pptr),work(qptr),
     &            work(Dptr), work(Eptr),
     &            work(auxptr), irc,icntl,cntl,info,rinfo)
*
       if (irc(1).eq.0) then
         icheck = 0
       endif
*
       return
       end
        subroutine init_zcg(icntl,cntl)
*
*  Purpose
*  =======
*    Set default values for the parameters defining the characteristics
* of the Conjugate Gradient algorithm.
*
*  Arguments
*  =========
*
* icntl    (output) INTEGER array. length 6
*            icntl(1) : stdout for error messages
*                       default value 6
*            icntl(2) : stdout for warnings
*                       default value 6
*            icntl(3) : stdout for convergence history
*                       default value 0
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       default value 0
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*                       default value 0
*            icntl(6) : maximum number of iterations
*                       default value -1, requiring the user to initialise it
*            icntl(7) : 0 - no estimation of the condition number
*                       1 - estimation of the condition number
*			    (actually estimations of the smallest and largest eigenvalues
*			     which ratio defines the condition number)
*                       default value 0
*
* cntl     (output) real*8 array, length 3
*            cntl(1) : tolerance for convergence
*                       default value 10-6
*            cntl(2) : scaling factor for normwise perturbation on A
*                       default value 0.0
*            cntl(3) : scaling factor for normwise perturbation on b
*                       default value 0.0
*
* Output variables
* ----------------
       integer icntl(*)
       real*8   cntl(*)
*
       icntl(1) = 6
       icntl(2) = 6
       icntl(3) = 0
       icntl(4) = 2
       icntl(5) = 0
       icntl(6) = -1
       icntl(7) = 0
* 
       cntl(1) = 1.0 d -6
       cntl(2) = 0.0 d 0
       cntl(3) = 0.0 d 0
*
       return
       end
