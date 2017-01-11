SUBROUTINE SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
!SUBROUTINE SupportDomain(gpos,x,ds,nv)
    !----------------------------------------------------------------------------
    ! This subroutine to determines nodes in the support domain of a Gauss point
    ! input--numnode: total number of field nodes;
    ! nx=2: for 2D problem;
    ! x(nx,numnode): coordinates of all field nodes;
    ! numgauss: number of Gauss points in a cell;
    ! gpos(2): x and y coordinate of a Gauss point;
    ! ds(nx,numnode): sizes of support domain;
    ! input and output-- ndex: when input ndex=0;
    ! when return ndex is the number of nodes in the support domain
    ! output--nv(ndex): No. of field nodes in the support domain
    !---------------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    dimension gpos(nx),x(nx,numnode),ds(nx,numnode),nv(numnode)


    !implicit none
!  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

!     REAL*8 :: gpos(:),ds(:,:)
!     INTEGER::x(:,:)
!
! REAL*8::nv(:)
! integer::numnode,ik,ndex
!REAL*8::dx,dy,eps

 !numnode=ubound(x,1)

  !  real:: gpos(nx),x(nx,numnode),ds(nx,numnode),nv(numnode)
    eps=1.E-16
    ndex=0
    nv=0

    !ds=1.0
    DO ik=1,numnode
    !Здесь временная подмена ds на 1 !!!!
        !dx=ds(1,ik)-dabs(gpos(1)-x(1,ik))
        dx=1-dabs(gpos(1)-x(1,ik))

        !dy=ds(2,ik)-dabs(gpos(2)-x(2,ik))
        dy=1-dabs(gpos(2)-x(2,ik))
        IF((dx.GE.eps).AND.(dy.GE.eps)) THEN
            ndex=ndex+1
            nv(ndex)=ik
        END IF
    ENDDO
    RETURN
END
