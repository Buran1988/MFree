SUBROUTINE Compute_Basis(gpos,gp,nx,mm)
    !------------------------------------------------------------------
    ! Compute basis functions and their derivatives
    ! Input:
    ! gpos - coordinates of the point of interest gpos(1)=x, gpos(2)=y;
    !  nx - problem dimension (2D -> nx=2);
    !  mm - number of monomials.
    ! Output-gp
    ! From 1 to 10 columns of gp: p,dpdx,dpdy,dpdxx,dpdxy,dpdyy,
    ! dpdxxx,dpdxxy,dpdxyy,dpdyyy
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION gpos(nx),gp(10,mm)


    gp=0.0

    gp(1,1)=1.0
    gp(1,2)=gpos(1)
    gp(1,3)=gpos(2)
    gp(1,4)=gpos(1)*gpos(1)
    gp(1,5)=gpos(1)*gpos(2)
    gp(1,6)=gpos(2)*gpos(2)
    gp(2,2)=1.0
    gp(2,4)=2.0*gpos(1)
    gp(2,5)=gpos(2)
    gp(3,3)=1.0
    gp(3,5)=gpos(1)
    gp(3,6)=2.0*gpos(2)
    gp(4,4)=2.0
    gp(5,5)=1.0
    gp(6,6)=2.0
    RETURN
END
