SUBROUTINE RPIM_ShapeFunc_2D(gpos,x,nv,phi,nx,numnode,ndex,&
    alfc,dc,q,nRBF, mbasis)

    !------------------------------------------------------------------
    ! Compute RPIM shape functions and their derivatives
    ! Input--gpos,x,nv,ds,alfc,dc,q,nx,numnode,ndex,mm,nRBF,nbasis
    ! nRBF=1: MQ; 2: EXP; 3: TSP
    ! Output--phi
    ! From 1 to 10 of the two dimension of phi denotes
    ! phi,dphix,dphiy,dphixx,dphixy,dphiyy
    ! dphidxxx,dphidxxy, dphidxyy, dphidyyy
    !------------------------------------------------------------------

    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION gpos(nx),x(nx,numnode),nv(ndex),rk(ndex+mbasis)
    DIMENSION phi(10,ndex),xv(nx,ndex),rr(10,ndex+mbasis)
    DIMENSION a(ndex+mbasis,ndex+mbasis),g0(ndex+mbasis,ndex+mbasis)

    IF(nrbf.EQ.1) THEN
        rc=alfc*dc ! For MQ;
    ENDIF
    IF(nrbf.EQ.2) THEN
        q=alfc/dc/dc ! For EXP;
    ENDIF
    ep=1.e-20 !required tolerance
    mg=ndex+mbasis

    g0=0.

    DO i=1,ndex !ndex - number of field nodes used in the support domain
        nn=nv(i)
        xv(1,i)=x(1,nn)
        xv(2,i)=x(2,nn)
    ENDDO
    ! ****************** Assemble the matrix of G0
    DO i=1,ndex
        nn=nv(i)
        CALL Compute_RadialBasis(x(1,nn),x(2,nn),xv,rr,ndex,rc,q,nRBF,mbasis)

            g0(i,:)=rr(1,:)

        IF(mbasis.GT.0) THEN
            g0(i,ndex+1)=1.
            g0(i,ndex+2)=x(1,nn)
            g0(i,ndex+3)=x(2,nn)
            g0(ndex+1,i)=1.
            g0(ndex+2,i)=x(1,nn)
            g0(ndex+3,i)=x(2,nn)
        ENDIF
    ENDDO
    ! ****************** Solve linear equation to obtain shape function phi

    a=g0

    CALL Compute_RadialBasis(gpos(1),gpos(2),xv,rr,ndex,rc,q,nRBF,mbasis)

    rk=rr(1,:)
    CALL GaussEqSolver_Sym(mg,mg,a,rk,ep,kwji)
    IF(kwji==1) THEN
        WRITE(*,*)'Fail...'
       ! PAUSE
    ENDIF

    phi(1,:)=rk

    ! ****************** Solve linear equation to obtain
    !    dphidx,dphidy,dphidxx,dphidxy,dphidyy
    ! dphidxxx,dphidxxy, dphidxyy, dphidyyy
    do i=2, 6
        a=g0

        rk=rr(i,:)

        CALL GaussEqSolver_Sym(mg,mg,a,rk,ep,kwji)

        phi(i,:)=rk

    end do

    RETURN
END
