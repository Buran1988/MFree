subroutine ComputeGamma(mm,a,gam,komponent)

    IMPLICIT REAL*8 (A-H,O-Z)
    integer komponent
    DIMENSION gpos(nx),x(nx,numnode),nv(numnode)
    DIMENSION ds(nx,numnode),xv(nx,ndex)
    DIMENSION gp(10,mm),gam(mm,10),a(mm,mm,10)
    DIMENSION b(mm,ndex,10),c(mm),aa(mm,mm),phi(10,ndex)



        ! ************* Compute dgamdx
    c=0.0
    DO 30 in=1,mm
        DO 30 jn=1,mm
            if (komponent=2) then
                c(in)=c(in)+a(in,jn,2)*gam(jn,1)
            else if (komponent=2) then

            end if

30      CONTINUE
return

!        c=gp(2,:)-c
!        aa=a(:,:,1)
!        CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
!        gam(:,2)=c
!
!        ! ************* Compute dgamdy
!        c=0.0
!        DO 50 in=1,mm
!            DO 50 jn=1,mm
!                c(in)=c(in)+a(in,jn,3)*gam(jn,1)
!50          CONTINUE
!            c=gp(3,:)-c
!            aa=a(:,:,1)
!            CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
!            gam(:,3)=c
!
!            ! ************* Compute dgamdxx
!            DO 70 in=1,mm
!                c(in)=0.
!                DO 70 jn=1,mm
!                    c(in)=c(in)+a(in,jn,4)*gam(jn,1)+2.0*a(in,jn,2)*gam(jn,2)
!70              CONTINUE
!                DO 75 kn=1,mm
!                    c(kn)=gp(4,kn)-c(kn)
!75              CONTINUE
!                DO 80 i1=1,mm
!                    DO 80 j1=1,mm
!                        aa(i1,j1)=a(i1,j1,1)
!80                  CONTINUE
!                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
!                    DO 85 k1=1,mm
!                        gam(k1,4)=c(k1)
!85                  CONTINUE
end subroutine ComputeGamma
