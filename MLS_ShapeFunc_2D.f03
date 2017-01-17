SUBROUTINE MLS_ShapeFunc_2D(gpos,x,nv,ds,phi,nx,numnode,ndex,mm)
    !------------------------------------------------------------------
    ! Compute MLS shape functions and their derivatives
    ! Input--gpos,x,nv,ds,nx,numnode,ndex,mm
    ! Output--phi
    ! From 1 to 10 of the two dimension of phi denotes
    ! phi,dphix,dphiy,dphixx,dphixy,dphiyy
    ! dphidxxx,dphidxxy, dphidxyy, dphidyyy
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION gpos(nx),x(nx,numnode),nv(numnode)
    DIMENSION ds(nx,numnode),xv(nx,ndex)
    DIMENSION gp(10,mm),gam(mm,10),a(mm,mm,10)
    DIMENSION b(mm,ndex,10),c(mm),aa(mm,mm),phi(10,ndex)

    gp=0.0

    CALL Compute_Basis(gpos,gp,nx,mm)
    CALL Compute_AB(gpos,x,nv,ds,a,b,nx,numnode,ndex,mm)
    ep=1.0E-20

    c=gp(1,:)


    aa=a(:,:,1)

    gam=0.0

    ! ************* Compute gam
    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
21  FORMAT(1x,' gam kwji=',i2)

    gam(:,1)=c

    ! ************* Compute dgamdx
    c=0.0
    DO 30 in=1,mm
        DO 30 jn=1,mm
            c(in)=c(in)+a(in,jn,2)*gam(jn,1)
30      CONTINUE

        c=gp(2,:)-c
        aa=a(:,:,1)
        CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
        gam(:,2)=c

        ! ************* Compute dgamdy
        c=0.0
        DO 50 in=1,mm
            DO 50 jn=1,mm
                c(in)=c(in)+a(in,jn,3)*gam(jn,1)
50          CONTINUE
            c=gp(3,:)-c
            aa=a(:,:,1)
            CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
            gam(:,3)=c

            ! ************* Compute dgamdxx
            DO 70 in=1,mm
                c(in)=0.
                DO 70 jn=1,mm
                    c(in)=c(in)+a(in,jn,4)*gam(jn,1)+2.0*a(in,jn,2)*gam(jn,2)
70              CONTINUE
                DO 75 kn=1,mm
                    c(kn)=gp(4,kn)-c(kn)
75              CONTINUE
                DO 80 i1=1,mm
                    DO 80 j1=1,mm
                        aa(i1,j1)=a(i1,j1,1)
80                  CONTINUE
                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                    DO 85 k1=1,mm
                        gam(k1,4)=c(k1)
85                  CONTINUE
                    ! ************* Compute dgamdxy
                    DO 90 in=1,mm
                        c(in)=0.
                        DO 90 jn=1,mm
                            c(in)=c(in)+a(in,jn,5)*gam(jn,1)+a(in,jn,2)*gam(jn,3)+ &
                                a(in,jn,3)*gam(jn,2)
90                      CONTINUE
                        DO 95 kn=1,mm
                            c(kn)=gp(5,kn)-c(kn)
95                      CONTINUE
                        DO 100 i1=1,mm
                            DO 100 j1=1,mm
                                aa(i1,j1)=a(i1,j1,1)
100                         CONTINUE
                            CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                            DO 105 k1=1,mm
                                gam(k1,5)=c(k1)
105                         CONTINUE
                            ! ************* Compute dgamdyy
                            DO 110 in=1,mm
                                c(in)=0.
                                DO 110 jn=1,mm
                                    c(in)=c(in)+a(in,jn,6)*gam(jn,1)+2.0*a(in,jn,3)*gam(jn,3)
110                             CONTINUE
                                DO 115 kn=1,mm
                                    c(kn)=gp(6,kn)-c(kn)
115                             CONTINUE
                                DO 120 i1=1,mm
                                    DO 120 j1=1,mm
                                        aa(i1,j1)=a(i1,j1,1)
120                                 CONTINUE
                                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                                    DO 125 k1=1,mm
                                        gam(k1,6)=c(k1)
125                                 CONTINUE
                                    ! ************* Compute dgamdxxx
                                    DO in=1,mm
                                        c(in)=0.
                                        DO jn=1,mm
                                            c(in)=c(in)+a(in,jn,7)*gam(jn,1)+3*a(in,jn,4)*gam(jn,2)+ &
                                                3*a(in,jn,2)*gam(jn,4)
                                        ENDDO
                                    ENDDO
                                    DO kn=1,mm
                                        c(kn)=gp(7,kn)-c(kn)
                                    ENDDO
                                    DO i1=1,mm
                                        DO j1=1,mm
                                            aa(i1,j1)=a(i1,j1,1)
                                        ENDDO
                                    ENDDO
                                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                                    DO k1=1,mm
                                        gam(k1,7)=c(k1)
                                    ENDDO
                                    ! ************* Compute dgamdxxy
                                    DO in=1,mm
                                        c(in)=0.
                                        DO jn=1,mm
                                            c(in)=c(in)+a(in,jn,8)*gam(jn,1)+ &
                                                a(in,jn,4)*gam(jn,3)+2*a(in,jn,5)*gam(jn,2)+ &
                                                2*a(in,jn,2)*gam(jn,5)+a(in,jn,3)*gam(jn,4)
                                        ENDDO
                                    ENDDO
                                    DO kn=1,mm
                                        c(kn)=gp(8,kn)-c(kn)
                                    ENDDO
                                    DO i1=1,mm
                                        DO j1=1,mm
                                            aa(i1,j1)=a(i1,j1,1)
                                        ENDDO
                                    ENDDO
                                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                                    DO k1=1,mm
                                        gam(k1,8)=c(k1)
                                    ENDDO
                                    ! ************* Compute dgamdxyy
                                    DO in=1,mm
                                        c(in)=0.
                                        DO jn=1,mm
                                            c(in)=c(in)+a(in,jn,9)*gam(jn,1)+ &
                                                a(in,jn,6)*gam(jn,2)+2*a(in,jn,5)*gam(jn,3)+ &
                                                2*a(in,jn,3)*gam(jn,5)+a(in,jn,2)*gam(jn,6)
                                        ENDDO
                                    ENDDO
                                    DO kn=1,mm
                                        c(kn)=gp(9,kn)-c(kn)
                                    ENDDO
                                    DO i1=1,mm
                                        DO j1=1,mm
                                            aa(i1,j1)=a(i1,j1,1)
                                        ENDDO
                                    ENDDO
                                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                                    DO k1=1,mm
                                        gam(k1,9)=c(k1)
                                    ENDDO
                                    ! ************* Compute dgamdyyy
                                    DO in=1,mm
                                        c(in)=0.
                                        DO jn=1,mm
                                            c(in)=c(in)+a(in,jn,10)*gam(jn,1)+ &
                                                3*a(in,jn,6)*gam(jn,3)+3*a(in,jn,3)*gam(jn,6)
                                        ENDDO
                                    ENDDO
                                    DO kn=1,mm
                                        c(kn)=gp(10,kn)-c(kn)
                                    ENDDO
                                    DO i1=1,mm
                                        DO j1=1,mm
                                            aa(i1,j1)=a(i1,j1,1)
                                        ENDDO
                                    ENDDO
                                    CALL GaussEqSolver_Sym(mm,mm,aa,c,ep,kwji)
                                    DO k1=1,mm
                                        gam(k1,10)=c(k1)
                                    ENDDO
                                    !! ************* Compute Phi and their derivatives
                                    DO 130 iph=1,ndex
                                        DO iiii=1,10
                                            phi(iiii,iph)=0.0
                                        ENDDO
                                        DO 130 jph=1,mm
                                            phi(1,iph)=phi(1,iph)+gam(jph,1)*b(jph,iph,1)
                                            phi(2,iph)=phi(2,iph)+gam(jph,2)*b(jph,iph,1)+ &
                                                gam(jph,1)*b(jph,iph,2)
                                            phi(3,iph)=phi(3,iph)+gam(jph,3)*b(jph,iph,1)+ &
                                                gam(jph,1)*b(jph,iph,3)
                                            phi(4,iph)=phi(4,iph)+gam(jph,4)*b(jph,iph,1)+ &
                                                2.0*gam(jph,2)*b(jph,iph,2)+gam(jph,1)*b(jph,iph,4)
                                            phi(5,iph)=phi(5,iph)+gam(jph,5)*b(jph,iph,1)+ &
                                                gam(jph,2)*b(jph,iph,3)+gam(jph,3)*b(jph,iph,2)+ &
                                                gam(jph,1)*b(jph,iph,5)
                                            phi(6,iph)=phi(6,iph)+gam(jph,6)*b(jph,iph,1)+ &
                                                2.0*gam(jph,3)*b(jph,iph,3)+gam(jph,1)*b(jph,iph,6)
                                            phi(7,iph)=phi(7,iph)+gam(jph,7)*b(jph,iph,1)+ &
                                                3.0*gam(jph,4)*b(jph,iph,2)+3*gam(jph,2)*b(jph,iph,4)+ &
                                                gam(jph,1)*b(jph,iph,7)
                                            phi(8,iph)=phi(8,iph)+gam(jph,8)*b(jph,iph,1)+ &
                                                2.0*gam(jph,5)*b(jph,iph,2)+2*gam(jph,2)*b(jph,iph,5)+ &
                                                gam(jph,1)*b(jph,iph,8)+gam(jph,4)*b(jph,iph,3)+ &
                                                gam(jph,3)*b(jph,iph,4)
                                            phi(9,iph)=phi(9,iph)+gam(jph,9)*b(jph,iph,1)+ &
                                                2.0*gam(jph,5)*b(jph,iph,3)+2*gam(jph,3)*b(jph,iph,5)+ &
                                                gam(jph,1)*b(jph,iph,9)+gam(jph,6)*b(jph,iph,2)+ &
                                                gam(jph,2)*b(jph,iph,6)
                                            phi(10,iph)=phi(10,iph)+gam(jph,10)*b(jph,iph,1)+ &
                                                3.0*gam(jph,6)*b(jph,iph,3)+3*gam(jph,3)*b(jph,iph,6)+ &
                                                gam(jph,1)*b(jph,iph,10)
130                                     CONTINUE
                                        RETURN
                                    END
