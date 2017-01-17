SUBROUTINE Compute_AB(gpos,x,nv,ds,a,b,nx,numnode,ndex,mm)
    !------------------------------------------------------------------
    ! Compute A matrix and B matrix and their derivatives
    ! input--gpos,x,nv,dm,nx,numnode,ndex,mm
    ! output--a,b
    ! From 1 to 10 of the third dimension of a denotes
    ! a,dax,day,daxx,daxy,dayy, dadxxx,dadxxy,dadxyy,dadyyy
    ! From 1 to 10 of the third dimension of b denotes
    ! b,dbx,dby,dbxx,dbxy,dbyy,dbdxxx,dbdxxy,dbdxyy,dbdyyy
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION gpos(nx),x(nx,numnode),nv(numnode),ds(nx,numnode)
    DIMENSION a(mm,mm,10),b(mm,ndex,10)
    DIMENSION xv(nx,ndex),dif(nx,ndex),w(ndex,10),p(6,ndex),pp(mm,mm)

    DO i=1,ndex
        nn=nv(i)
        xv(1,i)=x(1,nn)
        xv(2,i)=x(2,nn)
        p(1,i)=1.0
        p(2,i)=xv(1,i)
        p(3,i)=xv(2,i)
        p(4,i)=xv(1,i)*xv(1,i)
        p(5,i)=xv(1,i)*xv(2,i)
        p(6,i)=xv(2,i)*xv(2,i)
        dif(1,i)=gpos(1)-xv(1,i)
        dif(2,i)=gpos(2)-xv(2,i)
    ENDDO

    CALL Weight_W1(dif,nv,ds,w,nx,ndex,numnode)
    ! ************* Compute b and its derivatives
    DO 20 ii=1,mm
        DO 20 jj=1,ndex
            DO 20 kk=1,10
                b(ii,jj,kk)=p(ii,jj)*w(jj,kk)
20          CONTINUE
            ! ************* Compute a and its derivatives

            a=0.

            DO 30 iii=1,ndex
                DO 40 ik=1,mm
                    DO 40 jk=1,mm
                        pp(ik,jk)=p(ik,iii)*p(jk,iii)
40                  CONTINUE
                    DO 50 ikk=1,mm
                        DO 50 jkk=1,mm
                            DO 50 kkk=1,10
                                a(ikk,jkk,kkk)=a(ikk,jkk,kkk)+w(iii,kkk)*pp(ikk,jkk)
50                          CONTINUE
30                      CONTINUE
                        RETURN
                    END
