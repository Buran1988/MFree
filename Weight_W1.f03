SUBROUTINE Weight_W1(dif,nv,ds,w,nx,ndex,numnode)
    !------------------------------------------------------------------
    ! Cubic spline weight function
    ! input--dif,nv,ds,nx,ndex,numnode
    ! output--w
    ! from 1 to 10 column of w denotes w,dwdx,dwdy,dwdxx,dwdxy,dwdyy
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION dif(nx,ndex),nv(numnode),ds(nx,numnode),w(ndex,10)
    ep=1.0E-20

    DO 10 i=1,ndex
        nn=nv(i)
        difx=dif(1,i)
        dify=dif(2,i)
        IF(dabs(difx).LE.ep) THEN
            drdx=0.
        ELSE
            drdx=(difx/dabs(difx))/ds(1,nn)
        END IF
        IF (dabs(dify).LE.ep) THEN
            drdy=0.
        ELSE
            drdy=(dify/dabs(dify))/ds(2,nn)
        END IF
        rx=dabs(dif(1,i))/ds(1,nn)
        ry=dabs(dif(2,i))/ds(2,nn)
        IF(rx.GT.0.5) THEN
            wx=(4./3.)-4.*rx+4.*rx*rx-(4./3.)*rx*rx*rx
            dwxdx=(-4.+8.*rx-4.*rx*rx)*drdx
            dwxdxx=(8.-8.*rx)*drdx*drdx
            dwxdxxx=(-8.)*drdx*drdx*drdx
        ELSE IF(rx.LE.0.5) THEN
            wx=(2./3.)-4.*rx*rx+4.*rx*rx*rx
            dwxdx=(-8.*rx+12.*rx*rx)*drdx
            dwxdxx=(-8.+24.*rx)*drdx*drdx
            dwxdxxx=(24.)*drdx*drdx*drdx
        END IF
        IF(ry.GT.0.5) THEN
            wy=(4./3.)-4.*ry+4.*ry*ry-(4./3.)*ry*ry*ry
            dwydy=(-4.+8.*ry-4.*ry*ry)*drdy
            dwydyy=(8.-8.*ry)*drdy*drdy
            dwydyyy=(-8.)*drdy*drdy*drdy
        ELSE IF(ry.LE.0.5) THEN
            wy=(2./3.)-4.*ry*ry+4.*ry*ry*ry
            dwydy=(-8.*ry+12.*ry*ry)*drdy
            dwydyy=(-8.+24.*ry)*drdy*drdy
            dwydyyy=(24.)*drdy*drdy*drdy
        END IF
        w(i,1)=wx*wy
        w(i,2)=wy*dwxdx
        w(i,3)=wx*dwydy
        w(i,4)=wy*dwxdxx
        w(i,5)=dwxdx*dwydy
        w(i,6)=wx*dwydyy
        w(i,7)=wy*dwxdxxx
        w(i,8)=dwxdxx*dwydy
        w(i,9)=dwxdx*dwydyy
        w(i,10)=wx*dwydyyy
10  CONTINUE
    RETURN
END
