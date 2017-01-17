SUBROUTINE Weight_W2(dif,nv,ds,w,nx,ndex,numnode)
    !------------------------------------------------------------------
    ! Quartic spline weight function
    ! input--dif,nv,ds,nx,ndex,numnode
    ! output--w
    ! from 1 to 10 column of w denotes w,dwdx,dwdy,dwdxx,dwdxy,dwdyy
    !------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    dimension dif(nx,ndex),nv(numnode),ds(nx,numnode),w(ndex,10)
    ep=1.0e-20
    do 10 i=1,ndex
        nn=nv(i)
        difx=dif(1,i)
        dify=dif(2,i)
        if(dabs(difx).le.ep) then
            drdx=0.
        else
            drdx=(difx/dabs(difx))/ds(1,nn)
        end if
        if (dabs(dify).le.ep) then
            drdy=0.
        else
            drdy=(dify/dabs(dify))/ds(2,nn)
        end if
        rx=dabs(dif(1,i))/ds(1,nn)
        ry=dabs(dif(2,i))/ds(2,nn)
        wx=1.-6.*rx*rx+8.*rx*rx*rx-3.*rx*rx*rx*rx
        dwxdx=(-12.*rx+24.*rx*rx-12.*rx*rx*rx)*drdx
        dwxdxx=(-12.+48.*rx-36.*rx*rx)/(ds(1,nn)*ds(1,nn))
        dwxdxxx=(48.-72*rx)*drdx**3
        wy=1.-6.*ry*ry+8.*ry*ry*ry-3.*ry*ry*ry*ry
        dwydy=(-12.*ry+24.*ry*ry-12.*ry*ry*ry)*drdy
        dwydyy=(-12.+48.*ry-36.*ry*ry)/(ds(2,nn)*ds(2,nn))
        dwydyyy=(48.-72*ry)*drdy**3
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
10  continue
    return
end
