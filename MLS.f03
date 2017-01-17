PROGRAM MLS

    !------------------------------------------------------------------
    ! Main program for testing the MLS shape function.
    ! Call Subroutine MLS_ShapeFunc_2D( ).
    ! 25 field nodes (5X5) in domain [x,y]--[-1,1;-1,1].
    ! 61X61 interpolation points are used to plot 2-D MLS shape Func.
    !------------------------------------------------------------------


    IMPLICIT REAL*8 (A-H,O-Z)
    PARAMETER(nx=2,numnode=25)
    DIMENSION x(nx,numnode),nv(numnode), gpos(nx)
    DIMENSION phi(10,numnode),ds(nx,numnode)
    OPEN(2,FILE='phi.dat') ! Output file
    WRITE(2,50)
    mm=3 ! Number of basis
    xlength=2.
    ylength=2.
    ndivx=4
    ndivy=4
    xstep=xlength/ndivx
    ystep=ylength/ndivy
    nn=0
    DO i=1,ndivx+1
        DO j=1,ndivy+1
            nn=nn+1
            x(1,nn)=-1.+(i-1)*xstep !x coordinates of field nodes
            x(2,nn)=-1.+(j-1)*ystep !y coordinates of field nodes
        ENDDO
    ENDDO
    DO i=1,numnode
        nv(i)=i
        ds(1,i)=0.
        ds(2,i)=0.
    ENDDO
    ndex=25
    DO j=1,numnode
        xn=x(1,j)
        yn=x(2,j)
        rx0=abs(xn-1)
        IF(rx0.LT.abs(xn+1)) rx0=abs(xn+1)
        ry0=abs(yn-1)
        IF(ry0.LT.abs(yn+1)) ry0=abs(yn+1)
        ds(1,j)=rx0 ! rw for weight function (support domain)
        ds(2,j)=ry0
    ENDDO
    nce=numnode/2+1 ! the node in the centre of 25 field nodes
    nm=61
    ste=2./(nm-1)
    DO ix=1,nm

        phi=0
        gpos(1)=-1.+ste*(ix-1)
        DO j=1,nm
            gpos(2)=-1.+ste*(j-1)
            IF((abs(gpos(1)).LE.1).AND.(abs(gpos(2)).LE.1)) THEN
                CALL MLS_ShapeFunc_2D(gpos,x,nv,ds,phi,nx,numnode,ndex,mm)
            ELSE
            ENDIF
            ! ***********Output MLS shape function
            IF((abs(gpos(1)).LE.1.0).AND.(abs(gpos(2)).LE.1.0)) THEN
                DO kk=1,ndex
                    nd=nv(kk)
                    WRITE(2,100)nv(kk),x(1,nd),x(2,nd),phi(1,kk), &
                        phi(2,kk),phi(3,kk),phi(4,kk),phi(6,kk),gpos(1),gpos(2)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    WRITE(2,150)
50  FORMAT(1x,'Node', 5x,'x', 7x,'y', 8x,'Phi', 6x,'dPhidx', &
        5x,'dPhidy', 4x, 'dPhidxx', 4x,'dPhidyy',5x,'GposX',5x,'GposY',    /,80('-'))
100 FORMAT(1x,i4, 2f8.3, 7f11.5)
150 FORMAT(80('-'))


    WRITE(*,*)'Ready!'
END PROGRAM MLS
