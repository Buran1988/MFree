PROGRAM RPIM
    !------------------------------------------------------------------
    ! Main program for testing the RPIM shape function.
    ! Call Subroutine RPIM_ShapeFunc_2D( ).
    ! 25 field nodes (5X5) in domain [x,y]--[-1,1;-1,1].
    ! 61X61 sampling points are used to plot 2-D RPIM shape Func.
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    PARAMETER(nx=2,numnode=121)
    DIMENSION  x(nx,numnode),nv(numnode), gpos(nx),phi(10,numnode)

    OPEN(2,FILE='phi.dat') ! Output file
    WRITE(2,50)
    nRBF=1 ! Using MQ-RBF
    !nRBF=3
    !q=0.5
    q=0.98
    alfc=1.0
    dc=1.0
    mbasis=3 ! Number of basis
    xlength=10.
    ylength=10.
    ndivx=10
    ndivy=10
    xstep=xlength/ndivx
    ystep=ylength/ndivy
    nn=0

    DO i=1,ndivx+1
        DO j=1,ndivy+1
            nn=nn+1
            x(1,nn)= (i-1)*xstep !x coordinates of field nodes
            x(2,nn)= (j-1)*ystep !y coordinates of field nodes
        ENDDO
    ENDDO

    ! x (1:4) = (/ (sqrt (real (i)), i = 1, 4) /)
    !  x(1,nn)=(/(i, i=0,nn)/)

    DO i=1,numnode
        nv(i)=i ! Field nodes in support domain
    ENDDO

    ndex=25! number of field nodes in the support domain
    !nce=1
    !numnode/2+1 ! the node in the centre of 25 field nodes
    nm=10 ! interpolation points number
    ste=1. !9./(nm-1)

    !loop_through_points_x:
    DO ix=1,nm

        phi=0

        gpos(1)=0.4 + ste*(ix-1)

        !loop_through_points_y:

        DO j=1,nm

            gpos(2)=0.4 + ste*(j-1)
            !gpos(2)=0.4 + jxx-1

            IF((abs(gpos(1)).LE.10.0).AND.(abs(gpos(2)).LE.10.0)) THEN

                !  call DetectDomainNodes(x,nx,numnode)
                ds=1.0
                CALL SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)



                CALL RPIM_ShapeFunc_2D(gpos,x,nv,phi,nx,numnode,ndex,&
                    alfc,dc,q,nRBF, mbasis)

                write (*,*) 'Defined shape function for point: ', gpos(1), gpos(2)
                end if

            ! ***********Output RPIM shape function
            IF((abs(gpos(1)).LE.10.0).AND.(abs(gpos(2)).LE.10.0)) THEN
                DO kk=1,ndex
                    nd=nv(kk)
                    WRITE(2,100)nv(kk),x(1,nd),x(2,nd),phi(1,kk), &
                        phi(2,kk),phi(3,kk),phi(4,kk),phi(6,kk),gpos(1),gpos(2)
                ENDDO
            ENDIF

        ENDDO !loop_through_points_y
    ENDDO !loop_through_points_x

    WRITE(2,150)
50  FORMAT(1x,'Node', 5x,'x', 7x,'y', 8x,'Phi', 6x,'dPhidx', &
        5x,'dPhidy', 4x, 'dPhidxx', 4x,'dPhidyy',5x,'GposX',5x,'GposY',    /,80('-'))
100 FORMAT(1x,i4, 2f8.3, 7f11.5)
150 FORMAT(80('-'))
    WRITE(*,*)'Ready!'
END PROGRAM RPIM
