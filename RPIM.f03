PROGRAM RPIM
    !------------------------------------------------------------------
    ! Main program for testing the RPIM shape function.
    ! Call Subroutine RPIM_ShapeFunc_2D( ).
    ! 25 field nodes (5X5) in domain [x,y]--[-1,1;-1,1].
    ! 61X61 sampling points are used to plot 2-D RPIM shape Func.
    !------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    parameter(nx=2,numnode=25)
    dimension x(nx,numnode),nv(numnode), gpos(nx),phi(10,numnode)
    open(2,file='phi.dat') ! Output file
    write(2,50)
    nRBF=1 ! Using MQ-RBF

    q=2.03
    alfc=0.03
    dc=0.5
    mbasis=0 ! Number of basis
    xlength=2.
    ylength=2.
    ndivx=4
    ndivy=4
    xstep=xlength/ndivx
    ystep=ylength/ndivy
    nn=0
    do i=1,ndivx+1
        do j=1,ndivy+1
            nn=nn+1
            x(1,nn)=-1.+(i-1)*xstep !x coordinates of field nodes
            x(2,nn)=-1.+(j-1)*ystep !y coordinates of field nodes
        enddo
    enddo
    do i=1,numnode
        nv(i)=i ! Field nodes in support domain
    enddo
    ndex=25
    nce=numnode/2+1 ! the node in the centre of 25 field nodes
    nm=61
    ste=2./(nm-1)
    do ix=1,nm


        phi=0

        gpos(1)=-1.+ste*(ix-1)
        do j=1,nm
            gpos(2)=-1.+ste*(j-1)
            if((abs(gpos(1)).le.1).and.(abs(gpos(2)).le.1)) then
                call RPIM_ShapeFunc_2D(gpos,x,nv,phi,nx,numnode,ndex,&
                    alfc,dc,q,nRBF, mbasis)
            else
            endif
            ! ***********Output RPIM shape function
            if((abs(gpos(1)).le.1.0).and.(abs(gpos(2)).le.1.0)) then
                do kk=1,ndex
                    nd=nv(kk)
                    write(2,100)nv(kk),x(1,nd),x(2,nd),phi(1,kk), &
                        phi(2,kk),phi(3,kk),phi(4,kk),phi(6,kk),gpos(1),gpos(2)
                enddo
            endif
        enddo
    enddo
    write(2,150)
50  format(1x,'Node', 5x,'x', 7x,'y', 8x,'Phi', 6x,'dPhidx', &
        5x,'dPhidy', 4x, 'dPhidxx', 4x,'dPhidyy',5x,'GposX',5x,'GposY',    /,80('-'))
100 format(1x,i4, 2f8.3, 7f11.5)
150 format(80('-'))
    write(*,*)'Ready!'
END PROGRAM RPIM
