SUBROUTINE Compute_RadialBasis(x,y,xv,rk,ndex,R,q,nRBF,mbasis)
    !------------------------------------------------------------------
    ! Compute radial basis after added linear polynomial.
    ! Input: x,y - coordinates of the node considered
    !       xv[] - coordinates x and y for field nodes in the support domain
    !      ndex - number of field nodes in support domain
    !        r - RBF and its derivatives
    !        q - shape parameter
    !       nRBF- types of RBF (MQ, Exp, TSP)
    !       mbasis - number of monomials
    ! nRBF: 1: MQ; 2: Exp; 3: TPS
    !
    ! Output--rk[10,ndex+mbasis]
    ! From 1 to 10 denotes
    ! r,drx,dry,drdxx,drdxy,drdyy
    ! drdxxx,drdxxy, drdxyy, drdyyy
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION xv(2,ndex),rk(10,ndex+mbasis)

    rk=0

    DO 10 i=1,ndex
        rr2=(x-xv(1,i))**2+(y-xv(2,i))**2
        IF(nRBF.EQ.1) THEN ! MQ
            rk(1,i)=(rr2+R**2)**q
            rk(2,i)=2.*q*(rr2+R**2)**(q-1.)*(x-xv(1,i))
            rk(3,i)=2.*q*(rr2+R**2)**(q-1.)*(y-xv(2,i))
            rk(4,i)=2.*q*(rr2+R**2)**(q-1.)+4.*(q-1)*q* &
                (x-xv(1,i))**2*(rr2+R**2)**(q-2)
            rk(5,i)=4.*(q-1)*q*(x-xv(1,i))*(y-xv(2,i))* &
                (rr2+R**2)**(q-2)
            rk(6,i)=2.*q*(rr2+R**2)**(q-1.)+4.*q*(q-1)* &
                (y-xv(2,i))**2*(rr2+R**2)**(q-2)
        ENDIF
        IF(nRBF.EQ.2) THEN ! EXP
            rk(1,i)=exp(-q*rr2)
            rk(2,i)=-2.*q*exp(-q*rr2)*(x-xv(1,i))
            rk(3,i)=-2.*q*exp(-q*rr2)*(y-xv(2,i))
            rk(4,i)=-2*q*exp(-q*(rr2))+4*q*q*(x-xv(1,i))**2*exp(-q*rr2)
            rk(5,i)=4.*q*q*exp(-q*(rr2))*(y-xv(2,i))*(x-xv(1,i))
            rk(6,i)=-2*q*exp(-q*(rr2))+4*q*q*(y-xv(2,i))**2*exp(-q*rr2)
        ENDIF
        IF(nRBF.EQ.3) THEN ! TSP
            rk(1,i)=(rr2)**(0.5*q)
            rk(2,i)=q*(x-xv(1,i))*(rr2)**(0.5*q-1)
            rk(3,i)=q*(y-xv(2,i))*(rr2)**(0.5*q-1)
            rk(4,i)=q*(rr2)**(0.5*q-1)+2.*q*(0.5*q-1)*(x-xv(1,i))**2*(rr2) &
                **(0.5*q-2)
            rk(5,i)=q*(0.5*q-1)*(x-xv(1,i))*(y-xv(2,i))*(rr2)**(0.5*q-2)
            rk(6,i)=q*(rr2)**(0.5*q-1)+2.*q*(0.5*q-1)*(y-xv(2,i))**2*(rr2) &
                **(0.5*q-2)
        ENDIF
10  CONTINUE
    IF(mbasis.GT.0) THEN
        rk(1,ndex+1)=1.
        rk(1,ndex+2)=x
        rk(1,ndex+3)=y
        rk(2,ndex+2)=1.
        rk(3,ndex+3)=1.
    ENDIF
    RETURN
END
