SUBROUTINE GaussEqSolver_Sym(n,ma,a,b,ep,kwji)
    !------------------------------------------------------------------
    ! Solve sysmmetric linear equation ax=b by using Gauss elimination.
    ! If kwji=1, no solution;if kwji=0,has solution
    ! Input--n,ma,a(ma,n),b(n),ep,
    ! Output--b,kwji
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION a(ma,n),b(n),m(n+1)
    DO 10 i=1,n
10      m(i)=i
        DO 120 k=1,n
            p=0.0
            DO 20 i=k,n
                DO 20 j=k,n
                    IF(dabs(a(i,j)).GT.dabs(p)) THEN
                        p=a(i,j)
                        io=i
                        jo=j
                    ENDIF
20              CONTINUE
                IF(dabs(p)-ep) 30,30,35
30              kwji=1
                RETURN
35          CONTINUE
            IF(jo.EQ.k) GO TO 45
            DO 40 i=1,n
                t=a(i,jo)
                a(i,jo)=a(i,k)
                a(i,k)=t
40          CONTINUE
            j=m(k)
            m(k)=m(jo)
            m(jo)=j
45          IF(io.EQ.k) GO TO 55
            DO 50 j=k,n
                t=a(io,j)
                a(io,j)=a(k,j)
                a(k,j)=t
50          CONTINUE
            t=b(io)
            b(io)=b(k)
            b(k)=t
55          p=1./p
            in=n-1
            IF(k.EQ.n) GO TO 65
            DO 60 j=k,in
60              a(k,j+1)=a(k,j+1)*p
65              b(k)=b(k)*p
                IF(k.EQ.n) GO TO 120
                DO 80 i=k,in
                    DO 70 j=k,in
70                      a(i+1,j+1)=a(i+1,j+1)-a(i+1,k)*a(k,j+1)
80                      b(i+1)=b(i+1)-a(i+1,k)*b(k)
120                 CONTINUE
                    DO 130 i1=2,n
                        i=n+1-i1
                        DO 130 j=i,in
130                         b(i)=b(i)-a(i,j+1)*b(j+1)
                            DO 140 k=1,n
                                i=m(k)
140                             a(1,i)=b(k)
                                DO 150 k=1,n
150                                 b(k)=a(1,k)
                                    kwji=0
                                    RETURN
                                END
