Subroutine GaussEqSolver_Sym(n,ma,a,b,ep,kwji)
    !------------------------------------------------------------------
    ! Solve sysmmetric linear equation ax=b by using Gauss elimination.
    ! If kwji=1, no solution;if kwji=0,has solution
    ! Input--n,ma,a(ma,n),b(n),ep,
    ! Output--b,kwji
    !------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    dimension a(ma,n),b(n),m(n+1)
    do 10 i=1,n
10      m(i)=i
        do 120 k=1,n
            p=0.0
            do 20 i=k,n
                do 20 j=k,n
                    if(dabs(a(i,j)).gt.dabs(p)) then
                        p=a(i,j)
                        io=i
                        jo=j
                    endif
20              continue
                if(dabs(p)-ep) 30,30,35
30              kwji=1
                return
35          continue
            if(jo.eq.k) go to 45
            do 40 i=1,n
                t=a(i,jo)
                a(i,jo)=a(i,k)
                a(i,k)=t
40          continue
            j=m(k)
            m(k)=m(jo)
            m(jo)=j
45          if(io.eq.k) go to 55
            do 50 j=k,n
                t=a(io,j)
                a(io,j)=a(k,j)
                a(k,j)=t
50          continue
            t=b(io)
            b(io)=b(k)
            b(k)=t
55          p=1./p
            in=n-1
            if(k.eq.n) go to 65
            do 60 j=k,in
60              a(k,j+1)=a(k,j+1)*p
65              b(k)=b(k)*p
                if(k.eq.n) go to 120
                do 80 i=k,in
                    do 70 j=k,in
70                      a(i+1,j+1)=a(i+1,j+1)-a(i+1,k)*a(k,j+1)
80                      b(i+1)=b(i+1)-a(i+1,k)*b(k)
120                 continue
                    do 130 i1=2,n
                        i=n+1-i1
                        do 130 j=i,in
130                         b(i)=b(i)-a(i,j+1)*b(j+1)
                            do 140 k=1,n
                                i=m(k)
140                             a(1,i)=b(k)
                                do 150 k=1,n
150                                 b(k)=a(1,k)
                                    kwji=0
                                    return
                                END
