
module pseudo
    use matrixUtilities
    implicit None
    contains 


    function homier(A,iterMax,time) result(X)
        real,intent(out)::time
        real,dimension(:,:),intent(in)::A
        real,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=homier_aux(A,X0,I,m,n,iterMax,time)
      
    end function


    function homier_aux(A,X,I,m,n,IterMax,time) result(pseudo)
        real,intent(out)::time
        integer::m,n,k,IterMax
        real::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        real::normF
        real::init,finish
        
        call cpu_time(init)
        Xk=X
        
        do k=1,IterMax
            Y=matmul(A,Xk)
            Xk=matmul(Xk,(I+0.5*matmul((I-Y),(I+(2*I-Y)**2))))
            
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            print *,normF
            if(k>1 .and. normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function



    function midpoint(A,iterMax,time) result(X)
        real,intent(out)::time
        real,dimension(:,:),intent(in)::A
        real,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=midpoint_aux(A,X0,I,m,n,iterMax,time)
      
    end function


    function midpoint_aux(A,X,I,m,n,IterMax,time) result(pseudo)
        real,intent(out)::time
        integer::m,n,k,IterMax
        real::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        real::normF
        real::init,finish
        
        call cpu_time(init)
        Xk=X
        do k=1,IterMax
            Y=matmul(A,Xk)
            Xk=2*Xk-matmul(Xk,Y)
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            print *,normF
            if(normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function


end module pseudo

program part1
    use matrixUtilities
    use pseudo
    implicit none

    real,dimension(:,:),allocatable::mat,x
    real::times
    integer::iterMax

    iterMax=250

    mat=generatedExampleMatrix()
    x =midpoint(mat,iterMax,times)
    call printMatrix(x)

    

   
        

       
    

end program part1
