
module pseudo
    use matrixUtilities
    implicit None
    contains 


    function homier(A,iterMax,time,iterations,error) result(X)
        real,intent(out)::time,error
        integer,intent(out)::iterations
        real,dimension(:,:),intent(in)::A
        real,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=homier_aux(A,X0,I,m,n,iterMax,time,iterations,error)
      
    end function


    function homier_aux(A,X,I,m,n,IterMax,time,iterations,error) result(pseudo)
        real,intent(out)::time,error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        real::normF
        real::init,finish
        
        call cpu_time(init)
        Xk=X
        
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=matmul(Xk,(I+0.5*matmul((I-Y),(I+(2*I-Y)**2))))
            
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            error=normF
            if(k>1 .and. normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function



    function toutonian(A,iterMax,time,iterations,error) result(X)
        real,intent(out)::time,error
        integer,intent(out)::iterations
        real,dimension(:,:),intent(in)::A
        real,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=toutonian_aux(A,X0,I,m,n,iterMax,time,iterations,error)
      
    end function


    function toutonian_aux(A,X,I,m,n,IterMax,time,iterations,error) result(pseudo)
        real,intent(out)::time,error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        real::normF
        real::init,finish
        
        call cpu_time(init)
        Xk=X
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=0.5*matmul(Xk,(9*I-matmul(Y,(16*I-matmul(Y,(14*I-matmul(Y,(6*I-Y))))))))
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            error=normF
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
    real::times,error
    integer::iterations
    integer::iterMax

    iterMax=250

    mat=generatedExampleMatrix()
    x =toutonian(mat,iterMax,times,iterations,error)
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    call printMatrix(x)

    print *,"_______________________________________________________________________________________"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    x=homier(mat,iterMax,times,iterations,error)
    call printMatrix(x)

    

   
        

       
    

end program part1
