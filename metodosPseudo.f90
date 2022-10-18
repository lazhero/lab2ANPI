
module pseudo
    use matrixUtilities
    implicit None
    contains 

    !
    ! Descripcion: Método utilizado para calcular la aproximacion de la pseudoinversa de una matriz mediante el método de homier
    ! Entradas:
    !           * A: La matriz a la cual calcular la pseudoinversa
    !           * iterMax: el numero maximo de iteraciones necesarias
    ! Salidas:
    !           * time: El tiempo requerido para ejecutar la aproximacion
    !           * iteraciones: las iteraciones utilizadas para el calculo
    !           * error: el error resultante
    !           * X: la matriz pseudoinversa
    function homier(A,iterMax,time,iterations,error) result(X)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        real*16,dimension(:,:),intent(in)::A
        real*16,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m) !obtiene la matrix identidad
        X0=initialValue(A) ! calculo del valor inicia
        X=homier_aux(A,X0,I,m,n,iterMax,time,iterations,error) !llama a la funcion auxiliar
      
      
    end function

    !
    ! Descripcion: Método auxiliar utilizado para calcular la aproximacion de la pseudoinversa de una matriz mediante el método de homier
    ! Entradas:
    !           * A: La matriz a la cual calcular la pseudoinversa
    !           * X: La aproximacion inicial
    !           * I: Una matrix identidad de dimensiones mxm
    !           * m: el numero de filas de la matriz A
    !           * n: el numero de columnas de la matriz A
    !           * iterMax: el numero maximo de iteraciones necesarias
    ! Salidas:
    !           * time: El tiempo requerido para ejecutar la aproximacion
    !           * iteraciones: las iteraciones utilizadas para el calculo
    !           * error: el error resultante
    !           * pseudo: la matriz pseudo inversa resultante
    function homier_aux(A,X,I,m,n,IterMax,time,iterations,error) result(pseudo)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real*16::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        double precision::normF
        real*16::init,finish
        
        call cpu_time(init)
        Xk=X

        
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=matmul(Xk,(I+0.5*matmul((I-Y),(I+matmul((2*I-Y),(2*I-Y))))))
            
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



    function toutonian(A,iterMax,time,iterations,error) result(X)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        real*16,dimension(:,:),intent(in)::A
        real*16,dimension(:,:),allocatable::X,X0,I
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
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real*16::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        double precision::normF
        real*16::init,finish
        
        call cpu_time(init)
        Xk=X
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=0.5*matmul(Xk,(9*I-matmul(Y,(16*I-matmul(Y,(14*I-matmul(Y,(6*I-Y))))))))
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            error=normF
           
            if( normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function


    function newton(A,iterMax,time,iterations,error) result(X)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        real*16,dimension(:,:),intent(in)::A
        real*16,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=newton_aux(A,X0,I,m,n,iterMax,time,iterations,error)
      
    end function


    function newton_aux(A,X,I,m,n,IterMax,time,iterations,error) result(pseudo)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real*16::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        double precision::normF
        real*16::init,finish
        
        
        call cpu_time(init)
        Xk=X
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=2*Xk-matmul(Xk,Y)
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            error=normF
            !print *,normF
            if( normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function


    function chebyshev(A,iterMax,time,iterations,error) result(X)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        real*16,dimension(:,:),intent(in)::A
        real*16,dimension(:,:),allocatable::X,X0,I
        integer::m,n,iterMax
        integer,dimension(2)::dims
        dims=shape(A)
        m=dims(1)
        n=dims(2)
        I=getIdentity(m)
        X0=initialValue(A)
        X=chebyshev_aux(A,X0,I,m,n,iterMax,time,iterations,error)
    end function


    function chebyshev_aux(A,X,I,m,n,IterMax,time,iterations,error) result(pseudo)
        real*16,intent(out)::time
        double precision,intent(out)::error
        integer,intent(out)::iterations
        integer::m,n,k,IterMax
        real*16::A(m,n),X(n,m),I(m,m),pseudo(n,m),Xk(n,m),Y(m,m)
        double precision::normF
        real*16::init,finish
        
        
        call cpu_time(init)
        Xk=X
        do k=1,IterMax
            iterations=k
            Y=matmul(A,Xk)
            Xk=matmul( Xk,(3*I -  matmul(Y,(3*I-Y ))))
            Y=matmul(A,Xk)
            normF=norm(matmul(Y,A)-A)
            error=normF
            !print *,normF
            if( normF<=0.00001) then
                Exit
            end if
        end do
        pseudo=Xk
        call cpu_time(finish)
        time=finish-init


    end function


end module pseudo


