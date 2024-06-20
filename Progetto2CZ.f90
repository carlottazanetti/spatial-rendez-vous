module dati 
    implicit none 
    real*8,parameter:: h1=200.d0*10.d0**3 ,h2=600.d0*10.d0**3 ,Rearth=6378.39*10.d0**3 !m
    real*8,parameter:: r1=(200.d0+6378.39)*10.d0**3, r2=(600.d0+6378.39)*10.d0**3, Mearth=5.97e24, G=6.67e-11 !Kg e N*m^2*kg^-2
    real*8,parameter::Trv=acos(-1.d0)*sqrt(((200.d0+6378.39+600.d0+6378.39)*10.d0**(3))**3/(8.d0*6.67e-11*5.97e24)) !in secondi
    real*8,parameter::v1=sqrt(6.67e-11*5.97e24/((200.d0+6378.39)*10.d0**3)), v2=sqrt(6.67e-11*5.97e24/((600.d0+6378.39)*10.d0**3)) !in m/s
    real*8,parameter::v0=r2/trv
    save 
end module dati 

module RHS 
    use dati
    implicit none 

    contains 
    subroutine dydx(neq,x,y,f)
        integer, intent(in):: neq 
        real*8, intent(in):: x
        real*8,dimension(neq), intent(in):: y  
        real*8,dimension(neq), intent(out):: f 

        !y=(x,y,vx,vy)
        !f=(x',y',vx',vy')
        !scrivo il sistema adimensionalizzato
        f(1)=y(3)
        f(2)=y(4)
        f(3)=-G*Mearth*trv*y(1)/(sqrt((y(1)**2+y(2)**2))**3*v0*r2**2)
        f(4)=-G*Mearth*trv*y(2)/(sqrt((y(1)**2+y(2)**2))**3*v0*r2**2)

    end subroutine dydx 
end module RHS

module ODE_solver 
    use dati
    use RHS 
    implicit none 

    contains 
    subroutine rk4(neq,h,x,yold,ynew)
        integer, intent(in):: neq 
        real*8,intent(in):: h  , x
        real*8,dimension(neq),intent(in):: yold
        real*8,dimension(neq),intent(out):: ynew 
        real*8,dimension(neq):: k1,k2,k3,k4
        integer:: i 

        call dydx(neq,x,yold,k1) !k1=f(xn,yn) = f(yold)
        do i=1, neq 
            ynew(i)=yold(i)+h*0.5d0*k1(i) 
        end do
        call dydx(neq,x+h*0.5d0,ynew,k2) !k2=f(yn+0.5*h*k1, xn+0.5*h)= f(ynew)
        do i=1, neq 
            ynew(i)=yold(i)+h*0.5d0*k2(i) 
        end do
        call dydx(neq,x+h*0.5d0,ynew,k3) !k3= f(yn+0.5*h*k2, xn+0.5*h) = f(ynew)
        do i=1, neq 
            ynew(i)=yold(i)+h*k3(i)
        end do
        call dydx(neq,x+h,ynew,k4) !k4= f(yn+h*k3, xn+h) = f(ynew)
        do i=1, neq 
            ynew(i)=yold(i)+h*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))/6.d0 !y(n+1)=yn+h/6(k1+2k2+2k3+k4)
        end do
    end subroutine rk4
    
    subroutine save_results(filename,npoints,y,x)!subroutine che salva i risultati con 2 colonne 
        integer, intent(in):: npoints 
        character (len=*), intent(in):: filename 
        real*8,dimension(0:npoints),intent(in):: x, y 
        integer:: i 
 
        open(40,file=filename)
        do i=0,npoints 
            write(40,'(2(1pe20.10))') x(i), y(i)
        end do
        close(40)
    end subroutine save_results

    subroutine save_results2(filename,npoints,y,x,err)!subroutine che salva i risultati con 3 colonne
        integer, intent(in):: npoints 
        character (len=*), intent(in):: filename 
        real*8,dimension(0:npoints),intent(in):: x, y , err
        integer:: i 
    
        open(40,file=filename)
        do i=0,npoints 
            write(40,'(3(1pe20.10))') x(i), y(i), err(i)
        end do
        close(40)
    end subroutine save_results2
end module ODE_solver

subroutine shooting (neq,npoints,t,y,delta_v)
    use dati
    use ODE_solver 
    implicit none
    integer:: i , neq, npoints 
    real*8, dimension(0:npoints):: t
    real*8, dimension(neq,0:npoints):: y 
    real*8:: h , tmin, tmax , dist, toll 
    real*8:: w1, w2, y1, y2,y3, yn ,xn , w3 , delta_v

    !voglio integrare da 0 a 1 
    tmin=0.0d0
    tmax=1.d0
    h=(tmax-tmin)/npoints !passo
    do i=0,npoints 
        t(i)= tmin+i*h 
    end do
    !y=(x,y,vx,vy)
    !f=(x',y',vx',vy')
    !condizioni al contorno adimensionalizzate
    xn=-1.d0 !posizione della navicella target al tempo Trv sull'asse x
    yn=0.d0 !posizione della navicella target al tempo Trv sull'asse y
    y(1,0)=r1/r2 !x(0)
    y(2,0)=0.0d0  !y(0)
    y(3,0)=0.d0 !vx(0)
    y(4,0)=v1/v0 !condizione iniziale arbtitraria di vy(0)
    w1=y(4,0) !mi salvo il dato arbitrario
    do i=1,npoints 
        call rk4(neq,h,t(i),y(:,i-1),y(:,i)) !mettendo ':' ciclo tutti e 4 i vettori 
    end do !questo metodo prende in input yn e dà in outpt y(n+1)
    y1=y(2,npoints) !mi salvo il risultato (y(trv))

    y(1,0)=r1/r2 !x(0)
    y(2,0)=0.0d0  !y(0)
    y(3,0)=0.d0 !vx(0)
    y(4,0)=v1/v0+0.7 !condizione iniziale arbtitraria di vy(0)
    w2=y(4,0)
    do i=1,npoints 
        call rk4(neq,h,t(i),y(:,i-1),y(:,i))
    end do 
    y2=y(2,npoints) !salvo il risultato

    w3=w2+(w2-w1)/(y2-y1)*(yn-y2) 
    y(1,0)=r1/r2 !x(0)
    y(2,0)=0.0d0  !y(0)
    y(3,0)=0.d0 !vx(0)
    y(4,0)=w3
    do i=1,npoints 
        call rk4(neq,h,t(i),y(:,i-1),y(:,i))
    end do
    y3=y(2,npoints)
   
    dist=3.d10 !Valore alto per entrare nel ciclo 
    toll=3.d0/r2 !Tolleranza di 3m adimensionalizzata
    do while (toll<dist)
        w1=w2
        w2=w3
        y1=y2
        y2=y3
        w3=w2+(w2-w1)/(y2-y1)*(yn-y2)
        y(4,0)=w3
        do i=1,npoints 
            call rk4(neq,h,t(i),y(:,i-1),y(:,i))
        end do
        y3=y(2,npoints) 
        dist=abs(sqrt((y3-yn)**2+(y(1,npoints)-xn)**2)) !differenza di posizione tra le due navicelle adimensionalizzata
    end do
    delta_v=y(4,0)-v1/v0
    delta_v=delta_v*v0 !riporto nella giusta unità di misura
    y(1,:)=y(1,:)*r2
    y(2,:)=y(2,:)*r2
    y(3,:)=y(3,:)*v0
    y(4,:)=y(4,:)*v0
    t=t*trv
end subroutine shooting 

subroutine angle(y,neq,npoints,angolo) !subroutine per trovare l'angolo tra due vettori 
    implicit none 
    integer:: neq,npoints,i,j 
    integer,parameter::nd=2
    real*8, dimension(neq,0:npoints):: y 
    real*8, dimension(0:npoints):: angolo
    real*8, dimension(nd):: a,b 
    real*8:: norma_a, norma_b , num, cos_angolo
    
    do j=0,npoints
        norma_a= 0.d0
        norma_b=0.d0
        a(1)=y(1,j)!il vettore a ha come componenti x e y 
        a(2)=y(2,j)
        b(1)=y(3,j) !Il vettore b ha come componenti vx e vy
        b(2)=y(4,j)
        do i=1,nd
            norma_a=norma_a+a(i)**2
            norma_b=norma_b+b(i)**2
        end do
        norma_a=sqrt(norma_a)
        norma_b=sqrt(norma_b)
        num=0.d0
        do i=1,nd
            num=num+a(i)*b(i)
        end do
        cos_angolo=num/(norma_a*norma_b)
        angolo(j)=acos(cos_angolo)
    end do
end subroutine angle

subroutine distance(y,npoints,neq,delta_v2,dist)!subroutine per calcolare la distanza finale tra i due veicoli con delta v dversi
    use dati
    use ODE_solver
    implicit none 
    real*8, dimension(neq,0:npoints):: y 
    integer:: k , neq, npoints 
    real*8, dimension(0:npoints):: t
    real*8:: xn, yn, delta_v2, tmin, tmax, x , dist, h

    tmin=0.0d0
    tmax=1.d0
    h=(tmax-tmin)/npoints !passo
    do k=0,npoints 
        t(k)= tmin+k*h 
    end do
    !condizioni al contorno adimensionalizzate
    xn=-1.d0 !posizione della navicella target al tempo Trv sull'asse x
    yn=0.d0 !posizione della navicella target al tempo Trv sull'asse y
    y(1,0)=r1/r2 !x(0)
    y(2,0)=0.0d0  !y(0)
    y(3,0)=0.d0 !vx(0)
    y(4,0)=v1/v0+delta_v2/v0
    do k=1,npoints 
        call rk4(neq,h,t(k),y(:,k-1),y(:,k))
    end do 
    dist=abs(sqrt((y(2,npoints)-yn)**2+(y(1,npoints)-xn)**2)) !differenza di posizione tra le due navicelle adimensionalizzata
    dist=dist*r2 !metto la distanza in metri
end subroutine distance

program esame
    use dati
    use ODE_solver
    implicit none
    integer,parameter:: neq=4
    integer:: i , npoints, j , k , nd !npoints non è parameter perché verrà cambiato 
    real*8, allocatable:: t(:), energia(:), angolo(:), L(:) , errore_energia(:), errore_momento(:), tolleranza(:)
    real*8, allocatable:: y (:,:)
    character(len=14):: nome(0:9)
    real*8:: delta_v, delta_v2, dist, diff_delta_v(1000), distanze(1000), passo, toll_min, toll_max

    nome(0)='traiett_h1.txt' !nomi che mi servitanno per salvare i file nel ciclo 
    nome(1)='traiett_h2.txt'
    nome(2)='dist_h1.txt'
    nome(3)='dist_h2.txt'
    nome(4)='energia_h1.txt'
    nome(5)='energia_h2.txt'
    nome(6)='momento_h1.txt'
    nome(7)='momento_h2.txt'
    nome(8)='delta_v_h1.txt'
    nome(9)='delta_v_h2.txt'

    do i=0,1 !ciclo in cui il passo è inizialmente h1(definito nella subrouine dello shooting) e poi h2=0.1*h1
        npoints=50*10.d0**i !npoints è 50 al primo ciclo, e 500 al secondo 
        print*, 'I calcoli saranno effettuati con un passo di ', 1.d0/npoints
        if (allocated(t)) deallocate(t)
        if (allocated(energia)) deallocate(energia)
        if (allocated(angolo)) deallocate(angolo)
        if (allocated(L)) deallocate(L)
        if (allocated(errore_energia)) deallocate(errore_energia)
        if (allocated(errore_momento)) deallocate(errore_momento)
        if (allocated(y)) deallocate(y)
        if (allocated(tolleranza)) deallocate(tolleranza)
        allocate(t(0:npoints),energia(0:npoints),angolo(0:npoints),L(0:npoints))
        allocate(errore_energia(0:npoints),errore_momento(0:npoints),y(neq,0:npoints))

        !Pima e seconda parte: calcolo del valore di delta v 
        call shooting (neq,npoints,t,y,delta_v)
        print*,'Per ottenere un rendez-vous di successo,delta v risulta essere di ', delta_v, 'm/s'
        call save_results(nome(i),npoints,y(2,:)*0.001,y(1,:)*0.001) !salvo la traiettoria del veicolo 1 in km
        call save_results(nome(i+2),npoints, ((sqrt(y(1,:)**2+y(2,:)**2))-Rearth)*0.001 ,t/60.d0) !salvo la distanza del veicolo 1 in funzione del tempo in minuti 
        print*, 'Ho salvato la traiettoria del veicolo 1 e la sua distanza dalla superficie terrestre in km nei file ',nome(i),&
        &' e ',nome(i+2)
        print*, ' '

        !Terza parte: calcolo energia e quantità di moto
        energia=(0.5d0*(y(4,:)**2+y(3,:)**2))-(G*Mearth/sqrt(y(1,:)**2+y(2,:)**2)) 
        errore_energia=abs((energia(:)-energia(0))/energia(0))
        call save_results2(nome(i+4),npoints,energia,t/60.d0,errore_energia) !salvo l'energia con il tempo in minuti e l'errore
        print*, "Ho salvato i valori dell'energia in joule in funzione del tempo in minuti con relativo errore nel file ",nome(i+4)
        call angle(y,neq,npoints,angolo)
        L=sqrt(y(1,:)**2+y(2,:)**2)*sqrt(y(3,:)**2+y(4,:)**2)*sin(angolo)
        errore_momento=abs((L(:)-L(0))/L(0))
        call save_results2(nome(i+6),npoints,L,t/60.d0,errore_momento) ! Salvo il momento in funzione del tempo in minuti
        print*,'Ho salvato i valori della quantità di moto in kg*m/s in funzione del tempo in minuti con &
        & relativo errore nel file ',nome(i+6)
        print*, ' '

        !Quinta parte
        passo=(delta_v+0.01-delta_v+0.01)/1000.d0
        nd=0 !utilizzo questo indice per calcolare quanti valori di delta v mi fanno ottenere un rendez-vous di successo
        do j=1,1000
            delta_v2=delta_v-0.01+passo*j
            call distance(y,npoints,neq,delta_v2,dist)
            distanze(j)=dist !In metri
            diff_delta_v(j)=(delta_v2-delta_v)*100.d0 !in cm/s
            if (dist<3.d0) then
                nd=nd+1                
            end if
        end do
        allocate(tolleranza(nd))
        k=0 
        do j=1,1000
            if (distanze(j)<3.d0) then
                k=k+1
                tolleranza(k)=(diff_delta_v(j)/100.d0)+delta_v
            end if
        end do
        print*, 'Per ottenere un rendez-vous che abbia successo, delta v deve essere compreso tra', minval(tolleranza),' e ',&
        & maxval(tolleranza), 'm/s.'
        call save_results(nome(8+i),1000,distanze,diff_delta_v)!salvo le distanze in metri finale in funzione della differenza tra i delta v in cm/s
        print*,'Ho salvato i valori delle distanze finali tra i veicoli in metri in funzione &
        & della differenza tra i delta v in cm/s nel file ',nome(8+i)
        print*, ' '
        print*, ' '
        print*, ' '
    end do
    print*, 'the end'
end program esame
