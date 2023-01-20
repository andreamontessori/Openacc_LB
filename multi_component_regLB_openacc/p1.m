clear all
psi=load('psi.out');
psi=reshape(psi,101,101);
imagesc(psi)
psi_x=zeros(101,101);
psi_y=zeros(101,101);
location=zeros(100,100);
for i=2:100
    for j=2:100
    psi_x(i,j)=(3.0)*((1/9)*(psi(i+1,j)-psi(i-1,j)) + (1/36)*(psi(i+1,j+1) + psi(i+1,j-1)-psi(i-1,j+1)-psi(i-1,j-1)));
    psi_y(i,j)=(3.0)*((1/9)*(psi(i,j+1)-psi(i,j-1)) + (1/36)*(psi(i+1,j+1) - psi(i+1,j-1)+psi(i-1,j+1)-psi(i-1,j-1)));
    end
end

absPsi=sqrt(psi_x(:,:).^2+psi_y(:,:).^2);
for i=2:100
    for j=2:100
        if(absPsi(i,j)>0.001 && absPsi(i,j)>absPsi(i,j) && absPsi(i-1,j)>absPsi(i,j))
            location(i,j)=1;
            
        end
    end
end
