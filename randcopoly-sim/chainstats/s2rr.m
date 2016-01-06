function [valaa,valab,valba,valbb]=s2rr(N,NM,FA,LAM,k)

valeq=2*power(k,-2).*(-1+cos(k*NM)+NM*k.*sinint(k*NM));

valaa=FA*N*valeq;
valbb=(1-FA)*N*valeq;
valab=zeros(length(k),1);
valba=valab;


for J=2:N
    A=NM*(J-1)*k;
    
    valne=power(k,-2).*(-2*cos(A)+2*cos(k*NM).*cos(A)...
        -2*A.*sinint(A)+(A-k*NM).*sinint(A-k*NM)+(A+k*NM).*sinint(A+k*NM));

    valaa=valaa+2*(N-J+1)*(FA^2+FA*(1-FA)*LAM^(J-1))*valne;
    valab=valab+2*(N-J+1)*(FA*(1-FA)-FA*(1-FA)*LAM^(J-1))*valne;
    valba=valba+2*(N-J+1)*(FA*(1-FA)-FA*(1-FA)*LAM^(J-1))*valne;
    valbb=valbb+2*(N-J+1)*((1-FA)^2+FA*(1-FA)*LAM^(J-1))*valne;
end