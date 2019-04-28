function [V,H]=arnoldiMR(A,B,r)
% This is a program for model redcution of a first order system using
% Arnoldi factorization
% Given the first order system (A, B, C, D)
% This program outputs V such that
% Ar=V'AV  Br=V'B  Cr= CV  Dr=D


v=B/norm(B);

w=A*v; alpha=v'*w; f=w-v*alpha; H=alpha; V=v;

for i=2:r
    beta=norm(f); v=f/beta;
    
    q2=v-V*V'*v;
    q2=q2/norm(q2);
    
    ssr=V'*q2;
    q2=q2-V*ssr;
    
    V=[V q2];
    ek=eye(i-1);
    H=[H;beta*ek(i-1,:)];
    w=A*v;
    hhat=V'*w;
    H=[H hhat];
    f=w-V*hhat;
    i
end
