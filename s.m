clc, clear all
C1=[-30, 152, -254, 141];
a=1.5;
b=2;
A=Budan_Fourier(C1,a,b);
if A==0
    fprintf('no real root in the interval')
else
    fprintf('maximum real roots in (%2.2f, %2.2f] is %d',a,b, A)
end
function [m] = sign_chng(x)
c=0;
    for j=1:length(x)-1
        if ((x(j)>0)&&(x(j+1)<0))||((x(j)<0)&&(x(j+1)>0))
          c=c+1;
        else
            continue
        end
    end
    m=c;
end

function w = Budan_Fourier(Y,m,n)
    syms x
    P= poly2sym(Y, x);
    DP=diff(P);
    VL(1)=subs(P,x,m); 
    VL(2)=subs(DP,x,m);
    VR(1)=subs(P,x,n); 
    VR(2)=subs(DP,x,n);
 for i=3:length(Y)
    [R1,Q1] = polynomialReduce(P,DP);
    VL(i)=subs(-R1,x,m); 
    VR(i)=subs(-R1,x,n); 
    d=coeffs(R1);
    if length(d)~=1
        P=DP; DP=-R1;
    else
        break
    end
 end
    fprintf('VS(a):[%s]\n', join(string(VL),','));
    fprintf('VS(b):[%s]\n',join(string(VR),','))
    VS1=nonzeros(VL); VS2=nonzeros(VR);
    fprintf('n.s.c at %2.2f is %d\n, %2.2f is %d\n',m,sign_chng(VS1),n,sign_chng(VS2));
    w=sign_chng(VS1);S2=sign_chng(VS2);
end