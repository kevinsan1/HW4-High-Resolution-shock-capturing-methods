function f=limit(d1,d2,beta)
% 0 = first order, 1-2 = beta, 3= fromm, 4=vanleer, 5=vanalbada, 6=double
% minmod

if beta==0, %first order
    f=0;
elseif beta >= 1 & beta <= 2%beta: minmod (beta=1) and superbee (beta=2)
    if(d1*d2 < 0),
        f=0;
    else
        s=sign(d1);
        a=abs(d1);
        b=abs(d2);
        f=s*min(max([a b]),beta*min([a b]));
    end
elseif beta == 3 %Fromm
    f=0.5*(d1+d2);
elseif beta == 4 %vanleer
    if(d1*d2 <= 0),
        f=0;
    else
        f=2*d1*d2/(d1+d2);
    end
elseif beta ==  5 %vanalbada
     eps=1.e-20;
     f=(d1*(d2*d2+eps)+d2*(d1*d1+eps))/(d1*d1+d2*d2+2*eps);
elseif beta == 6 %double minmod
    if(d1*d2 < 0),
        f=0;
    else
        s=sign(d1);
        a=abs(d1);
        b=abs(d2);
        c=0.5*(a+b);
        f=s*min([2*a 2*b c]);
    end
end

