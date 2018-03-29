function [Grrr,Grtt,Grpp,Grtr,Gtrr,Gtpp,Gttr,Gttt]=Gammafunc(r1,t1,H,f,f1,AA)
                % this function produces christoffle symbols
                Grrr=(2*r1*H^2-f1)/(2*f);
                Grpp=-r1*f+(r1)^3*(H)^2;
                Grtr=(r1*(H)*(-2*r1*(H)^2+f1))/(2*sqrt(f));
                Grtt=1/2*((r1^2*(H)^2)*(2*r1*(H)^2-f1)+f*(-2*r1*(H)^2+f1)-2*r1*sqrt(f)*AA);
                %Gppr=1/r;
                Gtrr=(H)/(f)^(3/2);
                Gtpp=((r1^2)*H)/sqrt(f);
                Gttr=(-2*r1*(H)^2+f1)/(2*f);
                Gttt=(r1*H*(2*r1*(H)^2-f1))/(2*sqrt(f));
end
