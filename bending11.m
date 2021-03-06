%%%%%%%%%%%%%%%%%%%%%%%%
%initial conditions
format longe
tester1=0;
RSL=1000;%initial source to black hole distance
eps=.001;%accuracy checking parameter
b=1;
M=(4.779)*10^(-6);%Mass of the black hole
A(1)=-2*10^-8;%Acceleration of expanding universe
L=1.1; %angular momentum
hub=(7/6)*10^-4;%stepsize in Hubble parameter
angle=sqrt(4*M/2000);%initial shooting angle
m=1; q=1; x=1;
H0(1)=0; %Hubble parameter
tester(1)=0;
dA=.5*10^(-8);%acceleration stepsize
delta=10^-6;
ss=1;
s=1;
dr=0;
%%%%%%%%%%%%%%%%
%runge-kutta solution and equal distance condition
while ss<=9  
    initspeed=-abs(L/RSL)/tan(angle);
    if ss>=2
        initspeed=speedo;
    end
    while s<=3
        tester1=0;
        dr=0;
        %speedo=0;
         while (abs(tester1-RSL) > eps) %& b<20)
            if (abs(tester1-RSL)>=.002)
              ur1=initspeed+tanh(dr/50)/20;
            else
             ur1=initspeed+dr*abs(dr)/10;
            end
            ur=ur1;
            speedo=ur;
            initspeed=ur1;
            x=x+1;
            h=.00001;
            phi=0;
            r=1000;
            t=0;
            f=1-2*M/r;
            f1=2*M/(r)^2;
            uphi=L/(r^2);
            H=H0(s)+A(ss)*t;
            ur=ur1;
            ut=[-(r*H*ur/sqrt(f))+sqrt([r^2*H^2*ur^2/f]+([ur^2/f]+[L^2/r^2])*(f-r^2*H^2))]/(f-r^2*H^2);  
            while (phi<=pi)%condition for crossing horizontal line 
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Gamma Functions (determining the geometry)
                AA=A(ss);
                r1=r;
                t1=t;
                f=1-2*M/r;
                H=H0(s)+A(ss)*t;
                f1=2*M/(r)^2;
                uphi=L/(r^2);
                [Grrr,Grtt,Grpp,Grtr,Gtrr,Gtpp,Gttr,Gttt]=Gammafunc(r1,t1,H,f,f1,AA);
                %%%%%%%%%%%%%
                k1f=ur*h;
                p1f=-(Grrr*(ur)^2+2*Grtr*ur*ut+Grtt*(ut)^2+Grpp*uphi^2)*h;
                l1f=ut*h;
                n1f=-(Gttt*(ut)^2+2*Gttr*ut*ur+Gtrr*(ur)^2+Gtpp*(uphi)^2)*h;
                q1f=uphi*h;
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                r1=r+k1f/2;
                t1=t+l1f/2;
                f=1-2*M/(r1);
                H=H0(s)+A(ss)*(t1);
                f1=2*M/(r1)^2;
                uphi=L/((r1)^2);
               [Grrr,Grtt,Grpp,Grtr,Gtrr,Gtpp,Gttr,Gttt]=Gammafunc(r1,t1,H,f,f1,AA);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                k2f=(ur+p1f/2)*h;
                p2f=-(Grrr*(ur+p1f/2)^2+2*Grtr*(ur+p1f/2)*(ut+n1f/2)+Grtt*(ut+n1f/2)^2+Grpp*uphi^2)*h;
                l2f=(ut+n1f/2)*h;
                n2f=-(Gttt*(ut+n1f/2)^2+2*Gttr*(ut+n1f/2)*(ur+p1f/2)+Gtrr*(ur+p1f/2)^2+Gtpp*(uphi)^2)*h;
                q2f=uphi*h;
                %%%%%%%%%%%%%%%%%%%
                r1=r+k2f/2;
                t1=t+l2f/2;
                f=1-2*M/(r1);
                H=H0(s)+A(ss)*(t1);
                f1=2*M/(r1)^2;
                uphi=L/((r1)^2);
                [Grrr,Grtt,Grpp,Grtr,Gtrr,Gtpp,Gttr,Gttt]=Gammafunc(r1,t1,H,f,f1,AA);
                %%%%%%%%%%%%%%%%%%%
                k3f=(ur+p2f/2)*h;
                p3f=-(Grrr*(ur+p2f/2)^2+2*Grtr*(ur+p2f/2)*(ut+n2f/2)+Grtt*(ut+n2f/2)^2+Grpp*uphi^2)*h;
                l3f=(ut+n2f/2)*h;
                n3f=-(Gttt*(ut+n2f/2)^2+2*Gttr*(ut+n2f/2)*(ur+p2f/2)+Gtrr*(ur+p2f/2)^2+Gtpp*(uphi)^2)*h;
                q3f=uphi*h;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                r1=r+k3f;
                t1=t+l3f;
                f=1-2*M/(r1);
                H=H0(s)+A(ss)*(t1);
                f1=2*M/(r1)^2;
                uphi=L/((r1)^2);
                [Grrr,Grtt,Grpp,Grtr,Gtrr,Gtpp,Gttr,Gttt]=Gammafunc(r1,t1,H,f,f1,AA);
                %%%%%%%%%%%%%%%%%%%
                k4f=(ur+p3f)*h;
                p4f=-(Grrr*(ur+p3f)^2+2*Grtr*(ur+p3f)*(ut+n3f)+Grtt*(ut+n3f)^2+Grpp*uphi^2)*h;
                l4f=(ut+n3f)*h;
                n4f=-(Gttt*(ut+n3f)^2+2*Gttr*(ut+n3f)*(ur+p3f)+Gtrr*(ur+p3f)^2+Gtpp*(uphi)^2)*h;
                q4f=uphi*h;
                %%%%%%%%%%%%%%%%%%%%
                kf=(k1f+2*k2f+2*k3f+k4f)/6 ; pf=(p1f+2*p2f+2*p3f+p4f)/6 ; lf=(l1f+2*l2f+2*l3f+l4f)/6 ; nf=(n1f+2*n2f+2*n3f+n4f)/6 ; qf=(q1f+2*q2f+2*q3f+q4f)/6 ;
                r=r+kf ; ur=ur+pf ; t=t+lf ; ut=ut+nf ; phi=phi+qf ; uphi=L/r^2;
                
                f=1-2*M/r;
                H=H0(s)+A(ss)*t;
                r1=0;
                t1=0;
                m=m+1;
                J=[kf,pf,lf,nf,qf];
                adapter=max(J);
                h=.0001*exp(-adapter^2);
            end
            tester(x)=r;
            rspeed=ur;
            pspeed=uphi;
            fpos=r;
            fff=f;
            tspeed=ut;
            tester1=tester(x);
            dr=RSL-tester1;
            phi;
            HH=H;
            numb=m-1;
             m=1;
             expectedut=[-(r*H*ur/sqrt(f))+sqrt([r^2*H^2*ur^2/f]+([ur^2/f]+[L^2/r^2])*(f-r^2*H^2))]/(f-r^2*H^2);
         end
        %%%%%%%%%%%%%%%%%%%%%%% 
        %bending angle calculations
             hubble(s,ss)=H0(s)+A(ss)*t;
             H0(s+1)=H0(s)+hub;
             expectation(s,ss)=expectedut;
             position(s,ss)=fpos;
             timespeed(s,ss)=tspeed;
             fvalue(s,ss)=fff;
             speed(s,ss)=rspeed;
             angularspeed(s,ss)=pspeed;
             psi=(180/pi)*atan(fpos*pspeed/rspeed);
             euclid(s,ss)=psi;
             W=[1,sqrt(fff)*(HH*fpos+sqrt(fff)),0,0];
             K=[tspeed,rspeed,0,pspeed];
             Ustat=[1/sqrt(fff-HH^2*fpos^2),0,0,0];
             Ucomov=[1/sqrt(fff),fpos*HH,0,0];
             metric=[-(fff-HH^2*fpos^2),-HH*fpos/sqrt(fff),0,0;-HH*fpos/sqrt(fff),1/fff,0,0;0,0,0,0;0,0,0,fpos^2];
             delta(s,ss)=abs(tspeed-expectedut)
             kdotw=0;usdotk=0;ucdotk=0;usdotw=0;ucdotw=0;i=1;
             while i<5
                  j=1;
                  while j<5
                      kdotw=kdotw+metric(i,j)*K(i)*W(j);
                      usdotk=usdotk+metric(i,j)*K(i)*Ustat(j);
                      ucdotk=ucdotk+metric(i,j)*K(i)*Ucomov(j);
                      usdotw=usdotw+metric(i,j)*W(i)*Ustat(j);
                      ucdotw=ucdotw+metric(i,j)*W(i)*Ucomov(j);
                      j=j+1;
                  end
                  i=i+1;
             end
             ctetmstat=(kdotw/(usdotk*usdotw))+1;
             ctetmcomov=(kdotw/(ucdotk*ucdotw))+1;
             tetstat=acos(ctetmstat)*180/pi;
             statpsi(s,ss)=tetstat;
             tetcomov=acos(ctetmcomov)*180/pi;
             comovpsi(s,ss)=tetcomov;
            b=1;
            tester1=0;
            numb;
            totnum(s,ss)=numb;
            s=s+1;
    end
    A(ss+1)=A(ss)+dA
    ss=ss+1;
    s=1;
end
 xlswrite('number04.xls',totnum);
 xlswrite('euclid04.xls',euclid);
 xlswrite('statpsi04.xls',statpsi);
 xlswrite('comovpsi04.xls',comovpsi);
 xlswrite('speed.xls',speed);
 xlswrite('tspeed.xls',timespeed);
 xlswrite('hubble.xls',hubble);
 xlswrite('position.xls',position);
 xlswrite('fvalue.xls',fvalue);
 xlswrite('angular.xls',angularspeed);
 xlswrite('delta.xls',delta);