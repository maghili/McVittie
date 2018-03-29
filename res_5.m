%initial conditions
r(1)=1;%radial coordinate
t(1)=0;%time coordinate
ur(1)=0;%radial velocity
l=.5;%angular momentum
a=0;%acceleration
H(1)=.000;%Hubble parameter
M=.1;%Mass of the black hole
h=.001;%step size
p=.005;
lambda(1)=0;
phi(1)=pi/2;% initial position in azimuthal coordinate
u=1;
for j=1:1:100
    H(1)=H(1)+.001;
    fh(j)=H(1);
    ut(1)=l/(r(1)*((1-2*M/r(1))-((H(1)+a*t(1))*r(1))^2)^(1/2));
    for i=1:1:100000%Runge-Kutta order 4
        dt=ut(i)*h;
        dr=ur(i)*h;
        dut=((-((H(1)+a*t(i))*r(i)*(r(i)*(H(1)+a*t(i))^2-(M/r(i)^2)))...
            *ut(i)^2/((1-(2*M/r(i)))^(1/2))-((2*M/r(i)^2)-2*((H(1)+a*t(i))^2)...
            *r(i))*ut(i)*ur(i)/(1-2*M/r(i))-(H(1)+a*t(i))*ur(i)^2/(1-2*M/r(i))^(3/2)-(H(1)+a*t(i))*l^2/(r(i)^2*(1-2*M/r(i))^(1/2))))*h;
        wrr(i)=(M/r(i)^2-(H(1)+a*t(i))^2*r(i));
        wrt(i)=r(i)*(H(1)+a*t(i))*(-r(i)*(H(1)+a*t(i))^2+M/r(i)^2)/((1-2*M/r(i))^(1/2));
        wtt(i)=(1/2)*(r(i)^2*(H(1)+a*t(i))^2*(2*r(i)*(H(1)+a*t(i))^2-2*M/r(i)^2)+(1-2*M/r(i))*(-2*r(i)*(H(1)+a*t(i))^2+2*M/r(i)^2)-2*a*(1-2*M/r(i))^(1/2));
        wpp(i)=(-(1-2*M/r(i))+r(i)^2*(H(1)+a*t(i))^2)*l^2/r(i)^3;
        dur=(wrr(i)*ur(i)^2-wrt(i)*ur(i)*ut(i)-wtt(i)*ut(i)^2-wpp(i))*h;
        dphi=(l/r(i)^2)*h;
        uphi(i)=(l/r(i)^2);
        phi(i+1)=phi(i)+dphi;
        r(i+1)=r(i)+dr;
        t(i+1)=t(i)+dt;
        ut(i+1)=ut(i)+dut;
        ur(i+1)=ur(i)+dur;
        lambda(i+1)=lambda(1)+i*h;
    end
    q=1;
    while phi(q)<pi
        q=q+1;
    end
    psi(j)=atan(r(q-1)*uphi(q-1)/ur(q-1));
end
fig=figure()
plot(fh,psi)
xlabel('H_0')
ylabel('\psi (rad)')
print(fig,'H const','-dpng')
%clear