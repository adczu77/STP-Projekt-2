clear;
Ko = 7.8;
To = 5;
T1 = 2.2;
T2 = 4.36;
Tp = 0.5;

Gs = tf(Ko,[T1*T2 T1+T2 1],'InputDelay',To);
Gz = c2d(Gs,Tp,'zoh');

% figure;
% hold on;
% step(Gs);
% step(Gz);
% title('Porównanie odpowiedzi skokowych transmitancji ciągłej i dyskretnej:');
% legend('G(s)','G(z)');
% hold off;
%print('stp_2_1wykres1.png','-dpng','-r400');
GzDenominator=Gz.Denominator{1,1};
GzNumerator=Gz.Numerator{1,1};


% figure;
% H = feedback(Gz,0.2594);
% step(H);
% xlim([0 100])
%print('stp_2_3wykres1.png','-dpng','-r400');

Kr = 0.6 * 0.2594;
Ti = 0.5 * 21.5;
Td = 0.12 * 21.5;

r2=(Kr*Td)/Tp;
r1=Kr*(Tp/(2*Ti)-2*(Td/Tp)-1);
r0=Kr*(Tp/(2*Ti)+Td/Tp+1);

kk=200;
u(1:12)=0; y(1:12)=0; 
yzad(1:12)=0; yzad(13:kk)=1; 
e(1:12)=0; 

for k=13:kk
 y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
 e(k)=yzad(k)-y(k); 
 u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1); 
end
p=(0:0.5:kk/2);
% figure; stairs(p(1:kk),u); 
% title('u'); xlabel('k'); ylabel('u');
% print('stp_2_4wykres1.png','-dpng','-r400');
% figure; stairs(p(1:kk),y); hold on; stairs(p(1:kk),yzad,':'); 
% title('yzad, y'); xlabel('k'); ylabel('y');
% legend('y','yzad');
% print('stp_2_4wykres2.png','-dpng','-r400');

uPID = u;
yPID = y;


% N=100
D=100;
N=100;
Nu=100;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad100=[];
ydmc100(1:D)=0;
yzad100(1:12)=0;
yzad100(13:D)=1;
dUp(1:D-1,1)=0;
udmc100(1:D)=0;
for k=13:D
    ydmc100(k)=-GzDenominator(2)*ydmc100(k-1)-GzDenominator(3)*ydmc100(k-2)+GzNumerator(2)*udmc100(k-11)+GzNumerator(3)*udmc100(k-12);
    dU=K*(ones(N,1).*yzad100(k)-ones(N,1).*ydmc100(k)-Mp*dUp);
    udmc100(k)=dU(1)+udmc100(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% figure; stairs(p(1:length(udmc100)),udmc100); 
% title('u'); xlabel('k'); 
% % print('stp_2_4wykres3.png','-dpng','-r400');
% figure; stairs(p(1:length(ydmc100)),ydmc100); hold on; stairs(p(1:length(yzad100)),yzad100,':');
% legend('y','yzad')
% title('yzad, y'); xlabel('k');
% % print('stp_2_4wykres4.png','-dpng','-r400');


% N=50
D=100;
N=50;
Nu=50;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad50=[];
ydmc50(1:D)=0;
yzad50(1:12)=0;
yzad50(13:D)=1;
dUp(1:D-1,1)=0;
udmc50(1:D)=0;
for k=13:D
    ydmc50(k)=-GzDenominator(2)*ydmc50(k-1)-GzDenominator(3)*ydmc50(k-2)+GzNumerator(2)*udmc50(k-11)+GzNumerator(3)*udmc50(k-12);
    dU=K*(ones(N,1).*yzad50(k)-ones(N,1).*ydmc50(k)-Mp*dUp);
    udmc50(k)=dU(1)+udmc50(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=25
D=100;
N=25;
Nu=25;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad25=[];
ydmc25(1:D)=0;
yzad25(1:12)=0;
yzad25(13:D)=1;
dUp(1:D-1,1)=0;
udmc25(1:D)=0;
for k=13:D
    ydmc25(k)=-GzDenominator(2)*ydmc25(k-1)-GzDenominator(3)*ydmc25(k-2)+GzNumerator(2)*udmc25(k-11)+GzNumerator(3)*udmc25(k-12);
    dU=K*(ones(N,1).*yzad25(k)-ones(N,1).*ydmc25(k)-Mp*dUp);
    udmc25(k)=dU(1)+udmc25(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);


% N=15
D=100;
N=15;
Nu=15;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15=[];
ydmc15(1:D)=0;
yzad15(1:12)=0;
yzad15(13:D)=1;
dUp(1:D-1,1)=0;
udmc15(1:D)=0;
for k=13:D
    ydmc15(k)=-GzDenominator(2)*ydmc15(k-1)-GzDenominator(3)*ydmc15(k-2)+GzNumerator(2)*udmc15(k-11)+GzNumerator(3)*udmc15(k-12);
    dU=K*(ones(N,1).*yzad15(k)-ones(N,1).*ydmc15(k)-Mp*dUp);
    udmc15(k)=dU(1)+udmc15(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=13
D=100;
N=13;
Nu=13;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad13=[];
ydmc13(1:D)=0;
yzad13(1:12)=0;
yzad13(13:D)=1;
dUp(1:D-1,1)=0;
udmc13(1:D)=0;
for k=13:D
    ydmc13(k)=-GzDenominator(2)*ydmc13(k-1)-GzDenominator(3)*ydmc13(k-2)+GzNumerator(2)*udmc13(k-11)+GzNumerator(3)*udmc13(k-12);
    dU=K*(ones(N,1).*yzad13(k)-ones(N,1).*ydmc13(k)-Mp*dUp);
    udmc13(k)=dU(1)+udmc13(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% figure;
% hold on;
% stairs(p(1:length(udmc100)),udmc100, '-');
% stairs(p(1:length(udmc50)),udmc50, '--');
% stairs(p(1:length(udmc25)),udmc25,':');
% stairs(p(1:length(udmc15)),udmc15,'-.');
% stairs(p(1:length(udmc13)),udmc13,'-');
% hold off;
% title('u'); xlabel('k');
% legend('N=100','N=50','N=25','N=15','N=13')
% % print('stp_2_5wykres1.png','-dpng','-r400');
% figure;
% hold on;
% stairs(p(1:length(ydmc100)),ydmc100,'-');
% stairs(p(1:length(ydmc50)),ydmc50,'--');
% stairs(p(1:length(ydmc25)),ydmc25,':');
% stairs(p(1:length(ydmc15)),ydmc15,'-.');
% stairs(p(1:length(ydmc13)),ydmc13,'-');
% stairs(p(1:length(yzad100)),yzad100,':');
% hold off;
% legend('y N=100','y N=50','y N=25','y N=15','y N=13','yzad')
% title('yzad, y'); xlabel('k');
% % print('stp_2_5wykres2.png','-dpng','-r400');

% N=15, Nu=15
D=100;
N=15;
Nu=15;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_15=[];
ydmc15_15(1:D)=0;
yzad15_15(1:12)=0;
yzad15_15(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_15(1:D)=0;
for k=13:D
    ydmc15_15(k)=-GzDenominator(2)*ydmc15_15(k-1)-GzDenominator(3)*ydmc15_15(k-2)+GzNumerator(2)*udmc15_15(k-11)+GzNumerator(3)*udmc15_15(k-12);
    dU=K*(ones(N,1).*yzad15_15(k)-ones(N,1).*ydmc15_15(k)-Mp*dUp);
    udmc15_15(k)=dU(1)+udmc15_15(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=10
D=100;
N=15;
Nu=10;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_10=[];
ydmc15_10(1:D)=0;
yzad15_10(1:12)=0;
yzad15_10(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_10(1:D)=0;
for k=13:D
    ydmc15_10(k)=-GzDenominator(2)*ydmc15_10(k-1)-GzDenominator(3)*ydmc15_10(k-2)+GzNumerator(2)*udmc15_10(k-11)+GzNumerator(3)*udmc15_10(k-12);
    dU=K*(ones(N,1).*yzad15_10(k)-ones(N,1).*ydmc15_10(k)-Mp*dUp);
    udmc15_10(k)=dU(1)+udmc15_10(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=5
D=100;
N=15;
Nu=5;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_5=[];
ydmc15_5(1:D)=0;
yzad15_5(1:12)=0;
yzad15_5(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_5(1:D)=0;
for k=13:D
    ydmc15_5(k)=-GzDenominator(2)*ydmc15_5(k-1)-GzDenominator(3)*ydmc15_5(k-2)+GzNumerator(2)*udmc15_5(k-11)+GzNumerator(3)*udmc15_5(k-12);
    dU=K*(ones(N,1).*yzad15_5(k)-ones(N,1).*ydmc15_5(k)-Mp*dUp);
    udmc15_5(k)=dU(1)+udmc15_5(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=3
D=100;
N=15;
Nu=3;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_3=[];
ydmc15_3(1:D)=0;
yzad15_3(1:12)=0;
yzad15_3(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_3(1:D)=0;
for k=13:D
    ydmc15_3(k)=-GzDenominator(2)*ydmc15_3(k-1)-GzDenominator(3)*ydmc15_3(k-2)+GzNumerator(2)*udmc15_3(k-11)+GzNumerator(3)*udmc15_3(k-12);
    dU=K*(ones(N,1).*yzad15_3(k)-ones(N,1).*ydmc15_3(k)-Mp*dUp);
    udmc15_3(k)=dU(1)+udmc15_3(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=1
D=100;
N=15;
Nu=1;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_1=[];
ydmc15_1(1:D)=0;
yzad15_1(1:12)=0;
yzad15_1(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_1(1:D)=0;
for k=13:D
    ydmc15_1(k)=-GzDenominator(2)*ydmc15_1(k-1)-GzDenominator(3)*ydmc15_1(k-2)+GzNumerator(2)*udmc15_1(k-11)+GzNumerator(3)*udmc15_1(k-12);
    dU=K*(ones(N,1).*yzad15_1(k)-ones(N,1).*ydmc15_1(k)-Mp*dUp);
    udmc15_1(k)=dU(1)+udmc15_1(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% figure;
% hold on;
% stairs(p(1:length(udmc15_15)),udmc15_15, '-');
% stairs(p(1:length(udmc15_10)),udmc15_10, '--');
% stairs(p(1:length(udmc15_5)),udmc15_5,':');
% stairs(p(1:length(udmc15_3)),udmc15_3,'-.');
% stairs(p(1:length(udmc15_1)),udmc15_1,'-');
% hold off;
% title('u'); xlabel('k');
% legend('Nu=15','Nu=10','Nu=5','Nu=3','Nu=1')
% print('stp_2_5wykres3.png','-dpng','-r400');
% figure;
% hold on;
% stairs(p(1:length(ydmc15_15)),ydmc15_15,'-');
% stairs(p(1:length(ydmc15_10)),ydmc15_10,'--');
% stairs(p(1:length(ydmc15_5)),ydmc15_5,':');
% stairs(p(1:length(ydmc15_3)),ydmc15_3,'-.');
% stairs(p(1:length(ydmc15_1)),ydmc15_1,'-');
% stairs(p(1:length(yzad100)),yzad100,':');
% hold off;
% legend('y Nu=15','y Nu=10','y Nu=5','y Nu=3','y Nu=1','yzad')
% title('yzad, y'); xlabel('k');
% print('stp_2_5wykres4.png','-dpng','-r400');

% N=15, Nu=1, lambda=1
D=100;
N=15;
Nu=1;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=1;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_1_1=[];
ydmc15_1_1(1:D)=0;
yzad15_1_1(1:12)=0;
yzad15_1_1(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_1_1(1:D)=0;
for k=13:D
    ydmc15_1_1(k)=-GzDenominator(2)*ydmc15_1_1(k-1)-GzDenominator(3)*ydmc15_1_1(k-2)+GzNumerator(2)*udmc15_1_1(k-11)+GzNumerator(3)*udmc15_1_1(k-12);
    dU=K*(ones(N,1).*yzad15_1_1(k)-ones(N,1).*ydmc15_1_1(k)-Mp*dUp);
    udmc15_1_1(k)=dU(1)+udmc15_1_1(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=1, lambda=5
D=100;
N=15;
Nu=1;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=5;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_1_5=[];
ydmc15_1_5(1:D)=0;
yzad15_1_5(1:12)=0;
yzad15_1_5(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_1_5(1:D)=0;
for k=13:D
    ydmc15_1_5(k)=-GzDenominator(2)*ydmc15_1_5(k-1)-GzDenominator(3)*ydmc15_1_5(k-2)+GzNumerator(2)*udmc15_1_5(k-11)+GzNumerator(3)*udmc15_1_5(k-12);
    dU=K*(ones(N,1).*yzad15_1_5(k)-ones(N,1).*ydmc15_1_5(k)-Mp*dUp);
    udmc15_1_5(k)=dU(1)+udmc15_1_5(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);

% N=15, Nu=1, lambda=10
D=100;
N=15;
Nu=1;

s=[];
y(1:12)=0;
u(1:12)=0;
u(13:500)=1;

for k=13:500
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=10;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_1_10=[];
ydmc15_1_10(1:D)=0;
yzad15_1_10(1:12)=0;
yzad15_1_10(13:D)=1;
dUp(1:D-1,1)=0;
udmc15_1_10(1:D)=0;
for k=13:D
    ydmc15_1_10(k)=-GzDenominator(2)*ydmc15_1_10(k-1)-GzDenominator(3)*ydmc15_1_10(k-2)+GzNumerator(2)*udmc15_1_10(k-11)+GzNumerator(3)*udmc15_1_10(k-12);
    dU=K*(ones(N,1).*yzad15_1_10(k)-ones(N,1).*ydmc15_1_10(k)-Mp*dUp);
    udmc15_1_10(k)=dU(1)+udmc15_1_10(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

p = (0:0.5:D/2);


% figure;
% hold on;
% stairs(p(1:length(udmc15_1_1)),udmc15_1_1, '-');
% stairs(p(1:length(udmc15_1_5)),udmc15_1_5, '--');
% stairs(p(1:length(udmc15_1_10)),udmc15_1_10,':');
% hold off;
% title('u'); xlabel('k');
% legend('lambda=1','lambda=5','lambda=10')
% print('stp_2_5wykres5.png','-dpng','-r400');
% figure;
% hold on;
% stairs(p(1:length(ydmc15_1_1)),ydmc15_1_1,'-');
% stairs(p(1:length(ydmc15_1_5)),ydmc15_1_5,'--');
% stairs(p(1:length(ydmc15_1_10)),ydmc15_1_10,':');
% stairs(p(1:length(yzad100)),yzad100,':');
% hold off;
% legend('y lambda=1','y lambda=5','y lambda=10','yzad')
% title('yzad, y'); xlabel('k');
% print('stp_2_5wykres6.png','-dpng','-r400');
% udmc15_1_5(101:200)=udmc15_1_5(100);
% ydmc15_1_5(101:200)=ydmc15_1_5(100);

% p = (0:0.5:D);
% figure;
% hold on;
% stairs(p(1:kk),udmc15_1_5, '-');
% stairs(p(1:kk),uPID);
% legend("DMC",'PID')
% hold off;
% title('u'); xlabel('k');
% ylim([-0.1 1])
% % print('stp_2_6wykres2.png','-dpng','-r400');
% figure;
% hold on;
% stairs(p(1:kk),ydmc15_1_5,'-');
% stairs(p(1:kk),yPID);
% stairs(p(1:length(yzad100)),yzad100,':');
% legend("DMC",'PID','yzad')
% hold off;
% title('yzad, y'); xlabel('k');
% % print('stp_2_6wykres3.png','-dpng','-r400');

clear;
To_Tonom=[1;1.1;1.2;1.3;1.4;1.5;1.6;1.7;1.8;1.9;2];
Ko = 7.8;
To = 5;
T1 = 2.2;
T2 = 4.36;
Tp = 0.5;

Gs = tf(7.8,[T1*T2 T1+T2 1],'InputDelay',5);
Gz = c2d(Gs,Tp,'zoh');

GzDenominator=Gz.Denominator{1,1};
GzNumerator=Gz.Numerator{1,1};
% Kr = 0.6 * 0.2594;
% Ti = 0.5 * 21.5;
% Td = 0.12 * 21.5;
% 
% r2=(Kr*Td)/Tp;
% r1=Kr*(Tp/(2*Ti)-2*(Td/Tp)-1);
% r0=Kr*(Tp/(2*Ti)+Td/Tp+1);
% 
% kk=200;
% u(1:22)=0; y(1:22)=0; 
% yzad(1:22)=0; yzad(23:kk)=1; 
% e(1:22)=0; 
% 
% for k=23:kk
%  y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-21)+GzNumerator(3)*u(k-22);
%  e(k)=yzad(k)-y(k); 
%  u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1); 
% end
K=[12.35;11.8;11.4;11.1;10.7;10.4;10.1;9.85;9.65;9.4;9.15];
K_Knom = K./Ko;
figure
plot(To_Tonom,K_Knom);
title('Krzywa stabilności PID:')
xlabel('To/To_{nom}')
ylabel('K/Knom')
print('stp_2_6wykres1.png','-dpng','-r400');

% N=15, Nu=1, lambda=5
D=2000;
N=15;
Nu=1;
p = (0:0.5:D/2);
s=[];
y(1:12)=0;
u(1:12)=0;
u(13:4000)=1;

for k=13:4000
    y(k)=-GzDenominator(2)*y(k-1)-GzDenominator(3)*y(k-2)+GzNumerator(2)*u(k-11)+GzNumerator(3)*u(k-12);
    s(k)=y(k);
end
s(1:13)=[];

Gs = tf(1.07,[T1*T2 T1+T2 1],'InputDelay',10);
Gz = c2d(Gs,Tp,'zoh');

GzDenominator=Gz.Denominator{1,1};
GzNumerator=Gz.Numerator{1,1};
Mp=zeros(N,D-1);
for i=1:D-1
    for j = 1:N
        Mp(j,i)=(s(j+i)-s(i));
    end
end

M=zeros(N,Nu);
for i=1:Nu
    M(i:N,i)=s(1:(N-i+1))';
end

YY=eye(N);
lambda=5;
A=eye(Nu).*lambda;
K=((M'*YY*M+A)^-1)*M'*YY;

yzad15_1_5=[];
ydmc15_1_5(1:D)=0;
yzad15_1_5(1:22)=0;
yzad15_1_5(23:D)=0.5;
dUp(1:D-1,1)=0;
udmc15_1_5(1:D)=0;
for k=23:D
    ydmc15_1_5(k)=-GzDenominator(2)*ydmc15_1_5(k-1)-GzDenominator(3)*ydmc15_1_5(k-2)+GzNumerator(2)*udmc15_1_5(k-21)+GzNumerator(3)*udmc15_1_5(k-22);
    dU=K*(ones(N,1).*yzad15_1_5(k)-ones(N,1).*ydmc15_1_5(k)-Mp*dUp);
    udmc15_1_5(k)=dU(1)+udmc15_1_5(k-1);
    dUp(2:D-1)=dUp(1:D-2);
    dUp(1)=dU(1);
end

% figure;
% hold on;
% stairs(p(1:length(udmc15_1_5)),udmc15_1_5, '-');
% hold off;
% title('u'); xlabel('k');
% figure;
% hold on;
% stairs(p(1:length(ydmc15_1_5)),ydmc15_1_5,'-');
% stairs(p(1:length(yzad15_1_5)),yzad15_1_5,':');
% hold off;
% title('yzad, y'); xlabel('k');
% print('stp_2_5wykres6.png','-dpng','-r400');

K=[15.4;15.2;14.2;7;3;1.7;1;0.75;0.7;0.8;1.07];
K_Knom = K./Ko;
figure
plot(To_Tonom,K_Knom);
title('Krzywa stabilności DMC:')
xlabel('To/To_{nom}')
ylabel('K/Knom')
print('stp_2_6wykres4.png','-dpng','-r400');