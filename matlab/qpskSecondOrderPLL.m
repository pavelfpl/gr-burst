function [ y ] = qpskSecondOrderPLL( x, alfa, beta )

I_sig = zeros(1,length(x));
Q_sig = zeros(1,length(x));
y = zeros(1,length(x));
  
buff_1=0;
buff_3=0; 

for n=1:1:length(x)
   I_sig(n)=real(x(n))*cos(buff_3)-imag(x(n))*sin(buff_3);
   Q_sig(n)=real(x(n))*sin(buff_3)+imag(x(n))*cos(buff_3); 
   
   y(n) = I_sig(n) + 1j*Q_sig(n);

   if I_sig(n) > 0
      dI=1;
   else
      dI=-1; 
   end
    
   if I_sig(n) == 0
      dI=0;
   end
       
   if Q_sig(n) > 0
      dQ=1;
   else
      dQ=-1; 
   end
    
    if Q_sig(n) == 0
      dQ=0;
   end
   
   q=Q_sig(n)*dI-I_sig(n)*dQ;
    
   buff_1=buff_1+q*beta;
   buff_2=buff_1+q*alfa;
   buff_3=buff_3-buff_2;  % - no carrier - only phase 2*pi*fc/fs;
    
   if(buff_3>2*pi)
      buff_3=buff_3-2*pi;
   end
end