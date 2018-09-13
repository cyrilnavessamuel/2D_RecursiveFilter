%Input Image
input = imread('castle.jpg');
%Convert Image to GrayScale
input_gray=rgb2gray(input);
[rows, columns ,depth] = size(input_gray);
x=input_gray;
%Declare coefficients

%Value of alpha
alpha =1.0;

%Value of K:
knum = (1-exp(-alpha))^2;
kden= 1+(2*alpha*exp(-alpha))-exp(-2*alpha);
k = knum/kden;

%Value of other coefficients:
a1=0; a5=k;
a2=1; a6= k* exp(-alpha)*(alpha-1);
a3=-1; a7=k*exp(-alpha)*(alpha+1);
a4=0; a8 = -k*exp(-2*alpha);

b1=2*exp(-alpha);b2=-exp(-2*alpha);

c1 = -((1-exp(-alpha))^2); c2=1;

%Start of Recursive Filtering Scheme
%Derivative of X

%Along Horizontal
%Step1: Calculate y1

%Boundary Conditions
%x(m,0) = 0;
%yl(m, 0), y1(m, -1) = 0;
y1=x;
for n=1:columns
    for m=1:rows
        if n==1
            y1(m,n)=a1*x(m,n)+a2*0+b1*0+b2*0;
        elseif n==2
            y1(m,n)=a1*x(m,n)+a2*x(m,n-1)+b1*y1(m,n-1)+b2*0;
        else
            y1(m,n)=a1*x(m,n)+a2*x(m,n-1)+b1*y1(m,n-1)+b2*y1(m,n-2);
        end
    end
end

%Step2: Calculate y2
%Boundary Conditions
%x(m,N+1),x(m,N+2)= 0
%y2(m,N+1),y2(m,N+2) = 0
y2=x;
for n=columns:-1:n
    for m=1:rows
        if n==columns
            y2(m,n)=a3*0+a4*0+b1*0+b2*0;
        else
            y2(m,n)=a3*x(m,n+1)+a4*x(m,n+2)+b1*y2(m,n+1)+b2*y2(m,n+2);
        end
    end
end
%Step3:Calculate R(m,n)
r=x;
for n=1:columns
    for m=1:rows
        r(m,n) =c1*(y1(m,n)+y2(m,n));
    end
end
%Along Vertical
%Step4: Calculate y3
%Boundary Conditions
%r(0, n) = 0;
%y1(0, n) ,y1(-1,n) = 0;
y3=x;
for m=1:rows
    for n=1:columns
        if m==1
            y3(m,n)=a5*r(m,n)+a6*0+b1*0+b2*0;
        elseif m==2
            y3(m,n)=a5*r(m,n)+a6*r(m-1,n)+b1*y1(m-1,n)+b2*0;
        else
            y3(m,n)=a5*r(m,n)+a6*r(m-1,n)+b1*y1(m-1,n)+b2*y1(m-2,n);
        end
    end
end
%Step5: Calculate y4
%Boundary Conditions
%r(M + 1, n) = r(M + 2, n) = 0 and
%y2(M + 1, n) = y2(M + 2, n) = 0
y4=x;
for m=rows:-1:m
    for n=1:columns
        if m==rows
            y4(m,n)=a7*0+a8*0+b1*0+b2*0;
        else
            y4(m,n)=a7*r(m+1,n)+a8*r(m+2,n)+b1*y2(m+1,n)+b2*y2(m+2,n);
        end
    end
end

%Step5: Calculate Xfinal
xfinal=x;
for n=1:columns
    for m=1:rows
        xfinal(m,n)=c2*(y3(m,n)+y4(m,n));
    end
end

%---------------------------------------------------------------------------------------------%

%Derivative of Y

%Reassigning value of other coefficients with swapping ai=ai+4 , c1=c2:
ya5=0; ya1=k;
ya6=1; ya2= k* exp(-alpha)*(alpha-1);
ya7=-1; ya3=k*exp(-alpha)*(alpha+1);
ya8=0; ya4 = -k*exp(-2*alpha);

yb1=2*exp(-alpha);yb2=-exp(-2*alpha);

yc1 = 1; yc2=-((1-exp(-alpha))^2);

%Along Horizontal
%Step1: Calculate yy1

%Boundary Conditions
%x(m,0) = 0;
%yl(m, 0), y1(m, -1) = 0;
yy1=x;
for n=1:columns
    for m=1:rows
        if n==1
            yy1(m,n)=ya1*x(m,n)+ya2*0+yb1*0+yb2*0;
        elseif n==2
            yy1(m,n)=ya1*x(m,n)+ya2*x(m,n-1)+yb1*yy1(m,n-1)+yb2*0;
        else
            yy1(m,n)=ya1*x(m,n)+ya2*x(m,n-1)+yb1*yy1(m,n-1)+yb2*yy1(m,n-2);
        end
    end
end

%Step2: Calculate y2
%Boundary Conditions
%x(m,N+1),x(m,N+2)= 0
%y2(m,N+1),y2(m,N+2) = 0
yy2=x;
for n=columns:-1:n
    for m=1:rows
        if n==columns
            yy2(m,n)=ya3*0+ya4*0+yb1*0+yb2*0;
        else
            yy2(m,n)=ya3*x(m,n+1)+ya4*x(m,n+2)+yb1*yy2(m,n+1)+yb2*yy2(m,n+2);
        end
    end
end

%Step3:Calculate R(m,n)
ry=x;
for n=1:columns
    for m=1:rows
        ry(m,n) =yc1*(yy1(m,n)+yy2(m,n));
    end
end

%Along Vertical
%Step4: Calculate y3
%Boundary Conditions
%r(0, n) = 0;
%y1(0, n) ,y1(-1,n) = 0;
yy3=x;
for m=1:rows
    for n=1:columns
        if m==1
            yy3(m,n)=ya5*ry(m,n)+ya6*0+yb1*0+yb2*0;
        elseif m==2
            yy3(m,n)=ya5*ry(m,n)+ya6*ry(m-1,n)+yb1*yy1(m-1,n)+yb2*0;
        else
            yy3(m,n)=ya5*ry(m,n)+ya6*ry(m-1,n)+yb1*yy1(m-1,n)+yb2*yy1(m-2,n);
        end
    end
end

%Step5: Calculate y4
%Boundary Conditions
%r(M + 1, n) = r(M + 2, n) = 0 and
%y2(M + 1, n) = y2(M + 2, n) = 0
yy4=x;
for m=rows:-1:m
    for n=1:columns
        if m==rows
            yy4(m,n)=ya7*0+ya8*0+yb1*0+yb2*0;
        else
            yy4(m,n)=ya7*ry(m+1,n)+ya8*ry(m+2,n)+yb1*yy2(m+1,n)+yb2*yy2(m+2,n);
        end
    end
end

%Step5: Calculate Yfinal
yfinal=x;
for n=1:columns
    for m=1:rows
        yfinal(m,n)=yc2*(yy3(m,n)+yy4(m,n));
    end
end


%Calculate Gradient Magnitude

xfinal=double(xfinal);
yfinal =double(yfinal);
im_edge=double(x);
angle = double(x);

for m=1:rows
    for n=1:columns
        im_edge(m,n)= ((xfinal(m,n)^2+yfinal(m,n)^2)).^(0.5);
        angle(m,n)=atan2(xfinal(m,n),yfinal(m,n));
    end
end

%Calculate Non Maxima Suppression and use threshold hysterisis to get the final output of the image edges.
O = nonMaximaSuppression(im_edge,angle);
O = threshold( O, 30);
colormap(gray);
imagesc(O);

%Calculate the threshold points which are greater than the specified threshold for us being 30
function O = threshold( I, t)
  O = ones(size(I,1),size(I,2));
  for i = 1 : size(I,1)
    for j = 1 : size(I,2)
        if(I(i,j) >= t)
            O(i,j) = 255;
        end
    end
  end
  
end

%Calculate Non Maximum Suppression where we calculate the direction of the gradient %
function O = nonMaximaSuppression( E , a)
ang =a;
    O = zeros(size(E,1),size(E,2));
    for row=2:size(E,1)-1
        for col=2:size(E,2)-1            
            if(localmaxima(E(row-1:row+1,col-1:col+1),ang(row,col)))
                O(row,col) = E(row,col);
            else
                O(row,col) = 0;
            end
        end
    end
    
 %Calculates by taking 9 points neighborhood with gradient values and angle.
function tf = localmaxima( p , a )

v1 = [2 2]' + [cos(a) sin(a)]';
v2 = [2 2]' + [cos(a+pi) sin(a+pi)]';

p1 = bipolar(p,v1);
p2 = bipolar(p,v2);

tf = ~(p(2,2) < p1 || p(2,2) < p2);

end
        
    mnO = min(min(O));
    mxO = max(max(O));
    Oout = zeros(size(O,1),size(O,2));
    
    for i=1:size(O,1)
        for j=1:size(O,2)
            Oout(i,j) = (255/(mxO-mnO+1))*(O(i,j)-mnO+1);
        end
    end

%Bilinear Interpolation of values in the direction specified plus specified angle
function val = bipolar( I, p )
  val = 1;
  val1 = -1;
  val2 = -1;
  
  if(size(p,2) == 1 && size(p,1) == 2)
      p = p';
  end
  
  a = floor(p);
  b = a + [0 1];
  c = a + [1 1];
  d = a + [1 0];
  
  t = p - a;
  
  if(a(2) > 0 && a(2) <= size(I,1) && a(1) > 0 && a(1) <= size(I,2))
    if(d(1) <= size(I,2))
      val1 = (1-t(1))*I(a(2),a(1)) + t(1)*I(d(2),d(1));
    else
      val1 = I(a(2),a(1));
    end
  else
    if(d(2) > 0 && d(2) <= size(I,1) && d(1) > 0 && d(1) <= size(I,2))
      val1 = I(d(2),d(1));
    end
  end
  if(b(2) > 0 && b(2) <= size(I,1) && b(1) > 0 && b(1) <= size(I,2))
    if(c(1) <= size(I,2))
      val2 = (1-t(1))*I(b(2),b(1)) + t(1)*I(c(2),c(1));
    else
      val2 = I(b(2),b(1));
    end
  else
    if(c(2) > 0 && c(2) <= size(I,1) && c(1) > 0 && c(1) <= size(I,2))
      val2 = I(c(2),c(1));
    end
  end
  if(val1 == -1 || val2 == -1)
    if(val1 == -1 && val2 == -1)
      val = 1;
    else
      if(val1 ~= -1 && val2 == -1)
        val = val1;
      else
          if(val1 == -1 && val2 ~= -1)
            val = val2;
          end
      end
    end
  else
    val = (1-t(2))*val1 + t(2)*val2;
  end
end
end






