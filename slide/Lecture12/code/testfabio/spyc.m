function out=spyc(A,range,startc,id)
%SPYC.m _______________________________ Visualize the sparsity structure of A:
%       USE spyc to add colors for different scales (automatically)
% Usage:  spyc(A,range,startc, id) % 3 optional parameters
%    1/   range specifies the values at which spy check matrix A (see id)
%    2/   startc is for starting color (either a color say 'g' or integer say 4)
%    3/   If id=2, then spyc uses range as verbatim otherwise as 10^range
% Examples:   x=0:pi/63:pi;  A=gallery('circul',sin(x)*100);
%             spyc(A)
%             spyc(A, -3:3)
%             spyc(A, -3:2, 'r')    % same as  spyc(A, -3:2, 7) 
%             syc(A, [0.1 0.7 1.3 1.9],4)   
%             spyc(A,[1 3 5 55 88 99],1,2)
%
%@Copyright (2001-2002) Ke Chen (V3)  SEE also    spym (3D)

if nargin<1, help spyc; return; else
  if nargin>1, rlen=length(range); else, rlen=0; end
end
if nargin<2 | rlen==0
  valu=max(max(abs(A)));
  if valu>10
   range = (10:-1:0);
  elseif valu>1
   range = (1:-1:-2);
  else
   range = (-1:-1:-10);
  end
   id=1; end
%---------------------
   clist ='wbycmgr'; 
   wlist ='wbycmgr'; XL = []; R = []; VA=[];
if nargin<3
   c = 1; startc=1;
elseif isstr(startc)
   c = find(clist==startc);
else
   c=startc;
end
if nargin<4, id=1; end

  rlen=length(range); % Data flip if necessary
  if rlen>1 
    if range(2)>range(1)
    range=range(rlen:-1:1);
    end
  end
  nnz_check=1; nnz_lastc=0; 

for k=range
   valu = 10^k; if id==2, valu=k; end  % valu = k;
  check=abs(A)>valu; nnz_check=nnz(check);
  if nnz_check>0 & nnz_check~=nnz_lastc  % repetition
      XL =[ XL nnz_check];
      VA =[ VA valu];
    nnz_lastc=nnz_check;
  end % if
end % for

  rlen=length(XL);
  if c>7,c=1; end; symbol='.'; count=0;

if rlen<1 %-----------------------------Rescue plot
  Note_A_max_abs = full( max(max(abs(A))) )
  Note_A_min_abs = full( min(min(abs(A))) )
  close;spyc(A,(Note_A_max_abs-eps)./(10.^(3:-2:1)), 6,2)
  return
else %----------------------------------Actual plot

  XL=XL(rlen:-1:1);
  VA=VA(rlen:-1:1);
 for active=1:rlen
  k=VA(active); check=abs(A)>k; 
      spy(check,[clist(c) symbol]); hold on
    Ex = 0;
    if abs(k)>0, Ex=round( log(k)/log(10) ); end
  if k==round(k)
    if Ex==0
    R1=sprintf('%s %d,',wlist(c),k); R = [R R1];
    else
    s=sprintf('10^{%d}',Ex); 
    R1=sprintf('%s %s,',wlist(c),s); R = [R R1];
    end
  else
    if Ex==0
    R1=sprintf('%s %2.1f,',wlist(c),k); R = [R R1];
    else
    s=sprintf('10^{%d}',Ex); 
    R1=sprintf('%s %s,',wlist(c),s); R = [R R1];
    end
  end
  c=c+1; if c>7,c=1; end
  count=count+1; if count>7, symbol='o'; end
                if count>14, symbol='d'; end
                if count>21, symbol='+'; end
                if count>28, symbol='x'; end
                if count>35, symbol='*'; end
 end % for
end % IF %------------------------------------------

  rlen=length(R);  R=R(1:rlen-1); % Last comma
  hold off
  s=sprintf('Matrix %d x %d using colors/levels  ',size(A));
  s=sprintf('%s %s',s,R); h=title(s);
  x=sprintf('%d, ',XL); rlen=length(x); x=x(1:rlen-2);
  LL=xlabel(['spyc.m - NNZ for each color = ' x], 'Color','r');
  ylabel('|A|>values (as specified)', 'color','b');
  if length(s)>140 % too long
      zz=get(h,'Fontsize'); set(h,'Fontsize',zz*0.45);
     zz=get(LL,'Fontsize');set(LL,'Fontsize',zz*0.45);
  elseif length(s)>70 % too long
      zz=get(h,'Fontsize'); set(h,'Fontsize',zz*0.75);
     zz=get(LL,'Fontsize');set(LL,'Fontsize',zz*0.75);
  end
% clen=length(VA);
  set(gca,'Color',[0 0 0])
  set(gca,'YAxisLocation','right')
