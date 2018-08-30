t=cputime;
m=1508;%number of rows in A
n=2071; %number of coloumns in A
c=556; % random number=K, rank of the matrix A
P=rand(1,n); % for calculating coloumn Sums
A3=A.^2; % A3 has squared values of A
sum1=sum(sum(A3)); % Sum of all the squared values of A
for i=1:n 
A2=A(:,i).^2;
P(:,i)=sum(A2)/sum1; %coloumn sum/net sum
end;
Columns=zeros(c,1); %for storing what coloumns are being retrieved
C=rand(m,c); %C is of size mxc 
for i=1:c 
r = rand;
prob = P.'; %this code will pick up coloumns with greater probability
%j = sum(r >= cumsum([0, P])); 
j = sum(r>=norm(P,'fro'));
%j=mod(j,2071);
if(j==0)
    j=j+1;
end;

Columns(i,1)=j; %writing coloumns numbers that are being retrieved for C 
C(:,i)=A(:,j)/(sqrt(c*P(j))); %store the whole coloumn of A indexed by j
end;

Rows=zeros(c,1);  %for storing rows that are being retrieved
P=rand(1,m); % for calculating row Sums
A3=A.^2; % A3 has squared values of A
sum1=sum(sum(A3)); % Sum of all the squared values of A
for i=1:m 
A2=A(i,:).^2; %row sum/net sum
P(:,i)=sum(A2)/sum1;
end;
R=rand(c,n); %R is of size cxn 
for i=1:c 
r = rand;
prob = P.';%this code will pick up rows with greater probability
%j = sum(r >= cumsum([0, P]));
j = sum(r>=norm(P,'fro'));
%j=mod(j,1508);
if(j==0)
    j=j+1;
end;
Rows(i,1)=j; %writing row numbers that are being retrieved for R
R(i,:)=A(j,:)/(sqrt(c*P(j)));
end;
W=rand(c,c);
Columns=sort(Columns);
Rows=sort(Rows);
for i=1:c
for j=1:c
       W(i,j) =A(Rows(i,1),Columns(j,1));
end;
end;
%W is the required matrix
%W=zeros(2,2);
%W(1,1)=5;
%W(2,2)=5;
[X,Z,Y]=svd(W);
BB=X*Z*Y';
Z=pinv(Z);
Z=Z';
Z=Z.^2;
%this is finding sigma squared;
U=Y*Z*(X');
%finding U for Sigma;
Res=C*U*R; %this is the CUR matrix required.
%FB=Res-A;
%FB=FB*FB.';
FB=norm(Res-A,'fro'); %FB=2,315, FB=2,196
e=cputime-t;







