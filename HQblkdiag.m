function D=HQblkdiag(A,cnt)

[n,m]=size(A);
k=length(cnt)-1;
D=zeros(n,k*m);
for i=1:k
    D(cnt(i)+1:cnt(i+1),(i-1)*m+1:i*m)=A(cnt(i)+1:cnt(i+1),:);
end