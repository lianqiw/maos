clear a b c
%Test non-cell version
a{1}=rand(4,4);
b{1}=rand(4,4);
%try
%c=aolib('dcellmm2',a{1},b{1},'nn');
%c-a{1}*b{1}
%end
c=aolib('dcellmm2',a,b,'nn');
c{1}-a{1}*b{1}
c=aolib('dcellmm2',a,b,'nt');
c{1}-a{1}*b{1}'
c=aolib('dcellmm2',a,b,'tt');
c{1}-a{1}'*b{1}'
c=aolib('dcellmm2',a,b,'tn');
c{1}-a{1}'*b{1}

a{1}=rand(4,4);
b{1}=sprand(4,4,.3);
c=aolib('dcellmm2',a,b,'nn');
c{1}-a{1}*b{1}
c=aolib('dcellmm2',a,b,'nt');
c{1}-a{1}*b{1}'
c=aolib('dcellmm2',a,b,'tt');
c{1}-a{1}'*b{1}'
c=aolib('dcellmm2',a,b,'tn');
c{1}-a{1}'*b{1}


a{1}=sprand(4,4,.3);
b{1}=rand(4,4);
c=aolib('dcellmm2',a,b,'nn');
c{1}-a{1}*b{1}
c=aolib('dcellmm2',a,b,'nt');
c{1}-a{1}*b{1}'
c=aolib('dcellmm2',a,b,'tt');
c{1}-a{1}'*b{1}'
c=aolib('dcellmm2',a,b,'tn');
c{1}-a{1}'*b{1}

a{1}=sprand(4,4,.3);
b{1}=sprand(4,4,.3);
c1=aolib('dcellmm2',a,b,'nn');
c1{1}-a{1}*b{1}
c2=aolib('dcellmm2',a,b,'nt');
c2{1}-a{1}*b{1}'
c3=aolib('dcellmm2',a,b,'tt');
c3{1}-a{1}'*b{1}'
c4=aolib('dcellmm2',a,b,'tn');
c4{1}-a{1}'*b{1}


%Test ell version
clear a b c

a=mat2cell(rand(16,16),[8 8],[8 8]);
b=mat2cell(rand(16,16),[8 8],[8 8]);
for i=1:4
    if i==2
        for j=1:4
            a{j}=sprand(8,8,0.4);
        end
    end
    if i==3
        for j=1:4
            b{j}=sprand(8,8,0.4);
        end
    end
    c=aolib('dcellmm2',a,b,'nn');
    norm(full(cell2mat(c))-full(cell2mat(a)*cell2mat(b)))
    c=aolib('dcellmm2',a,b,'nt');
    norm(full(cell2mat(c))-full(cell2mat(a)*cell2mat(b)'))
    c=aolib('dcellmm2',a,b,'tt');
    norm(full(cell2mat(c))-full(cell2mat(a)'*cell2mat(b)'))
    c=aolib('dcellmm2',a,b,'tn');
    norm(full(cell2mat(c))-full(cell2mat(a)'*cell2mat(b)))

end
