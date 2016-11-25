% writes the root system parameters in a way CRootBox can read it

p = completeParameters(p);

fname = 'test.rparam';
fid = fopen([pwd '/' fname],'w');

for i  = 1 : length(p);
    p1=p{i};
    fprintf(fid,['# Parameter set for type ' '(by Rootbox)\n']);
    fprintf(fid,'type\t%g\n',p1.type);
    name = p1.name; % delete spaces
    name(name==' ')=[];    
    fprintf(fid,'name\t%s\n',name);
    fprintf(fid,'lb\t%g\t%g\n',p1.lb(1),p1.lb(2));
    fprintf(fid,'la\t%g\t%g\n',p1.la(1),p1.la(2));
    fprintf(fid,'ln\t%g\t%g\n',p1.ln(1),p1.ln(2));
    fprintf(fid,'nob\t%g\t%g\n',p1.nob(1),p1.nob(2));
    fprintf(fid,'r\t%g\t%g\n',p1.r(1),p1.r(2));
    fprintf(fid,'a\t%g\t%g0\n',p1.a(1),p1.a(2));
    fprintf(fid,'color\t%g\t%g\t%g\n',p1.color(1),p1.color(2),p1.color(3));
    fprintf(fid,'tropism\t%g\t%g\t%g\n',p1.tropism(1),p1.tropism(2),p1.tropism(3));
    fprintf(fid,'dx\t%g\n',p1.dx(1));
    fprintf(fid,'successors');
    fprintf(fid,'\t%g',size(p1.successor,1));
    for i = 1 : size(p1.successor,1)
        fprintf(fid,'\t%g',p1.successor(i,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'successorP');
    fprintf(fid,'\t%g',size(p1.successor,1));
    for i = 1 : size(p1.successor,1)
        fprintf(fid,'\t%g',p1.successor(i,2));
    end
    fprintf(fid,'\n');
    fprintf(fid,'theta\t%g\t%g\n',p1.theta(1),p1.theta(2));
    if isinf(p1.rlt(1))
        fprintf(fid,'rlt\t%g\t%g\n',1e9,0); % replace Inf by large number
    else
        fprintf(fid,'rlt\t%g\t%g\n',p1.rlt(1),p1.rlt(2));
    end
    fprintf(fid,'gf\t%g\n',p1.gf(1));
end

fclose(fid);