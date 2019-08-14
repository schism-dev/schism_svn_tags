function fg = hgrid2fg(hgrid);

fid = fopen(hgrid);
if fid<0
    disp('Fuhgeddaboudit . . .')
    fg = [];
    return
end
out = fgets(fid);
out = fscanf(fid,'%i %i',2);
elems = out(1); nodes = out(2);

out = fscanf(fid,'%i %f %f %f\n',[4,nodes]);
x=out(2,:)';y=out(3,:)';z=out(4,:)';

out = fscanf(fid,'%i %i %i %i %i\n',[5,elems]);
e=out(3:5,:)';

num = 0;
ob_num = fscanf(fid,'%i',1); fgets(fid);
if ~isempty(ob_num) % Then this is a Lat/Lon gridfile
%  _______________

ob_ntot = fscanf(fid,'%i',1); fgets(fid);
for ob = 1:ob_num
    ob_nod = fscanf(fid,'%i',1); fgets(fid);
    node = fscanf(fid,'%i',1);
    for nod = 2:ob_nod
        num = num + 1;
        bnd(num,1) = node;
        bnd(num,2) = fscanf(fid,'%i',1);
        node = bnd(num,2);
    end
end

l_num = fscanf(fid,'%i',1); fgets(fid);
l_ntot = fscanf(fid,'%i',1); fgets(fid);
isles = 0;
for l = 1:l_num
    l_nod = fscanf(fid,'%i ',1); outstr = fgets(fid);
    if strfind(outstr,'island')|strcmp(outstr(1),'1')
        isle = 1;
        isles = isles+1;
    else
        isle = 0;
    end
    if isle==0
        node = fscanf(fid,'%i',1);
        for nod = 2:l_nod
            num = num + 1;
            bnd(num,1) = node;
            bnd(num,2) = fscanf(fid,'%i',1);
            node = bnd(num,2);
        end
    else
        node = fscanf(fid,'%i',1);
        for nod = 2:l_nod
            num = num + 1;
            bnd(num,1) = node;
            bnd(num,2) = fscanf(fid,'%i',1);
            node = bnd(num,2);
        end
        if node~=bnd(num-(l_nod-2),1)
            num = num + 1;
            bnd(num,1) = node;
            bnd(num,2) = bnd(num-(l_nod-1),1);
        end
    end
end
    
fg.bnd = bnd;

ind = 1;
taken = ones(size(bnd))==0;
bndpth(ind) = NaN;
ind = ind+1;
next = 1;
while any(~taken)
    status = 'open';
    bndpth(ind) = bnd(next,1);
    taken(next,1) = true;
    ind = ind+1;
    bndpth(ind) = bnd(next,2);
    taken(next,2) = true;
    ind = ind+1;
    while strcmp(status,'open')
        next = find(bnd(:,1)==bndpth(ind-1) & ~taken(:,2));
        if ~isempty(next)
            bndpth(ind) = bnd(next,2);
            ind = ind + 1;
            taken(next,2) = true;
        else
            next = find(bnd(:,2)==bndpth(ind-1) & ~taken(:,1));
            if length(next)>1
                next
            end
            if ~isempty(next)
                bndpth(ind) = bnd(next,1);
                ind = ind + 1;
                taken(next,1) = true;
            else
                bndpth(ind) = NaN;
                ind = ind+1;
                status = 'closed';
                avail = find(~taken(:,1));
                if any(avail)
                    next = avail(1);
                end
            end
        end
    end
end
 
fg.bndpth = bndpth';

%_______________
else
    fg.bndpth = 'NA';
    fg.bnd = 'NA';
end

fg.name = 'selfe_grid';
fg.x = x;
fg.y = y;
fg.z = -z;
fg.e = e;

x1 = x(e(:,1)); x2 = x(e(:,2)); x3 = x(e(:,3)); 
y1 = y(e(:,1)); y2 = y(e(:,2)); y3 = y(e(:,3)); 
% Grant's back-o'-envelope method (works!)
%m = (y2-y1) ./ (x2-x1); % slope 1-2
%b = y2 - m.*x2;
%c = y3 + (1./m).*x3;
%X = (c-b)./(m+1./m);
%Y = m.*X + b;
%base = sqrt((x(e(:,2)) - x(e(:,1))).^2 + (y(e(:,2)) - y(e(:,1))).^2);
%height = sqrt((x(e(:,3)) - X).^2 + (y(e(:,3)) - Y).^2);
%ar = 1/2 .* base .* height;
%Heron's method???
%K=sqrt(s(s-a)(s-b)(s-c)) s=a/2+b/2+c/2
% Joseph's method (tri_area)
ar=((x1-x3).*(y2-y3)-(x2-x3).*(y1-y3))./2;

fg.ar = ar;