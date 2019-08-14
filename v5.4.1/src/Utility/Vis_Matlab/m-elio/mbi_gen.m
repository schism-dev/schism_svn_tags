function mbi = mbi_gen(st,fn,sStep)
%function mbi = mbi_gen(start,fin,sStep)
% generate member index (mbi) from start date to fin date
% star/fin as fin.yyyy, fin.ww, fin.dd, fin.nn 
%
%Sergey Frolov, september 2004


nc=1;
for y = st.yyyy:fn.yyyy,
	for w = st.ww:fn.ww,
        lld =1; uld=7;
        if w == st.ww
            lld = st.dd;
        end
        if w == fn.ww
            uld = fn.dd;
        end 
        
        for d = lld:uld,            
            lln =1; uln =96;
            if d == st.dd
                lln = st.nn;
            end
            if d == fn.dd
                uln = fn.nn;
            end
            
            for n = [lln:sStep:uln],
				mbi{nc}.yyyy  = y; 
                mbi{nc}.ww    = w; 
                mbi{nc}.dd    = d; 
                mbi{nc}.n     = n;
                nc=nc+1;
			end
            
		end
	end
end

%disp(['n: ' num2str(nc-1)])
