function [ob_st, ob_mn]  = ob_apply_ob2state(st, ob, stateInfo, dbInfo, stVar, obsType,dryFlag,normFlag)
%[ob_st, ob_mn]  = ob_apply_ob2state(st, ob, stateInfo, dbInfo, stVar, obsType, [dryFlag=0, normFlag=1])
% applies observation strategy "ob" to the state "st"
% stateInfo describes the structure of the state varaible
% dbInfo, describes layout of xyz grid (as in elcirc binary file) 
% stVar defines which varaible (e|s|t|u|v) in the state is observed
% obsType = {xy|xyz} defines wich observation opertors will be applied
% dryFlag =[0|1] if 1 output has special treatment of dry nodes
% normFlag=[0|1] if 1 output is in natural units, if 0 output is in units of the state vector
%
% ob_st - observed state, state is streemlined into a vector
% ob_mn - dimension of the ob_st.
%
% sergey frolov, december 2004

if nargin < 8
   normFlag = 1;
end
if nargin < 7
   dryFlag = 0;
end

if stVar    == 'e'
    strIdxVar   = stateInfo.idx.str.e;
    stVarVar    = stateInfo.var.e;
elseif stVar    == 's'
    strIdxVar   = stateInfo.idx.str.s;
    stVarVar    = stateInfo.var.s;
elseif stVar    == 't'
    strIdxVar   = stateInfo.idx.str.t;
    stVarVar    = stateInfo.var.t;
elseif stVar    == 'u'
    strIdxVar   = stateInfo.idx.str.u;
    stVarVar    = stateInfo.var.u;
elseif stVar    == 'v'
    strIdxVar   = stateInfo.idx.str.v;
    stVarVar    = stateInfo.var.v;
else
    error (['unknown stVar ' stVar])
end

% normFlag==0 - no normaliztion by stVarVar is needed
if normFlag==0
   stVarVar=1;
end

if stVar    == 'e'
    if strcmp(obsType,'xyz')
        error(['obsType and stVar are incompatible'])
    end
    stTmp = st(eval(strIdxVar))*stVarVar;
else 
    stTmp = st(eval(strIdxVar));
        %reorganize state vector by levles, 
        % vertical strucutre described in dbInfo.fheaders.s will work for t u v as well
    stTmp = map_eb2hts_nan(dbInfo.fheaders.s,stTmp);
    stTmp = stTmp*stVarVar;
end
    
if strcmp(obsType,'xy')
    ob_st = ob.xy.H*stTmp;
elseif strcmp(obsType,'xyz')
    ob_st = diag(ob.z.H*(ob.xy.H*stTmp)');
else
    error(['unknown obsType' obsType])
end

%in case dry nodes need to be considered
if dryFlag==1 & stVar ~= 'e'
	%get transect of elevation first
    e_IdxVar   	= stateInfo.idx.str.e;
    e_stVarVar  = stateInfo.var.e;
    e_stTmp 	= st(eval(e_IdxVar))*e_stVarVar;
    e_ob_st 	= ob.xy.H*e_stTmp;
%    e_ob_st_a1	= repmat(e_ob_st,1,62); %next line generalizes beyond corie, yet to be tested
    e_ob_st_a1  = repmat(e_ob_st,1,ob.gr.vgrid.nLevels);
    e_ob_st_a2  = repmat(ob.gr.vgrid.zLevel-ob.gr.vgrid.zMsl,1,size(ob.xy.H,1))';
    ob_st(find(e_ob_st_a1<e_ob_st_a2))=nan;
end

ob_mn	= size(ob_st);
ob_st 	=ob_st(:);


