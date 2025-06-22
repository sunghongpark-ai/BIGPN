function model = ParamReshape(model)

sizeUnet = model.sizeParam(1,:)  ; start = 1                     ; finish = prod(sizeUnet)         ; model.Unet = reshape(model.WeightParam(start:finish),sizeUnet)    ;
sizeWpath1 = model.sizeParam(2,:); start = start+prod(sizeUnet)  ; finish = finish+prod(sizeWpath1); model.Wpath1 = reshape(model.WeightParam(start:finish),sizeWpath1);
sizeWpath2 = model.sizeParam(3,:); start = start+prod(sizeWpath1); finish = finish+prod(sizeWpath2); model.Wpath2 = reshape(model.WeightParam(start:finish),sizeWpath2);
sizeWpath3 = model.sizeParam(4,:); start = start+prod(sizeWpath2); finish = finish+prod(sizeWpath3); model.Wpath3 = reshape(model.WeightParam(start:finish),sizeWpath3);
sizeWpath4 = model.sizeParam(5,:); start = start+prod(sizeWpath3); finish = finish+prod(sizeWpath4); model.Wpath4 = reshape(model.WeightParam(start:finish),sizeWpath4);
sizeBabt = model.sizeParam(6,:)  ; start = start+prod(sizeWpath4); finish = finish+prod(sizeBabt)  ; model.Babt = reshape(model.WeightParam(start:finish),sizeBabt)    ;
sizeBgfa = model.sizeParam(7,:)  ; start = start+prod(sizeBabt)  ; finish = finish+prod(sizeBgfa)  ; model.Bgfa = reshape(model.WeightParam(start:finish),sizeBgfa)    ;
sizeBnfl = model.sizeParam(8,:)  ; start = start+prod(sizeBgfa)  ; finish = finish+prod(sizeBnfl)  ; model.Bnfl = reshape(model.WeightParam(start:finish),sizeBnfl)    ;
sizeBtau = model.sizeParam(9,:)  ; start = start+prod(sizeBnfl)  ; finish = finish+prod(sizeBtau)  ; model.Btau = reshape(model.WeightParam(start:finish),sizeBtau)    ;

model.Dnet = diag(model.Unet);
model.Spath1 = 1./(1+(exp((-1)*model.Wpath1)));
model.Spath2 = 1./(1+(exp((-1)*model.Wpath2)));
model.Spath3 = 1./(1+(exp((-1)*model.Wpath3)));
model.Spath4 = 1./(1+(exp((-1)*model.Wpath4)));