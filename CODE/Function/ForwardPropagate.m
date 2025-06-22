function model = ForwardPropagate(model)

model.X0data = model.Xgene  ; H0data = zeros(model.NumPath,model.NumSubj); H0data(model.IdxPath0,:) = (model.Dnet+model.Lnet)\model.Dnet*model.X0data;
model.X1data = H0data       ; H1data = zeros(model.NumPath,model.NumSubj); H1data(model.IdxPath1,:) = (model.Spath1.*model.Mpath1)*model.X1data      ;
model.X2data = H0data+H1data; H2data = zeros(model.NumPath,model.NumSubj); H2data(model.IdxPath2,:) = (model.Spath2.*model.Mpath2)*model.X2data      ;
model.X3data = H0data+H2data; H3data = zeros(model.NumPath,model.NumSubj); H3data(model.IdxPath3,:) = (model.Spath3.*model.Mpath3)*model.X3data      ;
model.X4data = H0data+H3data; H4data = zeros(model.NumPath,model.NumSubj); H4data(model.IdxPath4,:) = (model.Spath4.*model.Mpath4)*model.X4data      ;

model.Hdata = H0data + H1data + H2data + H3data + H4data;

model.Pabt = 1./(1+exp((-1)*model.Babt'*model.Hdata));
model.Pgfa = 1./(1+exp((-1)*model.Bgfa'*model.Hdata));
model.Pnfl = 1./(1+exp((-1)*model.Bnfl'*model.Hdata));
model.Ptau = 1./(1+exp((-1)*model.Btau'*model.Hdata));