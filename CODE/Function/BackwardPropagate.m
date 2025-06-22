function model = BackwardPropagate(model)

dLdZabt = (1/model.NumTrain)*(model.Pabt(model.IdxTrain)-model.Yabt(model.IdxTrain));
dLdZgfa = (1/model.NumTrain)*(model.Pgfa(model.IdxTrain)-model.Ygfa(model.IdxTrain));
dLdZnfl = (1/model.NumTrain)*(model.Pnfl(model.IdxTrain)-model.Ynfl(model.IdxTrain));
dLdZtau = (1/model.NumTrain)*(model.Ptau(model.IdxTrain)-model.Ytau(model.IdxTrain));

dLdBabt = model.Hdata(:,model.IdxTrain)*dLdZabt'; gBabt = dLdBabt + 2*model.RegGamma*model.Babt;
dLdBgfa = model.Hdata(:,model.IdxTrain)*dLdZgfa'; gBgfa = dLdBgfa + 2*model.RegGamma*model.Bgfa;
dLdBnfl = model.Hdata(:,model.IdxTrain)*dLdZnfl'; gBnfl = dLdBnfl + 2*model.RegGamma*model.Bnfl;
dLdBtau = model.Hdata(:,model.IdxTrain)*dLdZtau'; gBtau = dLdBtau + 2*model.RegGamma*model.Btau;

dLdHdata = (model.Babt*dLdZabt)+(model.Bgfa*dLdZgfa)+(model.Bnfl*dLdZnfl)+(model.Btau*dLdZtau);

dLdHpath4 = dLdHdata(model.IdxPath4,:);
dLdSpath4 = dLdHpath4*model.X4data(:,model.IdxTrain)';
dLdWpath4 = dLdSpath4.*(model.Spath4.*(1-model.Spath4)).*model.Mpath4;
gWpath4 = dLdWpath4 + 2*model.RegGamma*model.Wpath4;

dLdHpath3 = dLdHdata+((model.Spath4.*model.Mpath4)'*dLdHpath4);
dLdHpath3 = dLdHpath3(model.IdxPath3,:);
dLdSpath3 = dLdHpath3*model.X3data(:,model.IdxTrain)';
dLdWpath3 = dLdSpath3.*(model.Spath3.*(1-model.Spath3)).*model.Mpath3;
gWpath3 = dLdWpath3 + 2*model.RegGamma*model.Wpath3;

dLdHpath2 = dLdHdata+((model.Spath3.*model.Mpath3)'*dLdHpath3);
dLdHpath2 = dLdHpath2(model.IdxPath2,:);
dLdSpath2 = dLdHpath2*model.X2data(:,model.IdxTrain)';
dLdWpath2 = dLdSpath2.*(model.Spath2.*(1-model.Spath2)).*model.Mpath2;
gWpath2 = dLdWpath2 + 2*model.RegGamma*model.Wpath2;

dLdHpath1 = dLdHdata+((model.Spath2.*model.Mpath2)'*dLdHpath2);
dLdHpath1 = dLdHpath1(model.IdxPath1,:);
dLdSpath1 = dLdHpath1*model.X1data(:,model.IdxTrain)';
dLdWpath1 = dLdSpath1.*(model.Spath1.*(1-model.Spath1)).*model.Mpath1;
gWpath1 = dLdWpath1 + 2*model.RegGamma*model.Wpath1;

dLdHpath0 = dLdHdata+((model.Spath1.*model.Mpath1)'*dLdHpath1);
dLdHpath0 = dLdHpath0(model.IdxPath0,:);
dLdUnet = dLdHpath0*model.X0data(:,model.IdxTrain)'/(model.Dnet+model.Lnet)*(eye(model.NumGene)-((model.Dnet+model.Lnet)\model.Dnet)).*eye(model.NumGene)*ones(1,model.NumGene)';
gUnet = dLdUnet + 2*model.RegGamma*model.Unet;

model.Gradient = [gUnet(:);gWpath1(:);gWpath2(:);gWpath3(:);gWpath4(:);gBabt(:);gBgfa(:);gBnfl(:);gBtau(:)];