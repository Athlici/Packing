PhiCompObj* regularPolygon(int n){
    double phi=2*M_PI/n;
    vector<point> v(n);
    for(int i=0;i<n;i++)
        v[i]=point(sin(phi*i),cos(phi*i));
    return new PhiPolygon(v);
}

circle boundCircMod(SmartPtr<IpoptApplication> app, PhiCompObj* A){
    Objective* obj = new FirstVar();
    vector<var> vars = {var(0,2e19),var(-2e19,2e19),var(-2e19,2e19),var(0,0)};

    Scale f = Scale(0);
    RotTrans g = RotTrans(1);
    PhiInfObj* C = new PhiCircCompl(circle(point(0,0),1));
    PhiFunc* phi = phiFunc(C,f,A,g);

    double x[] = {1,0,0,0};
    while(phi->eval(x)<0)
        x[0]*=2;

    bool newseg;
    do{
        SmartPtr<dNLP> nlp = new dNLP(obj,vars,phi->getIneqs(x),x);
        app->OptimizeTNLP(nlp);
        newseg = phi->getIneqs(x) != phi->getIneqs(nlp->res);
        for(int i=0;i<3;i++) x[i] = nlp->res[i];
    }while(newseg);
    
    return circle(point(x[1],x[2]),x[0]);
}

vector<point> circlePack(vector<circle> cs){
    int n = cs.size();
    vector<tuple<circle,int>> chain(n);
    chain[0] = {circle(point(-cs[0].r,0),cs[0].r),0};
    chain[1] = {circle(point( cs[1].r,0),cs[1].r),1};

    for(int i=2;i<n;i++){
        bool unplaced = true;
        for(int j=0;unplaced;j++){
            circle c1 = get<0>(chain[j]), c2 = get<0>(chain[(j+1)%i]);
            double x1=c1.p.x, y1=c1.p.y, r1=c1.r, x2=c2.p.x, y2=c2.p.y, r2=c2.r, r=cs[i].r;
            double dx=x1-x2, dy=y1-y2, dr=r1-r2, dx2=dx*dx, dy2=dy*dy, sr=r1+r2+2*r;
            double n1 = sqrt((dx2+dy2-dr*dr)*(sr*sr-dx2-dy2)), n2 = dx2+dy2;
            point pos = point((x1+x2-(n1*dy+dr*sr*dx)/n2)/2,(y1+y2+(n1*dx-dr*sr*dy)/n2)/2);

            bool placeable = true;
            for(int k=0;k<i&&placeable;k++){
                circle tmp = get<0>(chain[k]);
                dx = pos.x-tmp.p.x, dy = pos.y-tmp.p.y, sr = r+tmp.r;
                placeable = dx*dx+dy*dy-sr*sr > -1e-8;
            }
            if(placeable){
                for(int k=i-1;k>j;k--)
                    chain[k+1]=chain[k];
                chain[j+1] = {circle(pos,r),i};
                unplaced = false;
            }
        }
    }

    vector<point> res(n);
    for(int i=0;i<n;i++){
        tuple<circle,int> tmp = chain[i];
        res[get<1>(tmp)] = get<0>(tmp).p;
    }
    return res;
}
