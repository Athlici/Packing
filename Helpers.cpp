PhiPolygon regularPolygon(int n){
    double phi=2*M_PI/n;
    vector<point> v(n);
    for(int i=0;i<n;i++)
        v[i]=point(sin(phi*i),cos(phi*i));
    return PhiPolygon(v);
}

vector<double> boundCircMod(SmartPtr<IpoptApplication> app, PhiPolygon A){
    Objective* obj = new FirstVar();
    vector<var> vars = {var(0,2e19),var(-2e19,2e19),var(-2e19,2e19)};

    Scale f = Scale(0);
    Translate g = Translate(1);
    PhiCircCompl C = PhiCircCompl(circle(point(0,0),1));
    PhiFunc* phi = phiFunc(C,f,A,g);

    double x[] = {1,0,0};
    while(phi->eval(x)<0)
        x[0]*=2;

    bool newseg;
    do{
        SmartPtr<dNLP> nlp = new dNLP(obj,vars,phi->getIneqs(x),x);
        app->OptimizeTNLP(nlp);
        newseg = !(phi->getIneqs(x) == phi->getIneqs(nlp->res));
        for(int i=0;i<3;i++) x[i] = nlp->res[i];
    }while(newseg);
    
    return vector<double>(x,x+3);
}
