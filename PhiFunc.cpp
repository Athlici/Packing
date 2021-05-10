class PhiFuncPrim;

tuple<point,point> moveLine(point,point,double);

class PhiFunc{
    static double sign(double x){return ((x>=0) ? 1 : -1);}
  public:
    virtual double eval(const double* x) = 0;
    virtual vector<PhiFuncPrim*> getIneqs(const double* x) = 0;
    virtual vector<int> getIndices(const double* x) = 0;
    virtual vector<string> print(const double* x) = 0;

    static PhiFunc* phiFunc(point,point,RotTrans,point ,RotTrans);
    static PhiFunc* phiFunc(point,point,RotTrans,circle,RotTrans);
    static PhiFunc* phiFunc(circle,RotTrans,point ,RotTrans);
    static PhiFunc* phiFunc(circle,RotTrans,circle,RotTrans);
    static PhiFunc* phiFunc(circle,Scale,point ,RotTrans);
    static PhiFunc* phiFunc(circle,Scale,circle,RotTrans);
    static PhiFunc* phiFunc(point,point,Scale,point ,RotTrans);
    static PhiFunc* phiFunc(point,point,Scale,circle,RotTrans);
    static PhiFunc* phiFunc(circle,RotTrans,circle,point,RotTrans,double);
    static PhiFunc* phiFunc(circle,Scale   ,circle,point,RotTrans,double);
};

class PhiFuncPrim: public PhiFunc{
  public:
//0->no var, i/-i->(i-1)th var (cos)
//constant -> (0,0), linear -> (0,i), quadratic -> (i,j)
    virtual vector<tuple<int,int,double>> getF() = 0;
    virtual void getD1(const double* x,double* res) = 0;
    virtual vector<int> getD1ind() = 0;
    virtual void getD2(const double* x,double* res) = 0;
    virtual vector<tuple<int,int>> getD2ind() = 0;
    vector<PhiFuncPrim*> getIneqs (const double* x){return {this};}
    vector<int> getIndices (const double* x){return {};}
};

class PhiFuncNode: public PhiFunc{
    bool sense; //0 -> max,1 -> min
    vector<PhiFunc*> nodes;

    int argmax(const double* x){
        int ind = 0;
        double val = nodes[0]->eval(x);
        for(int i=1;i<nodes.size();i++){
            double tmp = nodes[i]->eval(x);
            if(tmp>val){
                val = tmp;
                ind = i;
            }
        }
        return ind;
    }

  public:
    PhiFuncNode(bool s,vector<PhiFunc*> n) : sense(s),nodes(n) {};

    double eval(const double* x){
        double res = nodes[0]->eval(x);
        for(int i=1;i<nodes.size();i++){
            double tmp = nodes[i]->eval(x);
            if(sense ^ (tmp>res))
                res=tmp;
        }
        return res;
    }

    vector<PhiFuncPrim*> getIneqs(const double* x){
        if(sense){
            vector<PhiFuncPrim*> res = {};
            for(int i=0;i<nodes.size();i++){
                vector<PhiFuncPrim*> tmp = nodes[i]->getIneqs(x);
                res.insert(res.end(),tmp.begin(),tmp.end());
            }
            return res;
        } else
            return nodes[argmax(x)]->getIneqs(x);
    }

    vector<int> getIndices(const double* x){
        if(sense){
            vector<int> res = {};
            for(int i=0;i<nodes.size();i++){
                vector<int> tmp = nodes[i]->getIndices(x);
                res.insert(res.end(),tmp.begin(),tmp.end());
            }
            return res;
        } else {
            int ind = argmax(x);
            vector<int> res = nodes[ind]->getIndices(x);
            res.push_back(ind);
            return res;
        }
    }

    vector<string> print(const double* x){
        vector<string> out = {string(sense?"Min":"Max") + "-Node : " + std::to_string(eval(x))};
        for(int i=0;i<nodes.size();i++){
            vector<string> tmp = nodes[i]->print(x);
            out.push_back("-" + tmp[0]);
            for(int j=1;j<tmp.size();j++)
                out.push_back("|" + tmp[j]);
        }
        return out;
    }
};

class PhiFuncLnScClRt : public PhiFuncPrim{
    RowVector3d trs;
    RowVector2d rot;
    double o;
    int i,j;

  public:
    PhiFuncLnScClRt(point l0,point l1,Scale f,circle c,RotTrans g){
        double x=c.p.x, y=c.p.y, dx=l0.x-l1.x, dy=l0.y-l1.y, det=l0.x*l1.y-l1.x*l0.y;
        trs << -det,-dy, dx; rot << dx*x+dy*y,-dy*x+dx*y;
        double n=1/sqrt(dx*dx+dy*dy);
        trs *= n; rot *= n;
        o = -c.r;
        i = f.ind; j = g.ind;
    }

    double eval(const double* x){
        Vector3d f(x[i],x[j],x[j+1]);
        Vector2d g(sin(x[j+2]),cos(x[j+2]));
        return trs.dot(f)+rot.dot(g)+o;
    }

    vector<tuple<int,int,double>> getF(){
        return {{0,i+1,trs[0]},{0,j+1,trs[1]},{0,j+2,trs[2]},{0,j+3,rot[0]},{0,-j-3,rot[1]},{0,0,o}};
    }

    vector<string> print(const double* x){
        return {"LnScClRt : " + std::to_string(eval(x))};
    }

    void getD1(const double* x,double* res){
        Vector2d g(cos(x[j+2]),-sin(x[j+2]));
        res[0] = trs[0]; res[1] = trs[1]; res[2] = trs[2];
        res[3] = rot*g;
    }

    vector<int> getD1ind(){
        return {i,j,j+1,j+2};
    }

    void getD2(const double* x,double* res){
        Vector2d g(sin(x[j+2]),cos(x[j+2]));
        res[0] = -rot*g;
    }

    vector<tuple<int,int>> getD2ind(){
        return {{j+2,j+2}};
    }
};

class PhiFuncCcScClRt : public PhiFuncPrim{
    double cd(circle c){  //Circle Distance
        return c.r*c.r-c.p.x*c.p.x-c.p.y*c.p.y;
    }
    Matrix3d A;
    double m,b;
    int i,j;

  public:
    PhiFuncCcScClRt(circle cc,Scale f,circle c,RotTrans g){
        double xc=cc.p.x, yc=cc.p.y, px=c.p.x, py=c.p.y;
        A << cd(cc)/2, xc, yc,
          px*yc-py*xc, py,-px,
          px*xc+py*yc,-px,-py;
        A *= 2; m = -2*cc.r*c.r; b = cd(c);
        i = f.ind; j = g.ind;
    }

    double eval(const double* x){
        RowVector3d f(x[i],sin(x[j+2]),cos(x[j+2]));
        Vector3d g(x[i],x[j],x[j+1]);
        return f*A*g-g[1]*g[1]-g[2]*g[2]+m*x[i]+b;
    }

    vector<tuple<int,int,double>> getF(){
        return {{ i+1,i+1,A(0,0)},{ i+1,j+1,A(0,1)},{ i+1,j+2,A(0,2)},
                { j+3,i+1,A(1,0)},{ j+3,j+1,A(1,1)},{ j+3,j+2,A(1,2)},
                {-j-3,i+1,A(2,0)},{-j-3,j+1,A(2,1)},{-j-3,j+2,A(2,2)},
                {j+1,j+1,-1},{j+2,j+2,-1},{0,i+1,m},{0,0,b}};
    }

    vector<string> print(const double* x){
        return {"CcScClRt : " + std::to_string(eval(x))};
    }

    vector<int> getD1ind(){
        return {i,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        RowVector3d f(x[i],sin(x[j+2]),cos(x[j+2])),fA=f*A;
        Vector3d g(x[i],x[j],x[j+1]),Ag=A*g;
        res[0] = fA[0]+Ag[0]+m;
        res[1] = fA[1]-2*g[1];
        res[2] = fA[2]-2*g[2];
        res[3] = f[2]*Ag[1]-f[1]*Ag[2];
    }

    vector<tuple<int,int>> getD2ind(){
        if(i<j)
          return {{i,i},{j,j},{j+1,j+1},{j,j+2},{j+1,j+2},{j+2,j+2},{i,j},{i,j+1},{i,j+2}};
        else
          return {{i,i},{j,j},{j+1,j+1},{j,j+2},{j+1,j+2},{j+2,j+2},{j,i},{j+1,i},{j+2,i}};
    }

    void getD2(const double* x,double* res){
        double s = sin(x[j+2]),c = cos(x[j+2]);
        Vector3d g(x[i],x[j],x[j+1]),Ag=A*g;
        res[0] = 2*A(0,0); res[1] = -2; res[2] = -2;
        res[3] = A(1,1)*c-A(2,1)*s;
        res[4] = A(1,2)*c-A(2,2)*s;
        res[5] = -Ag[1]*s-Ag[2]*c;
        res[6] = A(0,1); res[7] = A(0,2);
        res[8] = A(1,0)*c-A(2,0)*s;
    }
};

class PhiFuncHCcScClRt : public PhiFuncPrim{
    Eigen::Matrix<double,3,2> A;
    double b;
    int i,j;

  public:
    PhiFuncHCcScClRt(circle cc,Scale f,circle c,point p,RotTrans g,double s){
        double ox=cc.p.x, oy=cc.p.y, cx=c.p.x, cy=c.p.y;
        A << ox*(cx-p.x)+oy*(cy-p.y),ox*(cy-p.y)+oy*(p.x-cx),p.x-cx,p.y-cy,p.y-cy,cx-p.x;
        b = s*(cx*p.y-cy*p.x); A*=s;
        i = f.ind; j = g.ind;
    }

    double eval(const double* x){
        RowVector3d f(x[i],x[j],x[j+1]);
        Vector2d g(sin(x[j+2]),cos(x[j+2]));
        return f*A*g+b;
    }

    vector<tuple<int,int,double>> getF(){
        return {{i+1,j+3,A(0,0)},{i+1,-j-3,A(0,1)},
                {j+1,j+3,A(1,0)},{j+1,-j-3,A(1,1)},
                {j+2,j+3,A(2,0)},{j+2,-j-3,A(2,1)},{0,0,b}};
    }

    vector<string> print(const double* x){
        return {"HCcScClRt : " + std::to_string(eval(x))};
    }

    vector<int> getD1ind(){
        return {i,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        Vector2d g(sin(x[j+2]),cos(x[j+2])),dg(g[1],-g[0]);
        Vector3d Ag=A*g;
        RowVector3d f(x[i],x[j],x[j+1]);
        res[0] = Ag[0]; res[1] = Ag[1]; res[2] = Ag[2];
        res[3] = f*A*dg;
    }

    vector<tuple<int,int>> getD2ind(){
        if(i<j)
          return {{i,j+2},{j,j+2},{j+1,j+2},{j+2,j+2}};
        else
          return {{j+2,i},{j,j+2},{j+1,j+2},{j+2,j+2}};
    }

    void getD2(const double* x,double* res){
        Vector2d g(sin(x[j+2]),cos(x[j+2])),dg(g[1],-g[0]);
        Vector3d Adg=A*dg;
        RowVector3d f(x[i],x[j],x[j+1]);
        res[0] = Adg[0]; res[1] = Adg[1]; res[2] = Adg[2];
        res[3] = -f*A*g;
    }
};

class PhiFuncRtRtMdfg : public PhiFuncPrim{
    Matrix2d A,B;
    double c;
    int i,j;

  public:
    PhiFuncRtRtMdfg(Matrix2d Ai, Matrix2d Bi, double ci, RotTrans f, RotTrans g) :
        A(Ai),B(Bi),c(ci),i(f.ind),j(g.ind) {};

    double eval(const double* x){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2])),y(x[i]-x[j],x[i+1]-x[j+1]);
        Vector2d f2(sin(x[j+2]),cos(x[j+2]));
        return (y*A+f1*B)*f2+c;
    }

    vector<tuple<int,int,double>> getF(){
        return {{i+1,j+3,A(0,0)},{j+1,j+3,-A(0,0)},{i+1,-j-3,A(0,1)},{j+1,-j-3,-A(0,1)},
                {i+2,j+3,A(1,0)},{j+2,j+3,-A(1,0)},{i+2,-j-3,A(1,1)},{j+2,-j-3,-A(1,1)},
                {i+3,j+3,B(0,0)},{i+3,-j-3,B(0,1)},{-i-3,j+3,B(1,0)},{-i-3,-j-3,B(1,1)},{0,0,c}};
    }

    vector<string> print(const double* x){
        return {"RtRtMdfg : " + std::to_string(eval(x)) 
            + " " + std::to_string(i) + " " + std::to_string(j)};
    }

    vector<int> getD1ind(){
        return {i,i+1,i+2,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        Vector2d    f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        RowVector2d y(x[i]-x[j],x[i+1]-x[j+1]);
        Vector2d    Af2 = A*f2;
        double tmp[] = {Af2[0],Af2[1],f1d*B*f2,-Af2[0],-Af2[1],(y*A+f1*B)*f2d};
        for(int i=0;i<6;i++) res[i] = tmp[i];
    }

    vector<tuple<int,int>> getD2ind(){
        if(i<j)
            return {{i+2,i+2},{i,j+2},{i+1,j+2},{i+2,j+2},{j,j+2},{j+1,j+2},{j+2,j+2}};
        else
            return {{i+2,i+2},{j+2,i},{j+2,i+1},{j+2,i+2},{j,j+2},{j+1,j+2},{j+2,j+2}};
    }

    void getD2(const double* x,double* res){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        Vector2d    f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        RowVector2d y(x[i]-x[j],x[i+1]-x[j+1]);
        Vector2d    Af2d = A*f2d;
        double f1Bf2 = -f1*B*f2;
        double tmp[] {f1Bf2,Af2d[0],Af2d[1],f1d*B*f2d,-Af2d[0],-Af2d[1],f1Bf2-y*A*f2};
        for(int i=0;i<7;i++) res[i] = tmp[i];
    }
};

class PhiFuncClRtPtRt : public PhiFuncPrim{
    Matrix2d A,R1,R2;
    double b,s;
    int i,j;

  public:
    PhiFuncClRtPtRt(circle c, RotTrans f, point p, RotTrans g, double o = 1){
        double cx = c.p.x, cy = c.p.y, px = p.x, py = p.y;
        A  << -cx*px-cy*py, cy*px-cx*py, cx*py-cy*px, -cx*px-cy*py;
        R1 <<  py,-px,-px,-py; R2 << -cy, cx, cx, cy;
        b = cx*cx+px*px+cy*cy+py*py-c.r*c.r;
        s = o; A*=2*o; R1*=2*o; R2*=2*o; b*=o;
        i = f.ind; j = g.ind;
    }
    
    double eval(const double* x){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2]));
        Vector2d    f2(sin(x[j+2]),cos(x[j+2]));
        Vector2d     y(x[i]-x[j],x[i+1]-x[j+1]);
        return y.dot(s*y+R1*f2)+f1*(R2*y+A*f2)+b;
    }

    vector<tuple<int,int,double>> getF(){
        return {{i+1,i+1,s},{i+1,j+1,-2*s},{j+1,j+1,s},{i+2,i+2,s},{i+2,j+2,-2*s},{j+2,j+2,s},
                {i+1,j+3,R1(0,0)},{j+1,j+3,-R1(0,0)},{i+1,-j-3,R1(0,1)},{j+1,-j-3,-R1(0,1)},
                {i+2,j+3,R1(1,0)},{j+2,j+3,-R1(1,0)},{i+2,-j-3,R1(1,1)},{j+2,-j-3,-R1(1,1)},
                { i+3,i+1,R2(0,0)},{ i+3,j+1,-R2(0,0)},{ i+3,i+2,R2(0,1)},{ i+3,j+2,-R2(0,1)},
                {-i-3,i+1,R2(1,0)},{-i-3,i+1,-R2(1,0)},{-i-3,i+2,R2(1,1)},{-i-3,j+2,-R2(1,1)},
                {i+3,j+3,A(0,0)},{i+3,-j-3,A(0,1)},{-i-3,j+3,A(1,0)},{-i-3,-j-3,A(1,1)},{0,0,b}};
    }

    vector<string> print(const double* x){
        return {"ClRtPtRt : " + std::to_string(eval(x))};
    }

    vector<int> getD1ind(){
        return {i,i+1,i+2,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        Vector2d    f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        Vector2d     y(x[i]-x[j],x[i+1]-x[j+1]);
        Vector2d dy = 2*s*y+(f1*R1).transpose()+R2*f2;
        double tmp[] {dy[0],dy[1],f1d*(R1*y+A*f2),-dy[0],-dy[1],y.dot(R2*f2d)+f1*A*f2d};
        for(int i=0;i<6;i++) res[i] = tmp[i];
    }

    vector<tuple<int,int>> getD2ind(){
        if(i<j)
            return {{i,i},{i+1,i+1},{i,i+2},{i+1,i+2},{i+2,i+2},
                    {j,j},{j+1,j+1},{j,j+2},{j+1,j+2},{j+2,j+2},
                    {i,j},{i+1,j+1},{i,j+2},{i+1,j+2},{i+2,j},{i+2,j},{i+2,j+2}};
        else
            return {{i,i},{i+1,i+1},{i,i+2},{i+1,i+2},{i+2,i+2},
                    {j,j},{j+1,j+1},{j,j+2},{j+1,j+2},{j+2,j+2},
                    {j,i},{j+1,i+1},{j,i+2},{j+1,i+2},{j+2,i},{j+2,i},{j+2,i+2}};
    }

    void getD2(const double* x,double* res){
        RowVector2d f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        Vector2d    f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        Vector2d     y(x[i]-x[j],x[i+1]-x[j+1]);
        RowVector2d dyf1 = f1d*R1;
        Vector2d    dyf2 = R2*f2d;
        double fAf= -f1*A*f2;
        double tmp[] {2*s, 2*s, dyf1[0], dyf1[1],         -f1*R1*y+fAf,
                      2*s, 2*s, dyf2[0], dyf2[1],    -y.dot(R2*f2)+fAf,
                     -2*s,-2*s,-dyf1[0],-dyf1[1],-dyf2[0],-dyf2[1],fAf};
    }
};

PhiFunc* PhiFunc::phiFunc(point l0,point l1,RotTrans f,point p,RotTrans g){
    Matrix2d A,B;
    double px=p.x, py=p.y, dx=l1.x-l0.x, dy=l1.y-l0.y, c;
    double n=1/sqrt(dx*dx+dy*dy); dx *= n; dy *= n;
    A << dx, dy, dy, -dx;
    B << dy*px-dx*py, -dx*px-dy*py, dx*px+dy*py, dy*px-dx*py;
    c = n*(l1.x*l0.y-l0.x*l1.y);
    return new PhiFuncRtRtMdfg(A,B,c,g,f);
}

PhiFunc* PhiFunc::phiFunc(point p0,point p1,RotTrans f,circle c,RotTrans g){
    tuple<point,point> tmp = moveLine(p0,p1,c.r);
    return PhiFunc::phiFunc(get<0>(tmp),get<1>(tmp),f,c.p,g);
}

PhiFunc* PhiFunc::phiFunc(circle c,RotTrans f,point p,RotTrans g){
    return new PhiFuncClRtPtRt(c,f,p,g,sign(c.r));
}
//If c1<0 then it should be -c1>c2
PhiFunc* PhiFunc::phiFunc(circle c1,RotTrans f,circle c2,RotTrans g){
    return new PhiFuncClRtPtRt(circle(c1.p,c1.r+c2.r),f,c2.p,g,sign(c1.r));
}

//OK for scaling, otherwise negative circle radius determines complement
PhiFunc* PhiFunc::phiFunc(circle cc,Scale f,point p,RotTrans g){
    return new PhiFuncCcScClRt(cc,f,circle(p,0),g);
}

PhiFunc* PhiFunc::phiFunc(circle cc,Scale f,circle c,RotTrans g){
    return new PhiFuncCcScClRt(cc,f,c,g);
}

PhiFunc* PhiFunc::phiFunc(point p0,point p1,Scale f,point p,RotTrans g){
    return new PhiFuncLnScClRt(p0,p1,f,circle(p,0),g);
}

PhiFunc* PhiFunc::phiFunc(point p0,point p1,Scale f,circle c,RotTrans g){
    return new PhiFuncLnScClRt(p0,p1,f,c,g);
}

PhiFunc* PhiFunc::phiFunc(circle cc,RotTrans f,circle pc,point p,RotTrans g,double s){
    Matrix2d A,B;
    double ox=cc.p.x, oy=cc.p.y, dx=pc.p.x-p.x, dy=pc.p.y-p.y, c;
    A << dx, dy, dy, -dx;
    B << ox*dy-oy*dx, -ox*dy-oy*dx, ox*dy+oy*dx, ox*dy-oy*dx;
    c = pc.p.x*p.y-pc.p.y*p.x;
    A*= s; B*= s; c*= s;
    return new PhiFuncRtRtMdfg(A,B,c,f,g);
}

PhiFunc* PhiFunc::phiFunc(circle cc,Scale    f,circle pc,point p,RotTrans g,double s){
    return new PhiFuncHCcScClRt(cc,f,pc,p,g,s);
}
