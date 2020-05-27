class gQP{
  private:
    GRBModel model;
    vector<GRBVar> vars,vcos;
    bool* dual;
    int n;

    static int sign(int i){
        return (i!=0?(i>0?1:-1):0);
    }

    GRBQuadExpr parseExpr(vector<tuple<int,int,double>> f){
      GRBQuadExpr expr = GRBQuadExpr();
      expr.clear();
//      for(auto [i,j,d] : f)
      for(int l=0;l<f.size();l++){
        int i = get<0>(f[l]),j = get<1>(f[l]);
        double d = get<2>(f[l]);
        if(d>0.001 || d<-0.001)
        switch(3*sign(i)+sign(j)){
          case  4: expr.addTerm(d,vars[ i-1],vars[ j-1]); break;
          case  2: expr.addTerm(d,vars[ i-1],vcos[-j-1]); break;
          case -2: expr.addTerm(d,vcos[-i-1],vars[ j-1]); break;
          case -4: expr.addTerm(d,vcos[-i-1],vcos[-j-1]); break;
          case  1: expr.addTerm(d,vars[ j-1]);            break;
          case -1: expr.addTerm(d,vcos[-j-1]);            break;
          case  0: expr.addConstant(d);                   break;
        }
      
      }

      return expr;
    }

  public:
    double* res;

    gQP(GRBEnv env,Objective* fi,vector<var> varsi,vector<PhiFuncPrim*> phix/*,const double* x0i*/) : model(GRBModel(env)){
//      model = GRBModel(env);
      model.set(GRB_IntParam_NonConvex, 2); 
//      model.set(GRB_DoubleParam_TimeLimit, 60);
      
      n = varsi.size();
      vars.resize(n);vcos.resize(n);dual = new bool[n];
      for(int i=0;i<n;i++){
        var v = varsi[i];
        dual[i] = v.radian;
        if(v.radian){
          if(v.lb!=v.ub){
            vars[i] = model.addVar(-1,1,0,GRB_CONTINUOUS);
            vcos[i] = model.addVar(-1,1,0,GRB_CONTINUOUS);
            model.addQConstr(vars[i]*vars[i]+vcos[i]*vcos[i],GRB_GREATER_EQUAL,1); //overestimate
          }else{
            vars[i] = model.addVar(sin(v.lb),sin(v.lb),0,GRB_CONTINUOUS);
            vcos[i] = model.addVar(cos(v.lb),cos(v.lb),0,GRB_CONTINUOUS);
          }
        }else
          vars[i] = model.addVar(varsi[i].lb<-1e19?-GRB_INFINITY:varsi[i].lb,
                                 varsi[i].ub> 1e19? GRB_INFINITY:varsi[i].ub,0,GRB_CONTINUOUS);
      }

      model.setObjective(parseExpr(fi->getF()), GRB_MINIMIZE);

      for(int i=0;i<phix.size();i++)
        model.addQConstr(parseExpr(phix[i]->getF()),GRB_GREATER_EQUAL,0);
    }

    void optimize(){
      res = new double[n];
      model.write("test.lp");
      model.optimize();
std::cout << "{";
      for(int i=0;i<n;i++){
        if(dual[i]){
          double s = vars[i].get(GRB_DoubleAttr_X),c = vcos[i].get(GRB_DoubleAttr_X);
          res[i] = atan2(s,c);
        }else
          res[i] = vars[i].get(GRB_DoubleAttr_X);
std::cout << res[i] << ",";
      }
std::cout << "}" << std::endl;
    }
};
