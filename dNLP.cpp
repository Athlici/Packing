#include "IpTNLP.hpp"

using namespace Ipopt;

class dNLP : public TNLP{
  private:
    Objective* f;
    vector<var> vars;
    vector<PhiFuncPrim*> phix;
    const double* x0;
    vector<int> jacl;
    vector<vector<int>> hesmap;
    vector<tuple<int,int>> hesind;
    
  public:
    dNLP(Objective* fi,vector<var> varsi,vector<PhiFuncPrim*> phixi,const double* x0i) : f(fi),vars(varsi),phix(phixi),x0(x0i) {};

    bool get_nlp_info(Index& n, Index& m, Index& njac, Index& nhes, IndexStyleEnum& index_style){
      n = vars.size();
      m = phix.size();

      jacl.resize(m);
      njac = 0;
      for(int i=0;i<m;i++){
        jacl[i] = phix[i]->getD1ind().size();
        njac += jacl[i];
      }
  
      int hesnz[n*n] = {0};
      for(int i=0;i<m;i++){
        vector<tuple<int,int>> tmp = phix[i]->getD2ind();
        for(int j=0;j<tmp.size();j++)
          hesnz[n*get<0>(tmp[j])+get<1>(tmp[j])] = 1;
      }
      nhes = 0;
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
          if(hesnz[i*n+j]>0){
            hesnz[i*n+j]=++nhes;
            hesind.push_back({i,j});
          }
      hesmap.resize(m);
      for(int i=0;i<m;i++){
        vector<tuple<int,int>> tmp = phix[i]->getD2ind();
        hesmap[i].resize(tmp.size());
        for(int j=0;j<tmp.size();j++)
          hesmap[i][j] = hesnz[n*get<0>(tmp[j])+get<1>(tmp[j])]-1;
      }
  
      index_style = TNLP::C_STYLE;
      return true;
    }
  
    bool get_bounds_info(Index n, Number* xl, Number* xu, Index m, Number* gl, Number* gu){
      for(int i=0;i<n;i++){
        xl[i] = vars[i].lb;
        xu[i] = vars[i].ub;
      }
      for(int i=0;i<m;i++){
        gl[i] = 0;
        gu[i] = 2e19;
      }
      return true;
    }
  
    bool get_starting_point(Index n, bool initx, Number* x, bool initz, 
            Number* zL, Number* zU, Index m, bool init_lambda, Number* lambda){
      for(int i=0;i<n;i++)
        x[i] = x0[i];
      return true;
    }
  
    bool eval_f(Index n, const Number* x, bool newx, Number& objv){
      objv = f->eval(n,x);
      return true;
    }
  
    bool eval_grad_f(Index n, const Number* x, bool newx, Number* gradf){
      f->grad(n,x,gradf);
      return true;
    }
  
    bool eval_g(Index n, const Number* x, bool newx, Index m, Number* g){
      for(int i=0;i<m;i++)
        g[i] = phix[i]->eval(x);
      return true;
    }
  
    bool eval_jac_g(Index n, const Number* x, bool newx, Index m, Index njac, Index* iRow, Index *jCol, Number* values){
      int ind = 0;
      if(values == NULL){
        for(int i=0;i<m;i++){
          vector<int> tmp = phix[i]->getD1ind();
          for(int j=0;j<tmp.size();j++){
            iRow[ind] = i;
            jCol[ind] = tmp[j];
            ind++;
          }
        }
      } else {
        for(int i=0;i<m;i++){
          phix[i]->getD1(x,&values[ind]);
          ind += jacl[i];
        }
      }
      return true;
    }
  
    bool eval_h(Index n, const Number* x, bool newx, Number sigmaf, Index m, const Number* lambda,
                      bool newlambda, Index nhes, Index* iRow, Index* jCol, Number* values){
      //TODO add objective hessian
      if(values == NULL){
        for(int i=0;i<nhes;i++){
          iRow[i] = get<0>(hesind[i]);
          jCol[i] = get<1>(hesind[i]);
        }
      } else {
        for(int i=0;i<nhes;i++)
          values[i] = 0;
        for(int i=0;i<m;i++){
          int s = hesmap[i].size();
          double val[s];
          phix[i]->getD2(x,val);
          for(int j=0;j<s;j++)
            values[hesmap[i][j]] += lambda[i]*val[j];
        }
      }
      return true;
    }

  void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, 
        const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq){};

private:
  dNLP(const dNLP&);
  dNLP& operator=(const dNLP&);
};
