using namespace Ipopt;

class Objective{
    public:
        virtual Number eval(Index n, const Number* x) = 0;
        virtual void grad(Index n, const Number* x,  Number* grad) = 0;
    };

class FirstVar: public Objective{
    public:
        FirstVar(){};

        Number eval(Index n, const Number* x){
            return x[0];
        }

        void grad(Index n, const Number* x,  Number* grad){
            grad[0] = 1;
            for(int i=1;i<n;i++)
                grad[i] = 0;
        }
};
