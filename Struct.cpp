struct var{  //variables
    var() {};
    var(double lbi,double ubi) : lb(lbi),ub(ubi) {};
    double lb,ub;
    string name;
};

class point{  //points TODO: make Vector2d?
    public:
        point() {};
        point(double xi,double yi) : x(xi),y(yi) {};
        double x,y;

        void move(point p){
            x+=p.x; y+=p.y;
        }
};

class circle{  //circles
    public:
        circle() {};
        circle(point pi,double ri) : p(pi),r(ri) {};
        point p;
        double r;

        void move(point pd){
            p.move(pd);
        }
};
