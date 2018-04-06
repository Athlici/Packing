struct var{
    double lb,ub;
    string name;
    var(double lbi,double ubi) : lb(lbi),ub(ubi) {};
};

struct point{
    point() {};
    point(double xi,double yi) : x(xi),y(yi) {};
    double x,y;
};

struct circle{
    circle(point pi,double ri) : p(pi),r(ri) {};
    point p;
    double r;
};
