struct var{
    var() {};
    var(double lbi,double ubi) : lb(lbi),ub(ubi) {};
    double lb,ub;
    string name;
};

struct point{  //TODO: make Vector2d
    point() {};
    point(double xi,double yi) : x(xi),y(yi) {};
    double x,y;
};

struct circle{
    circle() {};
    circle(point pi,double ri) : p(pi),r(ri) {};
    point p;
    double r;
};
