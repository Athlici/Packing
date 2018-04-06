PhiPolygon regularPolygon(int n){
    double phi=2*M_PI/n;
    vector<point> v(n);
    for(int i=0;i<n;i++)
        v[i]=point(sin(phi*i),cos(phi*i));
    return PhiPolygon(v);
}
