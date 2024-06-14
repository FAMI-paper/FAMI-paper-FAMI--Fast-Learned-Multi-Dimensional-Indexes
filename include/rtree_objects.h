#ifndef __OBJECTS
#define __OBJECTS

#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <iostream>

template<typename T>
class Point{
    private:
    //2D - Coordinates of Point.  
    T X;
    T Y;

    public:
    //Define a constructor.
    Point(T x = 0, T y = 0): X(x), Y(y){};
    //Overload of operators:
    //A point is equal or less than other if its 'X' value and it 'Y' value
    //are equal or less than the 'X' and 'Y' value of the other Point.
    bool operator>=(const Point& p){
          if(this->X >= p.X && this->Y >=p.Y)
            return true;
          else
            return false;
    }
    //A point is equal or greater than other if its 'X' value and it 'Y' value
    //are equal or greater than the 'X' and 'Y' value of the other Point.
    bool operator<=(const Point& p){
          if(this->X <= p.X && this->Y <= p.Y)
            return true;
          else
            return false;
    }

    bool operator==(const Point& p){
          if(this->X == p.X && this->Y == p.Y)
            return true;
          else
            return false;
    }
    //setter & getter functions:
    T get_X(){return this->X;}
    T get_Y(){return this->Y;}
    void set_X(T x){this->X = x;}
    void set_Y(T y){this->Y = y;}
    void set_X_Y(T x, T y){
        this->X = x;
        this->Y = y;
    }
};

template<typename T>
class Polygon{
private: 
    //MBB of the Polygon is represented 
    //by its minimum Point & maximun Point
    Point<T> Pmin;
    Point<T> Pmax;
    int key;
    
    //Colection of Points that represent all vertices
    std::vector<Point<T>> vertices;
    void area_added(Polygon<T> &);
    void area_added(Polygon<T> &, Point<T> & , Point<T> & );
    
public:
    int corners; // number of vertices
    //Initialize the Vertices of Polygon
    Polygon(std::vector<Point<T>> ) ;
    //It represent a rectangle region.
    Polygon(Point<T> pmin, Point<T> pmax):Pmin(pmin),Pmax(pmax), corners(4), key(-1){};
    //Generic Polygon, a point is Polygon
    Polygon(Point<T>);

    bool operator==(const Polygon<T>& polypon){
        if(Pmax == polypon.Pmax && Pmin == polypon.Pmin){
          for(int i = 0; i<polypon.vertices.size(); i++){
            Point<T> point = polypon.vertices[i];
            auto iter = find(vertices.begin(),vertices.end(),point);
            if(iter == vertices.end()){
              return false;
            }
          }
          return true;
        }else{
          return false;
        }
    }

    //Cost of add a polygon.
    T cost_two_polygons(Polygon<T> & );
    Point<T> get_Pmax(){return Pmax;}
    Point<T> get_Pmin(){return Pmin;}
    void set_Polygon(Point<T> min, Point<T> max){this->Pmin = min; this->Pmax = max;}
    bool intersect_with_BB(Polygon<T> &);
    bool traverse_with(Polygon<T> &);
    bool is_Within_of(const Polygon<T> &);
    Polygon<T> get_mbb();
    //To get KNN query- Geometric distance.
    //template <class >
    T distance_geometric(Point<T> );
    //
    //template <class T>
    T max_distance_geometric(Point<T> );

    std::vector<Point<T>> get_vertices();
    bool set_key(int );
    int get_key(){return this->key;}

    ~Polygon();
};


template<typename T>
Polygon<T>::Polygon(std::vector<Point<T>> p):vertices(p), key(-1),corners(p.size()){
    Polygon<T> mbb = this->get_mbb();
    this->Pmin = mbb.Pmin;
    this->Pmax = mbb.Pmax;   
}

template<typename T>
Polygon<T>::Polygon(Point<T> P):corners(1), key(-1) {
    this->Pmin = P;
    this->Pmax = P;
    this->vertices.push_back(P);
}

template<typename T>
void Polygon<T>::area_added(Polygon<T> & reg){
    T x_max = this->Pmax.get_X();
    if(reg.Pmax.get_X() > x_max){
        x_max = reg.Pmax.get_X();
    }
    T x_min = this->Pmin.get_X();
    if(reg.Pmin.get_X() < x_min){
        x_min = reg.Pmin.get_X();
    }

    T y_max = this->Pmax.get_Y();
    if(reg.Pmax.get_Y() > y_max){
        y_max = reg.Pmax.get_Y();
    }
    T y_min = this->Pmin.get_Y();
    if(reg.Pmin.get_Y() < y_min){
        y_min = reg.Pmin.get_Y();
    }
        //This points max and minimung will welp to calculate the region of parents
    this->Pmin = Point<T>(x_min,y_min);
    this->Pmax = Point<T>(x_max,y_max);
}

template<typename T>
void Polygon<T>::area_added(Polygon<T> & reg, Point<T> & pmin, Point<T> & pmax){
    T x_max = this->Pmax.get_X();
    if(reg.Pmax.get_X() > x_max){
        x_max = reg.Pmax.get_X();
    }
    T x_min = this->Pmin.get_X();
    if(reg.Pmin.get_X() < x_min){
        x_min = reg.Pmin.get_X();
    }

    T y_max = this->Pmax.get_Y();
    if(reg.Pmax.get_Y() > y_max){
        y_max = reg.Pmax.get_Y();
    }
    T y_min = this->Pmin.get_Y();
    if(reg.Pmin.get_Y() < y_min){
        y_min = reg.Pmin.get_Y();
    }
        //This points max and minimung will welp to calculate the region of parents
    pmin = Point<T>(x_min,y_min);
    pmax = Point<T>(x_max,y_max);
}

template<typename T>
T Polygon<T>::cost_two_polygons(Polygon<T> & reg){
    Point<T> pmin, pmax;
    area_added(reg,pmin,pmax);
    T d = (pmax.get_X()-pmin.get_X())*(pmax.get_Y()-pmin.get_Y());
    d -= (reg.Pmax.get_X() - reg.Pmin.get_X())*(reg.Pmax.get_Y() - reg.Pmin.get_Y());
    d -= (this->Pmax.get_X() - this->Pmin.get_X())*(this->Pmax.get_Y() - this->Pmin.get_Y());
    return d;
}

template<typename T>
bool Polygon<T>::intersect_with_BB(Polygon<T> & pol){
    if(pol.traverse_with(*this))
        return true;
    if(pol.is_Within_of(*this))
        return true;
    if(this->Pmin <= pol.Pmax && this->Pmin >= pol.Pmin){
        return true;
    }
    if(this->Pmax <= pol.Pmax && this-> Pmax >= pol.Pmin){
        return true;
    }
    Point<T> myPoint_leftUP(this->Pmin.get_X(),this->Pmax.get_Y());
    Point<T> myPoint_right_L(this->Pmax.get_X(), this->Pmin.get_Y());
    Point<T> BB_point_leftUP(pol.Pmin.get_X(),pol.Pmax.get_Y());
    Point<T> BB_point_right_L(pol.Pmax.get_X(), pol.Pmin.get_Y());
    
    if(myPoint_leftUP <= BB_point_leftUP && myPoint_leftUP >= BB_point_right_L){
        return true;
    }
    if(myPoint_right_L <= BB_point_leftUP && myPoint_right_L >= BB_point_right_L){
        return true;
    }
    
    //if(this->traberse_with(pol))
    //    return true;
    
    return false;
}

template<typename T>
bool Polygon<T>::is_Within_of(const Polygon<T> & pol){
    if(this->Pmin >= pol.Pmin && this->Pmax <= pol.Pmax)
            return true;
    else
        return false;
}

template<typename T>
bool Polygon<T>::traverse_with(Polygon<T> & other){
    if(this->Pmin.get_X() <= other.Pmax.get_X() && this->Pmax.get_Y() <= other.Pmax.get_Y())
        return true;
    if(this->Pmin.get_Y() <= other.Pmax.get_Y() && this->Pmax.get_X() <= other.Pmax.get_X())
        return true;

    return false;
}

// build a rectangle that includes all points
template<typename T>
Polygon<T> Polygon<T>::get_mbb(){      
    T x_min= std::numeric_limits<T>::max(); 
    T y_min = x_min;
    T x_max = std::numeric_limits<T>::min();
    T y_max = x_max;
    for(int i = 0; i < this->corners; i++){
        if(this->vertices[i].get_X() < x_min)
            x_min = this->vertices[i].get_X();
        if(this->vertices[i].get_Y() < y_min)
            y_min = this->vertices[i].get_Y();
        if(this->vertices[i].get_X() > x_max)
            x_max = this->vertices[i].get_X();
        if(this->vertices[i].get_Y() > y_max)
            y_max = this->vertices[i].get_Y();
    }
    return Polygon<T>(Point<T>(x_min,y_min),Point<T>(x_max,y_max));
}

//template <class T>
//T Polygon::distance_geometric(Point q){return 0.0;}
template<typename T>
T Polygon<T>::distance_geometric(Point<T> q){
    T d_X_min = std::numeric_limits<T>::max();
    if(q.get_X() >= this->Pmin.get_X() && q.get_X() <= this->Pmax.get_X()){
        d_X_min = 0.0;
    }
    else{
        if(abs(q.get_X() - this->Pmin.get_X()) < d_X_min){
            d_X_min = abs(q.get_X() - this->Pmin.get_X());
        }
        if(abs(q.get_X() - this->Pmax.get_X()) < d_X_min){
            d_X_min = abs(q.get_X() - this->Pmax.get_X());
        }
    }
    T d_Y_min = std::numeric_limits<T>::max();
    if(q.get_Y() >= this->Pmin.get_Y() && q.get_Y() <= this->Pmax.get_Y()){
        d_Y_min = 0.0;
    }
    else{
        if(abs(q.get_Y() - this->Pmin.get_Y()) < d_Y_min){
            d_Y_min = abs(q.get_Y() - this->Pmin.get_Y());
        }
        if(abs(q.get_Y() - this->Pmax.get_Y()) < d_Y_min){
            d_Y_min = abs(q.get_Y() - this->Pmax.get_Y());
        }
    }
    T d = sqrt(d_X_min*d_X_min + d_Y_min*d_Y_min);
    return d;
}
template<typename T>
T Polygon<T>::max_distance_geometric(Point<T> q){
    T d_X_max = std::numeric_limits<T>::min();;
    if(abs(q.get_X() - this->get_Pmin().get_X()) > d_X_max){
        d_X_max = abs(q.get_X() - this->get_Pmin().get_X());
    }
    if(abs(q.get_X() - this->get_Pmax().get_X()) > d_X_max){
        d_X_max = abs(q.get_X() - this->get_Pmax().get_X());
    }
    T d_Y_max = std::numeric_limits<T>::min();;
    if(abs(q.get_Y() - this->get_Pmin().get_Y()) > d_Y_max){
        d_Y_max = abs(q.get_Y() - this->get_Pmin().get_Y());
    }
    if(abs(q.get_Y() - this->get_Pmax().get_Y()) > d_Y_max){
        d_Y_max = abs(q.get_Y() - this->get_Pmax().get_Y());
    }
    T d = sqrt(d_X_max*d_X_max+d_Y_max*d_Y_max);
    return d;
}

template<typename T>
std::vector<Point<T>> Polygon<T>::get_vertices(){
    return this->vertices;
}

template<typename T>
bool Polygon<T>::set_key(int k){
    if(this->key==-1){
        this->key = k;
        return true;
    }
    return false;
}

template<typename T>
Polygon<T>::~Polygon(){    
}
#endif