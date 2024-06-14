#ifndef __RTREE
#define __RTREE

#include "rtree_objects.h"
#include <string>
#include <iostream>

template<typename T>
struct d_leaf{
    public: 
    Polygon<T> * polygon;
    Polygon<T> * region;
    d_leaf(Polygon<T> *p = nullptr, Polygon<T> * r = nullptr):polygon(p),region(r){};
    void set_data(Polygon<T> * p , Polygon<T> * r ){polygon = p; region = r;}
    bool operator==(const d_leaf<T>& leaf){
        if(polygon == leaf.polygon && region == leaf.region){
            return true;
        }else{
            return false;
        }
    }
};

template<typename T>
struct d_internal_node;

template<typename T>
class RTree;

template<typename T>
class RTree_node{
private:
    int M; // max num of node
    int elements; // current num of node
    bool is_leaf;
    int at_level;
    RTree_node<T> * father;
    std::vector<d_leaf<T>> data_leafs;
    std::vector<d_internal_node<T>> data_internal_node;

    //choose the first nodes in quadratic split
    void choose_origin(int & , int & );
    //Get the MBB of the node.
    Polygon<T> mbb_node();
public:
    int get_level(){return at_level;}
    RTree_node(RTree_node<T> * f, bool _l, int _M, int _lvl);
    ~RTree_node();
    friend class RTree<T>;
};

template<typename T>
struct d_internal_node{
    public:
    Polygon<T> * region;
    RTree_node<T> * child;
    d_internal_node(Polygon<T> *r = nullptr, RTree_node<T> * c = nullptr):region(r),child(c){};
    void set_data(Polygon<T> * r, RTree_node<T> * c){region = r; child = c;}
};

template<typename T>
struct data_query_return{
    public:
    Polygon<T> * Pol;
    int lvl;
    data_query_return(Polygon<T> *p, int _lvl):Pol(p),lvl(_lvl){};
};

template<typename T>
class RTree{
private:
    RTree_node<T> * root;
    int M;
    int m; 
    int H; //Height of tree
    int indx;
    int items_cnt=0;
    //Insert Polygon internaly
    RTree_node<T> * insert_polygon(RTree_node<T> *, d_leaf<T>);
    //general function to insert
    void insert(RTree_node<T> *&, d_leaf<T>);
    //Insert element in an internal node.
    void insert_internal_region(RTree_node<T> *, d_internal_node<T>  );
    //Distribute polygons on leaf:
    void distribute_polygons(RTree_node<T> * , RTree_node<T> * , const std::vector<d_leaf<T>> &);
    //Distribute regions on internal nodes:
    void distribute_regions(RTree_node<T> *, RTree_node<T> *, const std::vector<d_internal_node<T>> &);
    //Quadratic split on leafs:
    RTree_node<T> * quadratic_split(RTree_node<T> *);
    //Quadratic split on internal nodes:
    RTree_node<T> *  quadratic_split_internal_nodes(RTree_node<T> *);
    //Select the leaf where the new Polygon try to be inserted
    RTree_node<T> * select_leaf(RTree_node<T> * , Polygon<T> * );
    //Adjust the tree after of insert a polygon
    void adjust_tree(RTree_node<T> *, RTree_node<T> * );
    //Recursive Range Search.
    void range_search_recursive(RTree_node<T> * , Polygon<T> & , std::vector<Polygon<T> *> &);
    // Find leaf that includes delete entry
    RTree_node<T> * find_leaf(RTree_node<T> * , d_leaf<T> );
    // condense tree
    void condense_tree(RTree_node<T> * ,std::vector<d_leaf<T>> &,std::vector<d_internal_node<T>> &,std::vector<int> &);
    // find internal node when condensing tree
    void find_internal_node(RTree_node<T> * ,std::vector<RTree_node<T>*> &, int );

    //Methods to get KNN

    void DFS_recursive(Point<T> , int , RTree_node<T> * , std::vector<d_leaf<T> *> & , std::vector<T> & ,T & );
    //insert sort
    
    void insert_sort(std::vector<T> & , std::vector<d_leaf<T>*> & );

    void insert_sort(std::vector<T> & , std::vector<RTree_node<T>*> &);

public:
    RTree(){};
    RTree(int _M): M(_M), m((_M+1)/2),root(nullptr), H(0), indx(0){};
    
    //Insert Polygon in Front-end
    void insert_polygon(Polygon<T> * , Polygon<T> *);

    void delete_polygon(Polygon<T> * , Polygon<T> *);

    //Range search in Front
    void range_search(Polygon<T> , std::vector<Polygon<T> *> &);
    
    
    //get the k-nearest neighbor Polygons.
    void k_NN_DFS(Point<T> , int , std::vector<d_leaf<T>*> &);

    void getAll_values(RTree_node<T> *,std::vector<Polygon<T> *> &);

    std::string show_values_JSON();

    void showAll_values_JSON(RTree_node<T> *, int , std::string &);

    void get_polygons_JSON(const std::vector<d_leaf<T>*> &, std::string &);

    std::string get_regions_JSON(RTree_node<T>* );

    void get_Range_Search_JSON(const std::vector<Polygon<T> *> &, std::string &);
    
    int get_items_cnt();

    // ~RTree();

};


template<typename T>
RTree_node<T>::RTree_node(RTree_node<T> * f, bool _l, int _M, int _lvl):is_leaf(_l), M(_M), elements(0), father(f),at_level(_lvl){
    if(this->is_leaf){
        this->data_leafs.resize(M+1); 
    }
    else{
        this->data_internal_node.resize(M+1);
    }
}

/*
    Use quadratic method to pick two rectangles as seeds when split occurs
*/
template<typename T>
void RTree_node<T>::choose_origin(int & a, int & b){
    int i = 0;
    T d = std::numeric_limits<T>::min();
    
    while(i < elements -1){
        for(int j = i+1; j < elements; j++){
            if(this->is_leaf){
                T costo = this->data_leafs[i].region->cost_two_polygons(*this->data_leafs[j].region);
                if( costo > d){
                    a = i;
                    b = j;
                    d = costo;
                }
            }
            else{
                T costo = this->data_internal_node[i].region->cost_two_polygons(*this->data_internal_node[j].region);
                if( costo > d){
                    a = i;
                    b = j;
                    d = costo;
                }
            }
        }
        i++;
    }
}

template<typename T>
Polygon<T> RTree_node<T>::mbb_node(){
    T x_min = std::numeric_limits<T>::max(); 
    T x_max = std::numeric_limits<T>::min();
    T y_min = x_min, y_max = x_max;
    if(this->is_leaf){
        for( int i = 0; i < this->elements; i++){
            if(this->data_leafs[i].region->get_Pmin().get_X() < x_min){
                x_min = this->data_leafs[i].region->get_Pmin().get_X();
            }
            if(this->data_leafs[i].region->get_Pmin().get_Y() < y_min){
                y_min = this->data_leafs[i].region->get_Pmin().get_Y();
            }
                
            if(this->data_leafs[i].region->get_Pmax().get_X() > x_max){
                x_max = this->data_leafs[i].region->get_Pmax().get_X();
            }
            if(this->data_leafs[i].region->get_Pmax().get_Y() > y_max){
                y_max = this->data_leafs[i].region->get_Pmax().get_Y();
            }
        }
    }
    else{
        for( int i = 0; i < this->elements; i++){
            if(this->data_internal_node[i].region->get_Pmin().get_X() < x_min){
                x_min = this->data_internal_node[i].region->get_Pmin().get_X();
            }
            if(this->data_internal_node[i].region->get_Pmin().get_Y() < y_min){
                y_min = this->data_internal_node[i].region->get_Pmin().get_Y();
            }
                
            if(this->data_internal_node[i].region->get_Pmax().get_X() > x_max){
                x_max = this->data_internal_node[i].region->get_Pmax().get_X();
            }
            if(this->data_internal_node[i].region->get_Pmax().get_Y() > y_max){
                y_max = this->data_internal_node[i].region->get_Pmax().get_Y();
            }
        }
    }
    
    return Polygon<T>(Point<T>(x_min,y_min),Point<T>(x_max,y_max));    
}

template<typename T>
RTree_node<T> * RTree<T>::insert_polygon(RTree_node<T> * node, d_leaf<T> data){
    //If there isn't any node, create one
    RTree_node<T> * posible_Brother = nullptr;
        //There's space in the leaf?
    if(node->elements < node->M){
        if(data.polygon->set_key(this->indx))
            this->indx++;
        node->data_leafs[node->elements] = data;
        node->elements++;            
    }
    else{
        node->data_leafs[node->elements] = data;
        node->elements++;
        posible_Brother = quadratic_split(node);
    }
    return posible_Brother;
}
/*
bool RTree::insert(RTree_node *& node, d_leaf data){
    RTree_node * posible_Brother = nullptr;
    RTree_node * leaf_chose = nullptr;
    if(this->root == nullptr){
        this->root = new RTree_node(true, this->M, 0);
    }
    leaf_chose = select_leaf(this->root, data.region);
    posible_Brother = insert_polygon(leaf_chose, data);
    adjust_tree(leaf_chose,posible_Brother);
    if(this->root->elements > M ){
        RTree_node * internal_brother = nullptr;
        RTree_node * node = this->root; 
        internal_brother = quadratic_split_internal_nodes(this->root);
        adjust_tree(node, internal_brother);
        this->H++;
    }
}*/

template<typename T>
void RTree<T>::insert(RTree_node<T> * &node, d_leaf<T> data){
    RTree_node<T> * posible_Brother = nullptr;
    RTree_node<T> * leaf_chose = nullptr;
    if(this->root == nullptr){
        this->root = new RTree_node<T>(nullptr, true, this->M, 0);
    }
    leaf_chose = select_leaf(this->root, data.region);
    /*
    std::string json = "";
    showAll_values_JSON(leaf_chose,0,json);
    std::cout << json << std::endl;
    */
    posible_Brother = insert_polygon(leaf_chose, data);
    adjust_tree(leaf_chose,posible_Brother);
    if(this->root->elements > M ){
        RTree_node<T> * internal_brother = nullptr;
        RTree_node<T> * node = this->root; 
        internal_brother = quadratic_split_internal_nodes(this->root);
        adjust_tree(node, internal_brother);
        this->H++;
    }
    
}

template<typename T>
RTree_node<T> * RTree<T>::quadratic_split(RTree_node<T> * node){
    int i, j;
    node->choose_origin(i,j);
    if(node->father == nullptr){
        node->father = new RTree_node<T>(nullptr,false,this->M, H+1);
        //Just insert elements.
        insert_internal_region(node->father, d_internal_node<T>(node->data_leafs[i].region, node));
        this->root = node->father;
        this->H++;
    }
    else{
        //father region must be equal to the node region 'i'
        for(int m = 0; m < node->father->elements; m++){
            if(node->father->data_internal_node[m].child == node){
                *node->father->data_internal_node[m].region = *node->data_leafs[i].region;
            }
        }
    }
    RTree_node<T> * brother = new RTree_node<T>(node->father, true, this->M, node->at_level);
    insert_internal_region(brother->father, d_internal_node<T>(node->data_leafs[j].region,brother));
    insert_polygon(brother, node->data_leafs[j]);
    std::vector<d_leaf<T>> tmp;
    for(int n = 0; n < node->elements; n++){
        if(n != i && n != j){
            tmp.push_back(node->data_leafs[n]);
        }
    }
    //Clear the node, just with its 'i' index data
    node->data_leafs[0] = node->data_leafs[i];
    for(int m = 1; m < node->elements; m++){
        node->data_leafs[m].set_data(nullptr,nullptr);
    }
    node->elements = 1;
    distribute_polygons(node, brother,tmp);
    return brother;
}

template<typename T>
RTree_node<T> *  RTree<T>::quadratic_split_internal_nodes(RTree_node<T> * node){
    int i, j;
    node->choose_origin(i,j);
    if(node->father == nullptr){
        node->father = new RTree_node<T>(nullptr,false,this->M, H+1);
        insert_internal_region(node->father, d_internal_node<T>(node->data_internal_node[i].region, node));
        this->root = node->father;
    }
    else{
        //father region must be equal to the node region 'i'
        for(int m = 0; m < node->father->elements; m++){
            if(node->father->data_internal_node[m].child == node){
                *node->father->data_internal_node[m].region = *node->data_internal_node[i].region;
            }
        }
    }
    
    RTree_node<T> * brother = new RTree_node<T>(node->father, false, this->M, node->at_level);
    insert_internal_region(brother->father, d_internal_node<T>(node->data_internal_node[j].region,brother));
    node->data_internal_node[j].child->father = brother;
    insert_internal_region(brother, node->data_internal_node[j]);

    std::vector<d_internal_node<T>> tmp;
    for(int n = 0; n < node->elements; n++){
        if(n != i && n != j){
            tmp.push_back(node->data_internal_node[n]);
        }
    }
    
    node->data_internal_node[0] = node->data_internal_node[i];
    for(int m = 1; m < node->elements; m++){
        node->data_internal_node[m].set_data(nullptr,nullptr);
    }
    node->elements = 1;

    distribute_regions(node, brother, tmp);
    return brother;
}

template<typename T>
void RTree<T>::insert_internal_region(RTree_node<T> * node, d_internal_node<T> data ){
    Polygon<T> * reg = new Polygon<T>(data.region->get_Pmin(),data.region->get_Pmax());
    if(reg->set_key(this->indx))
        this->indx++;
    data.region = reg;    
    node->data_internal_node[node->elements] = data;
    node->elements++;
}

template<typename T>
void RTree<T>::distribute_polygons(RTree_node<T> * node, RTree_node<T> * brother, const std::vector<d_leaf<T>> & data){
    int i, j ;
    //Search where is the node & brother in the father.
    for(int a = 0; a < node->father->elements; a++){
        if(node->father->data_internal_node[a].child == node){
            i = a;
        }
        else if(node->father->data_internal_node[a].child == brother){
            j = a;
        }
    }
    for(int m = 0; m < data.size(); m++){
        int cost_node = node->father->data_internal_node[i].region->cost_two_polygons(*data[m].region);
        int cost_brother = brother->father->data_internal_node[j].region->cost_two_polygons(*data[m].region);
        if(cost_node < cost_brother && node->elements < this->M - this->m + 1){
            insert_polygon(node, data[m]);
        }
        else if(cost_node > cost_brother && brother->elements < this->M - this->m + 1){
            insert_polygon(brother, data[m]);
        }
        else{
            if(node->elements < this->M - this->m + 1){
                insert_polygon(node, data[m]);
            }
            else{
                insert_polygon(brother, data[m]);
            }
        }
    }
}

template<typename T>
void RTree<T>::distribute_regions(RTree_node<T> * node, RTree_node<T> *brother, const std::vector<d_internal_node<T>> &data){
    int i, j ;
    //Search where is the node & brother in the father.
    for(int a = 0; a < node->father->elements; a++){
        if(node->father->data_internal_node[a].child == node){
            i = a;
        }
        else if(node->father->data_internal_node[a].child == brother){
            j = a;
        }
    }
    for(int m = 0; m < data.size(); m++){
        int cost_node = node->father->data_internal_node[i].region->cost_two_polygons(*data[m].region);
        int cost_brother = brother->father->data_internal_node[j].region->cost_two_polygons(*data[m].region);
        if(cost_node < cost_brother && node->elements < this->M - this->m + 1){
            insert_internal_region(node, data[m]);
            
        }
        else if(cost_node > cost_brother && brother->elements < this->M - this->m + 1){
            data[m].child->father = brother;
            insert_internal_region(brother, data[m]);
        }
        else{
            if(node->elements < this->M - this->m + 1){
                insert_internal_region(node, data[m]);
                
            }
            else{
                data[m].child->father = brother;
                insert_internal_region(brother, data[m]);
            
            }
        }
    }
}

template<typename T>
void RTree<T>::adjust_tree(RTree_node<T> * node, RTree_node<T> *brother){
 
    if(node->father != nullptr){
        for(int m = 0; m < node->father->elements; m++){
            if(node->father->data_internal_node[m].child == node){
                int key = node->father->data_internal_node[m].region->get_key();
                *node->father->data_internal_node[m].region = node->mbb_node();
                node->father->data_internal_node[m].region->set_key(key);
            }
            else if(node->father->data_internal_node[m].child == brother){
                int key = brother->father->data_internal_node[m].region->get_key();
                *brother->father->data_internal_node[m].region = brother->mbb_node();
                brother->father->data_internal_node[m].region->set_key(key);
            }
        }
        if(node->elements > this->M && !node->is_leaf){
                
            RTree_node<T> * internal_brother = nullptr;
            internal_brother = quadratic_split_internal_nodes(node);
            adjust_tree(node, internal_brother);
        }
        else{
            adjust_tree(node->father,nullptr);   
        }
            
    }
        
    //}
}

template<typename T>
void RTree<T>::insert_polygon(Polygon<T> * pol, Polygon<T> * reg){
    this->items_cnt++;
    RTree<T>::insert(this->root, d_leaf<T>(pol,reg));
}

template<typename T>
RTree_node<T> * RTree<T>::select_leaf(RTree_node<T> * node, Polygon<T> * p_region){
    if(node->is_leaf){
        return node;
    }
    else{
        int index_min;
        T min_cost = std::numeric_limits<T>::max();
        for(int i = 0; i < node->elements; i++){
            T child_cost = node->data_internal_node[i].region->cost_two_polygons(*p_region);
            if(child_cost < min_cost){
                min_cost = child_cost;
                index_min = i;
            }
        }
        return select_leaf(node->data_internal_node[index_min].child,p_region);
    } 
}

template<typename T>
void RTree<T>::delete_polygon(Polygon<T> * pol, Polygon<T> * reg){
    RTree_node<T>* tmp = find_leaf(this->root,d_leaf<T>(pol,reg));

    if(tmp != nullptr){
        for(int i = 0; i < tmp->elements; i++){
            if(tmp->data_leafs[i] == d_leaf<T>(pol,reg)){
                tmp->data_leafs.erase(tmp->data_leafs.begin()+i);
                tmp->elements--;
                break;
            }
        }
    }

    std::vector<d_leaf<T>> d_leaf_list;
    std::vector<d_internal_node<T>> d_internal_node_list;
    std::vector<int> d_internal_node_level_list;
    condense_tree(tmp,d_leaf_list,d_internal_node_list,d_internal_node_level_list);

}

template<typename T>
RTree_node<T> * RTree<T>::find_leaf(RTree_node<T> * node, d_leaf<T> data){
    if(node->is_leaf){
        for(int i = 0; i < node->elements; i++){
            if(node->data_leafs[i] == data){
                return node;
            }
        }
        
        return nullptr;
    }else{
        for(int i = 0; i< node->elements; i++){
            if(node->data_internal_node[i].region->intersect_with_BB(*data.region)){
                RTree_node<T>* tmp = find_leaf(node->data_internal_node[i].child, data);
                if(tmp != nullptr){
                    return tmp;
                } 
            }
        }
        return nullptr;
    }
}

template<typename T>
void RTree<T>::condense_tree(RTree_node<T> * node,std::vector<d_leaf<T>> &d_leaf_list,std::vector<d_internal_node<T>> &d_internal_node_list,std::vector<int> &d_internal_node_level_list){

    std::vector<RTree_node<T>*> res;

    if(node == this->root){
        if(!d_leaf_list.empty()){
            for(int i = 0; i < d_leaf_list.size(); i++){
                insert_polygon(d_leaf_list[i].polygon,d_leaf_list[i].region);
            }
        }
        
        if(!d_internal_node_list.empty()){
            for(int i = 0; i< d_internal_node_list.size(); i++){
                int level = d_internal_node_level_list[i];
                find_internal_node(this->root,res, level);
                if(!res.empty()){
                    insert_internal_region(res[0],d_internal_node_list[i]);
                }else{
                    std::cout << "Don't find internal node" << std::endl;
                }
                res.clear();
            }
        }
    }else{
        for(int i = 0; i < node->father->elements; i++){
            if(node->father->data_internal_node[i].child == node){
                if(node->elements < this->m){
                    if(node->is_leaf){
                        for(int j = 0; j < node->elements; j++){
                            d_leaf_list.push_back(node->data_leafs[j]);
                        }
                    }else{
                        for(int j = 0; j < node->elements; j++){
                            d_internal_node_list.push_back(node->data_internal_node[j]);
                            d_internal_node_level_list.push_back(node->at_level);
                        }
                    }
                    node->father->data_internal_node.erase(node->father->data_internal_node.begin()+i);
                    node->father->elements--;
                    condense_tree(node->father,d_leaf_list,d_internal_node_list,d_internal_node_level_list);
                }else{
                    *node->father->data_internal_node[i].region = node->mbb_node();
                    condense_tree(node->father,d_leaf_list,d_internal_node_list,d_internal_node_level_list);
                }
            }
        }

    }

}

template<typename T>
void RTree<T>::find_internal_node(RTree_node<T> * node,std::vector<RTree_node<T>*> &res, int level){
    if(node->at_level == level && node->elements < this->M){
        res.push_back(node);
        return;
    }else{
        for(int i = 0; i < node->elements; i++){
            if(!node->is_leaf){
                find_internal_node(node->data_internal_node[i].child,res, level);
            }
        }
    }
}

/*
    Range search.
    It returns all the path Regions & Polygons Pointers and its respective level.
    The level will help to draw the path query.
*/
template<typename T>
void RTree<T>::range_search_recursive(RTree_node<T> * node, Polygon<T> & query, std::vector<Polygon<T> *> & ans){
    if(node != nullptr){
        if(!node->is_leaf){
            for(int i = 0; i < node->elements;i++){
                if(node->data_internal_node[i].region->is_Within_of(query))
                {
                    getAll_values(node->data_internal_node[i].child,ans);
                    // continue;
                }
                else if(node->data_internal_node[i].region->intersect_with_BB(query) ){
                    //ans.push_back(data_query_return(node->data_internal_node[i].region,node->get_level()));
                    range_search_recursive(node->data_internal_node[i].child,query,ans);
                }
            }
        }
        else{
            for(int i = 0; i < node->elements; i++){
                if(node->data_leafs[i].region->is_Within_of(query)){
                    ans.push_back(node->data_leafs[i].polygon);
                }
            }
        }
    }
}
template<typename T>
void RTree<T>::range_search(Polygon<T> query, std::vector<Polygon<T> *> & ans){
    T x_min = query.get_Pmin().get_X();
    T y_min = query.get_Pmin().get_Y();
    if(query.get_Pmax().get_X() < x_min)
        x_min = query.get_Pmax().get_X();
    if(query.get_Pmax().get_Y() < y_min)
        y_min = query.get_Pmax().get_Y();
    T x_max = query.get_Pmin().get_X();
    T y_max = query.get_Pmin().get_Y();
    if(query.get_Pmax().get_X() > x_max)
        x_max = query.get_Pmax().get_X();
    if(query.get_Pmax().get_Y() > y_max)
        y_max = query.get_Pmax().get_Y();
    Polygon<T> q(Point<T>(x_min,y_min),Point<T>(x_max,y_max));
    range_search_recursive(this->root, q, ans);
}

template<typename T>
void RTree<T>::DFS_recursive(Point<T> q, int k, RTree_node<T> * node, std::vector<d_leaf<T> *> & L, std::vector<T> & ddk,T & poor){
	if(node->is_leaf){
        for(int i = 0; i < node->elements; i++){
            ddk.push_back(node->data_leafs[i].region->distance_geometric(q));
            L.push_back(&node->data_leafs[i]);
        }
        if(ddk.size()>k)
        {
	        insert_sort(ddk,L);
            ddk.resize(k);
            L.resize(k);
        	poor=ddk[k-1];
        }
    }
    else{
	    std::vector<T> branch_value;
		std::vector<RTree_node<T> *> branch;
        for(int i = 0; i < node->elements; i++){
	        branch_value.push_back(node->data_internal_node[i].region->distance_geometric(q));
            branch.push_back(node->data_internal_node[i].child);
        }
        insert_sort(branch_value,branch);
		for(int i = 0; i < node->elements; i++){
  	        if(branch_value[i]<=poor){
                DFS_recursive(q,k,branch[i],L,ddk,poor);
            }
        }
      }
}

template<typename T>
void RTree<T>::insert_sort(std::vector<T> & dtmp, std::vector<d_leaf<T>*> & chld){
    for(int i = 0;i < dtmp.size();i++){
		float cur_value = dtmp[i];
		d_leaf<T> * cur_child = chld[i];
		int j = i -1;
		while(j>=0 && dtmp[j] > cur_value){
			dtmp[j+1]= dtmp[j];
			chld[j+1] = chld[j];
			j--;
		}
    	dtmp[j+1] = cur_value;
		chld[j+1] = cur_child;
	}
}

template<typename T>
void RTree<T>::insert_sort(std::vector<T> & dtmp, std::vector<RTree_node<T>*> & chld){
    for(int i = 0;i < dtmp.size();i++){
		float cur_value = dtmp[i];
		RTree_node<T> * cur_child = chld[i];
		int j = i -1;
		while(j>=0 && dtmp[j] > cur_value){
			dtmp[j+1]= dtmp[j];
			chld[j+1] = chld[j];
			j--;
		}
    	dtmp[j+1] = cur_value;
		chld[j+1] = cur_child;
	}
}

template<typename T>
void RTree<T>::k_NN_DFS(Point<T> q, int k, std::vector<d_leaf<T>*> &L){
    std::vector<T> dk;
    T poor = std::numeric_limits<T>::max();
    
    if(this->root != nullptr)
        DFS_recursive(q, k, this->root,L,dk,poor);
}

template<typename T>
void RTree<T>::getAll_values(RTree_node<T> * node,std::vector<Polygon<T> *> & ans)
{
    if(node != nullptr){
        if(!node->is_leaf){
            for(int i = 0; i < node->elements;i++){
                getAll_values(node->data_internal_node[i].child,ans);
            }
        }
        else{
            for(int i = 0; i < node->elements; i++){
                ans.push_back(node->data_leafs[i].polygon);
            }
        }
    }
}

template<typename T>
std::string RTree<T>::show_values_JSON()
{
    std::string json = "";

    if(this->root != nullptr)
        showAll_values_JSON(this->root, 0, json);
    if(json.length() != 0)
        json.erase(json.length()-1);
    return "["+json+"]";
}

template<typename T>
void RTree<T>::showAll_values_JSON(RTree_node<T> *node, int level, std::string &json)
{
    if(!node->is_leaf)
    {
        for(int i=0; i<node->elements; i++)
        {
            json +="{";
            json += "\"level\":"+std::to_string(node->at_level)+",";
            json += "\"is_leaf\":"+std::to_string(0)+",";
            json += "\"key\":"+std::to_string(node->data_internal_node[i].region->get_key())+",";
            json += "\"elements\":[";
            json +="["+std::to_string(node->data_internal_node[i].region->get_Pmin().get_X())+","+std::to_string(node->data_internal_node[i].region->get_Pmin().get_Y())+"],";
            json +="["+std::to_string(node->data_internal_node[i].region->get_Pmax().get_X())+","+std::to_string(node->data_internal_node[i].region->get_Pmax().get_Y())+"]";
            json +="]}";
            json +=",\n";
            showAll_values_JSON(node->data_internal_node[i].child, level+1, json);

        }
    }
    else
    {
        for(int i=0; i<node->elements; i++)
        {
            json +="{";
            json += "\"level\":"+std::to_string(node->at_level)+",";
            json += "\"is_leaf\":"+std::to_string(1)+",";
            json += "\"elements\":[";
            std::vector<Point<T>> vertices = node->data_leafs[i].polygon->get_vertices();
            json += "[";
            for(int j=0; j<node->data_leafs[i].polygon->get_vertices().size(); j++)
            {
                json +="["+std::to_string(vertices[j].get_X())+","+std::to_string(vertices[j].get_Y())+"]";
            }
            json += "]";
            json +="]}";
            json +=",\n";
        }

        // std::cout << json <<std::endl;
    }

}

template<typename T>
 std::string RTree<T>::get_regions_JSON(RTree_node<T>* node){
    std::string json = "";
    for(int i=0; i<node->elements; i++)
    {
        json +="{";
        json += "\"level\":"+std::to_string(node->at_level)+",";
        json += "\"is_leaf\":"+std::to_string(0)+",";
        json += "\"key\":"+std::to_string(node->data_internal_node[i].region->get_key())+",";
        json += "\"elements\":[";
        json +="["+std::to_string(node->data_internal_node[i].region->get_Pmin().get_X())+","+std::to_string(node->data_internal_node[i].region->get_Pmin().get_Y())+"],";
        json +="["+std::to_string(node->data_internal_node[i].region->get_Pmax().get_X())+","+std::to_string(node->data_internal_node[i].region->get_Pmax().get_Y())+"]";
        json +="]}\n";

    }
    return json;
}

template<typename T>
void RTree<T>::get_polygons_JSON(const std::vector<d_leaf<T>*> & ans,std::string &json){
    json += "[";
    for(int i=0; i<ans.size(); i++){
        json +="[";    
        for(int j = 0; j < ans[i]->polygon->corners;j++){
           json +="[";
           json += std::to_string(ans[i]->polygon->get_vertices()[j].get_X());
           json += ",";
           json += std::to_string(ans[i]->polygon->get_vertices()[j].get_Y()); 
           json +="]";
           if((j+1) != ans[i]->polygon->corners)
                json +=",";
        }
        json +="]";
        if((i+1) != ans.size())
            json +=",";
    }
    json += "]";
}

template<typename T>
void RTree<T>::get_Range_Search_JSON(const std::vector<Polygon<T> *> & data, std::string &json){
    json += "[";
    for(int i=0; i<data.size(); i++){
        json +="[";
        json +="{";
        json +="\"elements\":[";    
        auto dat = data[i]->get_vertices();
        for(int j = 0; j < dat.size(); j++){
            json +="[";  
            
            json += std::to_string(dat[j].get_X());
            json += ",";
            json += std::to_string(dat[j].get_Y()); 
        
            json +="]";
            if((j+1) != dat.size())
                    json +=",";
        }            
    
        json +="],";
        json+="\"level\":"+std::to_string(0);
        json += "}";

        json +="]";
        if((i+1) != data.size())
            json +=",\n";
    }
    json += "]";
}

template<typename T>
int RTree<T>::get_items_cnt(){
    return this->items_cnt;
}

// template<typename T>
// RTree<T>::~RTree(){
//     if(this->root != nullptr)
//         delete this->root;
// }

template<typename T>
RTree_node<T>::~RTree_node(){
    if(!this->is_leaf){
        for(int i = 0; i < this->elements; i++){
            delete this->data_internal_node[i].child;
            this->data_internal_node[i].child = nullptr;
            delete this->data_internal_node[i].region;
            this->data_internal_node[i].region = nullptr;
        }
    }
    else{
        for(int i = 0; i < this->elements; i++){
            if(this->data_leafs[i].polygon!=nullptr){
                delete this->data_leafs[i].polygon;
                this->data_leafs[i].polygon = nullptr;
                delete this->data_leafs[i].region;
                this->data_leafs[i].region = nullptr;
            }
        }
    }
}


#endif