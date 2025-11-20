#include <memory>
#include "MinKnapsack.h"
#include <stdio.h>
#include "Point2D.h"
#include "UnionFind.cpp"

#define CIRCLE_CENTER_EPSILON 1.0e-7

long nR = 0;
double total;

enum direction{
    PLUS,
    MINUS,
    EQUAL
};

struct Special_vertex {
    Voronoi::NewDiagram::HalfEdgePtr edgein;
    Voronoi::NewDiagram::HalfEdgePtr edgeout;
    Voronoi::NewDiagram::FacePtr reg; //region at the right of the pending edge
    Voronoi::NewDiagram::SitePtr eti_minus;
    Voronoi::NewDiagram::SitePtr eti_plus;
    double d_minus;
    double d_plus;
    direction direz;
};

// Functions to help delete references
void clearSpecialVertexVector(std::vector<Special_vertex>& vec) {
    for (auto& sv : vec) {
        // --- Clean up edgein ---
        if (sv.edgein) {
            sv.edgein->twin.reset();
            sv.edgein->next.reset();
            sv.edgein->region.reset();
            sv.edgein->tail.reset();
            sv.edgein->head.reset();
        }

        // --- Clean up edgeout ---
        if (sv.edgeout) {
            sv.edgeout->twin.reset();
            sv.edgeout->next.reset();
            sv.edgeout->region.reset();
            sv.edgeout->tail.reset();
            sv.edgeout->head.reset();
        }

        // --- Clean up region ---
        if (sv.reg) {
            // Optionally reset edges inside region if needed
            sv.reg.reset();
        }

        // --- Clean up site pointers ---
        sv.eti_minus.reset();
        sv.eti_plus.reset();
    }

    // Finally clear the vector itself
    vec.clear();
}

void clearHalfEdgeVector(std::vector<Voronoi::NewDiagram::HalfEdgePtr>& vec) {
    for (auto& he : vec) {
        if (he) {
            // Break cycles
            he->twin.reset();
            he->next.reset();
            he->region.reset();
            he->tail.reset();
            he->head.reset();
            he->label.reset();
        }
    }
    vec.clear();
}

void clearFaceVector(std::vector<Voronoi::NewDiagram::FacePtr>& vec) {
    // Faces themselves don't own edges or other cycles.
    // Just reset the shared_ptrs.
    for (auto& f : vec) {
        f.reset();
    }
    vec.clear();
}



Voronoi::NewDiagram::FacePtr createRegion(std::list<Voronoi::NewDiagram::FacePtr>& R,Voronoi::NewDiagram::FacePtr&  r, Voronoi::NewDiagram::SitePtr lambda, Voronoi::NewDiagram::HalfEdgePtr e) {
    nR++;
    Voronoi::NewDiagram::FacePtr t_ptr = std::make_shared<Voronoi::NewDiagram::Face>();
    auto& t = *t_ptr;
    t.ID = nR;
    t.flag = 1;
    t.firstEdge = e;
    if ((r->weight + lambda->capacity)>total) {
        t.sites = r->sites;
        t.weight = r->weight;
        t.pivot = lambda;
    }
    else {
        t.sites = r->sites;
        t.sites.push_back(lambda);
        t.weight = r->weight+lambda->capacity;
        t.pivot = nullptr;
    }
    R.push_back(t_ptr);
    std::cout << "I added the region: " << R.back()->ID <<"\n";
    return t_ptr;
}

void createPendingEdge(Voronoi::NewDiagram::HalfEdgePtr e, Voronoi::NewDiagram::HalfEdgePtr f) {
    auto a = std::make_shared<Voronoi::NewDiagram::HalfEdge>();
    auto b = std::make_shared<Voronoi::NewDiagram::HalfEdge>();
    a->tail = e->head;
    a->head = Voronoi::NewDiagram::Vertex::getNullVertex(); //??? when i put a real vertex i have to be careful to not change the value of the static vertex!!!!!!!!!!!!
    b->tail = Voronoi::NewDiagram::Vertex::getNullVertex(); //???
    b->head = f->tail;
    a->label = e->twin->label;
    b->label = f->twin->label;
    a->twin = b;
    b->twin = a;
    b->next = f;
    e->next = a;
}


Point2D circumcenter(Point2D p1, Point2D p2, Point2D p3) {
    Point2D u1 = (p1 - p2).normalized(), u2 = (p3 - p2).normalized();

    double cross = crossProduct(u1, u2);

    // check if vectors are collinear
    if (fabs(cross) < CIRCLE_CENTER_EPSILON) {
        return false;
    }

    // get cental points
    Point2D pc1 = 0.5 * (p1 + p2), pc2 = 0.5 * (p2 + p3);

    // get free components
    double b1 = dotProduct(u1, pc1), b2 = dotProduct(u2, pc2);

    Point2D center;

    // calculate the center of a circle
    center.x = (b1 * u2.y - b2 * u1.y) / cross;
    center.y = (u1.x * b2 - u2.x * b1) / cross;

    return center;
}

Point2D circumcenter(std::vector<Voronoi::NewDiagram::SitePtr> s) {
    return circumcenter(s.at(0)->point, s.at(1)->point, s.at(2)->point);
}

double dist(Point2D p1, Point2D p2) {
    return sqrt(pow(p1.x - p2.x,2) + pow(p1.y - p2.y,2));
}

bool pointInsideRegion(Point2D p, Point2D p1, Point2D p2) {
    return (((p.x - p1.x) * (p2.y - p1.y)) 
        >= ((p.y - p1.y) * (p2.x - p1.x)));
}


void chooseDirection(Special_vertex* r) {
    r->direz = EQUAL;
    if ((1/r->d_minus)>(1/r->d_plus)) {
        r->direz = MINUS;
    }
    if ((1 / r->d_minus) < (1 / r->d_plus)) {
        r->direz = PLUS;
    }
}



int circumcenterInside(size_t& r_index, std::vector<Special_vertex>& L, int E, bool open) {
    auto v_ptr = std::make_shared<Voronoi::NewDiagram::Vertex>();
    auto& v = *v_ptr;

    size_t next_index = (r_index + 1) % L.size();
    size_t next_next_index = (r_index + 2) % L.size();

    auto& r = L[r_index];
    auto& next = L[next_index];
    auto& next_next = L[next_next_index];

    v.triplet = std::vector{ next.eti_minus, next.eti_plus, next_next.eti_plus };
    v.point = circumcenter(v.triplet);
    next.edgein->next->head = v_ptr;
    next.edgein->next->twin->tail = v_ptr;
    next_next.edgein->next->head = v_ptr;
    next_next.edgein->twin->tail = v_ptr;

    if (E == 3) {
        if (!open) {
            r.edgein->next->head = v_ptr;
            r.edgein->next->twin->tail = v_ptr;
        }
        else {
            createPendingEdge(next.edgein->next, next_next.edgein->next->twin);
            r.edgein->next = next.edgein->next->next->twin;
            next.edgein->next->next = r.edgeout;
            r.edgein->next->tail->infinite = true;
            r.edgein->next->twin->head->infinite = true;
        }
        next.edgein->next->next = r.edgein->next->twin;
        next_next.edgein->next->next = next.edgein->next->twin;
        r.edgein->next->next = next_next.edgein->next->twin;

        L.erase(L.begin() + r_index);
        L.erase(L.begin() + next_index - 1); // Adjust index after erasure
        L.erase(L.begin() + next_next_index - 2);
        E = 0;

    }else{
        next_next.edgein->next->next = next.edgein->next->twin;
        createPendingEdge(next.edgein->next, next_next.edgein->next->twin);

        Special_vertex p = Special_vertex{};
        p.edgein = next.edgein->next;
        p.edgeout = next_next.edgein->next->twin;
        p.eti_minus = p.edgein->label;
        p.eti_plus = p.edgeout->label;

        L.erase(L.begin() + next_index);
        L.erase(L.begin() + next_next_index - 1);
        size_t p_index = r_index + 1;
        L.insert(L.begin() + p_index, p);

        E--;

        Point2D P = circumcenter(r.eti_minus->point, p.eti_minus->point, p.eti_plus->point);
        if (pointInsideRegion(P, r.edgein->head->point, p.edgeout->head->point)) {
            r.d_plus = dist(P, r.edgeout->tail->point);
            p.d_minus = dist(P, p.edgein->head->point);
        }
        else {
            r.d_plus = -dist(P, r.edgeout->tail->point);
            p.d_minus = -dist(P, p.edgein->head->point);
        }

        Point2D P1 = circumcenter(p.eti_minus->point, p.eti_plus->point, L[(p_index + 1) % L.size()].eti_minus->point);
        if (pointInsideRegion(P1, p.edgein->head->point, L[(p_index + 1) % L.size()].edgeout->head->point)) {
            p.d_plus = dist(P, p.edgeout->tail->point);
            L[(p_index + 1) % L.size()].d_minus = dist(P, L[(p_index + 1) % L.size()].edgein->head->point);
        }
        else {
            p.d_plus = -dist(P, p.edgeout->tail->point);
            L[(p_index + 1) % L.size()].d_minus = -dist(P, L[(p_index + 1) % L.size()].edgein->head->point);
        }
        chooseDirection(&r);
        chooseDirection(&p);
        chooseDirection(&L[(p_index + 1) % L.size()]);
    }
    return E;
}


void partition(Voronoi::NewDiagram::FacePtr& r_ptr, std::list<Voronoi::NewDiagram::FacePtr>& R) {
    auto& r = *r_ptr;
    auto e = std::vector<Voronoi::NewDiagram::HalfEdgePtr>();
    auto lambda = std::vector<Voronoi::NewDiagram::SitePtr>();
    auto reg = std::vector<Voronoi::NewDiagram::FacePtr>();

    e.push_back(r.firstEdge);
    lambda.push_back(e.at(0)->twin->label);
    reg.push_back(createRegion(R, r_ptr, lambda.at(0), e.at(0)));

    while (e.at(0)->next != r.firstEdge && e.at(0)->next->twin->label == lambda.at(0)) {
        e.at(0) = e.at(0)->next;
    }

    if (e.at(0)->next->twin->label != lambda.at(0)) {
        lambda.push_back(e.at(0)->next->twin->label);
        e.push_back(e.at(0)->next);
        createPendingEdge(e.at(0), e.at(1));
        reg.push_back(createRegion(R, r_ptr, lambda.at(1), e.at(0)->next->twin));

        while (e.at(1)->next->twin->label == lambda.at(1)) {
            e.at(1) = (e.at(1))->next;
        }

        if ((e.at(1))->next->twin->label == lambda.at(0)) {
            e.at(0)->next->head = e.at(1)->head;
            e.at(0)->next->twin->tail = e.at(1)->head;
            e.at(0)->next->next = e.at(1)->next;
            e.at(1)->next = e.at(0)->next->twin;
        }
        else {
            int k = 1;
            lambda.push_back(e.at(1)->next->twin->label);
            e.push_back(e.at(1)->next);
            createPendingEdge(e.at(1), e.at(2));
            bool open = false;

            do {
                reg.push_back(createRegion(R, r_ptr, lambda.at(k + 1), e.at(k)->next->twin));
                k++;
                while (e.at(k)->next->twin->label == lambda.at(k)) {
                    e.at(k) = e.at(k)->next;
                }
                lambda.push_back(e.at(k)->next->twin->label);
                e.push_back(e.at(k)->next);

                if (!(e.at(k)->head->infinite)) {
                    createPendingEdge(e.at(k), e.back());
                }
                else {
                    open = true;
                }
            } while (lambda.back() != lambda.front());

            // Partition in 3 parts (first part)
            std::vector<Special_vertex> L{};
            int E = 0;
            for (int h = 0; h < k; h++) {
                ++E;
                Special_vertex rv{};
                rv.eti_minus = lambda.at(h);
                rv.eti_plus = lambda.at(h % (k + 1));
                rv.reg = reg.at(h % (k + 1));
                rv.edgein = e.at(h);

                if (rv.edgein->head->infinite) {
                    rv.edgeout = e.at(h)->next;
                }
                else {
                    rv.edgeout = e.at(h)->next->twin->next;
                }
                L.push_back(rv);
            }

            for (size_t l_iter = 0; l_iter < L.size(); ++l_iter) {
                size_t l_next = (l_iter + 1) % L.size();
                Point2D P = circumcenter(L[l_iter].eti_minus->point, L[l_iter].eti_plus->point, L[l_next].eti_plus->point);

                if (pointInsideRegion(P, L[l_iter].edgein->head->point, L[l_next].edgein->head->point)) {
                    if (!L[l_iter].edgeout->tail->infinite) {
                        L[l_iter].d_plus = dist(P, L[l_iter].edgeout->tail->point);
                    }
                    if (!L[l_next].edgeout->tail) {
                        L[l_next].d_minus = dist(P, L[l_next].edgeout->tail->point);
                    }
                }
                else {
                    if (!L[l_iter].edgeout->tail->infinite) {
                        L[l_iter].d_plus = -dist(P, L[l_iter].edgeout->tail->point);
                    }
                    if (!L[l_next].edgeout->tail) {
                        L[l_next].d_minus = -dist(P, L[l_next].edgeout->tail->point);
                    }
                }
            }

            for (auto& vertex : L) {
                chooseDirection(&vertex);
            }

            // Partition in 3 parts (second part)
            size_t l_last_iter = L.size() - 1;
            do {
                while ((((l_last_iter + 2) % L.size()) != 0) &&
                    (L[(l_last_iter + 1) % L.size()].direz != MINUS || L[(l_last_iter + 2) % L.size()].direz != PLUS || L[(l_last_iter + 1) % L.size()].d_plus < 0)) {
                    l_last_iter = (l_last_iter + 1) % L.size();
                }
                // Check if there exists a circumcenter inside
                if (L[(l_last_iter + 1) % L.size()].direz != MINUS &&
                    L[(l_last_iter + 2) % L.size()].direz != PLUS &&
                    L[(l_last_iter + 1) % L.size()].d_plus >= 0 &&
                    !(L[(l_last_iter + 1) % L.size()].edgein->head->infinite)) {

                    E = circumcenterInside(l_last_iter, L, E, open);
                }
                else {
                    for (size_t l_iter = 0; l_iter < L.size(); ++l_iter) {
                        L[(l_iter + 1) % L.size()].edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
                        L[l_iter].edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
                        L[(l_iter + 1) % L.size()].edgein->next->next = L[l_iter].edgein->next->twin;
                        L[l_iter].reg->firstEdge = L[l_iter].edgein->next->twin;
                    }
                    L[0].edgein->next->twin->tail->infinite = true;
                    L[1].edgein->next = L[0].edgein->next->twin;
                    L[0].reg->firstEdge = L[0].edgein->next->twin;
                    L[1].edgein->next->head->infinite = true;
                    L[1].edgein->next->next = L[0].edgeout;

                    // I delete the list making sure to not leave dandling references
                    clearSpecialVertexVector(L);
                    E = 0;
                }
            } while (E != 0);
        }
    }

}


void build_minKnapsack(Voronoi::NewDiagram& diagram, std::vector<std::pair<Point2D, double>>& points, double capacity)
{
	std::list<Voronoi::NewDiagram::FacePtr>& R = diagram.getFaces();
    total = capacity;
    nR = R.size() - 1;
    std::list<Voronoi::NewDiagram::FacePtr>::iterator it = R.begin();
    // I keep examining regions until I arrive to the end of the list
    // If I comment the while I test a single iteration!
    //while(it != R.end()){
        long last = nR;
        // I partition all the new regions and add new ones to the end of the list
        while (it != R.end() && (*it)->ID <= last) {
            std::cout << "I examine the region: " << (*it)->ID << " with pivot: "<< (*it)->pivot << "\n";
            //std::cout << "I examine: " << *(*it);
            //std::cout << "Point: " << ((*it)->sites).at(0)->point << "\n";
            if ((*it)->flag == 1 && (*it)->pivot == nullptr && ((*it)->weight < capacity)) {
                std::cout << "I divide the region: " << (*it)->ID << "\n";
                partition(*it, R);
                (*it)->flag = 0;
                // I erase the current region and pass to the next (putting the flag to 0 wasn't really necessary)
                it = R.erase(it);
            }
            else {
                ++it;
            }
        }
        std::list<Voronoi::NewDiagram::FacePtr>::iterator& firstNewRegion = it;

        // Create structure Union-Find
        std::vector<int> UF(nR-((*firstNewRegion)->ID));

        // All new edges created during the partitions must be assigned their region
        // I create a list of new edges (ListE)
        auto E = std::list<Voronoi::NewDiagram::HalfEdgePtr>();

        // I iterate over all new regions (it is already pointing to the first new region)
        while (it != R.end()){
            // I insert each new region in the Union-Find structure and each of their edges in the ListE
            // REMARK: this list will contain both new edges and old edges of the previous regions
            UF[nR-((*it)->ID)]; //TODO
            Voronoi::NewDiagram::HalfEdgePtr e = (*it)->firstEdge;
            do{
                e->region;
                e = e->next;
                E.push_back(e);
            }while(e!=(*it)->firstEdge);
            ++it;
        }

        std::list<Voronoi::NewDiagram::HalfEdgePtr>::iterator e = E.begin();
        while (e != E.end()){
            // I iterate over ListE: if an edge is actually a new edge (divides different regions) then I take it out of ListE
            // This means that ListE will contain all edges to be eliminated
            if((*e)->region->sites!=(*e)->twin->region->sites || (*e)->region->pivot>0 || (*e)->twin->region->pivot >0){
                e = E.erase(e);
            }else{
                // If an edge is NOT a new edge it means its neighbouring regions are actually the same region: I sign it as to be eliminated and apply Union-Find
                (*e)->region->flag=0;
                UF; //TODO merge together the regions
                ++e;
            }
        }

        int lastNewRegion = nR;
        for (int i =0; i< (lastNewRegion-(*firstNewRegion)->ID);++i){
            // I iterate over all new regions and use Union-Find to find the first of all components that form a same new region
            if (UF[i]==UF[i]){ //TODO fix find and next
                auto r = std::make_shared<Voronoi::NewDiagram::FacePtr>(); // TODO
                Voronoi::NewDiagram::FacePtr t_ptr = std::make_shared<Voronoi::NewDiagram::Face>();
                auto& T = *t_ptr;
                nR++;
                T.ID = nR;
                T.flag = 1;
                T.weight = (*r)->weight; //TODO
                T.sites = (*r)->sites; //TODO
                auto newRegionFirstEdge = Voronoi::NewDiagram::HalfEdgePtr(nullptr);
                do{
                    auto e = (*r)->firstEdge;
                    do
                    {
                        (*e).region = t_ptr; //TODO understand if t_ptr is pointing to updated version of T
                        if(newRegionFirstEdge==nullptr || e->tail==Voronoi::NewDiagram::Vertex::getNullVertex()){
                            newRegionFirstEdge=e;
                        }
                    } while (e!=nullptr);
                    
                } while (UF[i]); //TODO iterate until end of list
                

            }
            // If a region is composed of only one component I don't need to do anything,
            // otherwise I create a new record for the first component region: this will represent the whole fused region
            // I insert this region in the list of regions ListR

        }
        
        // I iterate over ListE and update their connections so that, when they wil be eliminated, all other edges will be correctly connected

        // I delete all edges in ListE

        it = firstNewRegion;
        // I delete regions that have been fused together
        while (it != R.end()) {
            if ((*it)->flag == 0) { 
                it = R.erase(it); 
            }
            else {
                ++it; 
            }
        }
        it = firstNewRegion;
    //}
}


// Now inside the subfunctions I am using vectors for simplicity, the commented versions are the ones that used lists instead!

//int circumcenterInside(std::list<Special_vertex*>::iterator& r_it,
//    std::list<Special_vertex*>& L, int E, bool open) {
//    
//    auto v_ptr = std::make_shared<Voronoi::NewDiagram::Vertex>();
//    auto& v = *v_ptr;
//    auto nextNext_it = std::next(r_it);
//    auto next = *nextNext_it;
//    auto r = *r_it;
//    auto next_it = std::next(r_it, 2);
//    auto nextNext = *next_it;
//    v.triplet = std::vector{ next->eti_minus, next->eti_plus, nextNext->eti_plus };
//    v.point = circumcenter(v.triplet);
//    next->edgein->next->head = v_ptr;
//    next->edgein->next->twin->tail = v_ptr;
//    nextNext->edgein->next->head = v_ptr;
//    nextNext->edgein->twin->tail = v_ptr;
//
//    if (E == 3) {
//        if (!open) {
//            r->edgein->next->head = v_ptr;
//            r->edgein->next->twin->tail = v_ptr;
//        }
//        else {
//            createPendingEdge(next->edgein->next, nextNext->edgein->next->twin);
//            r->edgein->next = next->edgein->next->next->twin;
//            next->edgein->next->next->next = r->edgeout;
//            //r->edgein->next->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
//            r->edgein->next->tail->infinite = true;
//            //r->edgein->next->twin->head = Voronoi::NewDiagram::Vertex::getNullVertex();
//            r->edgein->next->twin->head->infinite = true;
//        }
//        // Connect the three remaining pending edges
//        next->edgein->next->next = r->edgein->next->twin;
//        nextNext->edgein->next->next = next->edgein->next->twin;
//        r->edgein->next->next = nextNext->edgein->next->twin;
//        L.erase(r_it);
//        L.erase(next_it);
//        L.erase(nextNext_it);
//        E = 0;
//    }
//    else {
//        nextNext->edgein->next->next = next->edgein->next->twin;
//        createPendingEdge(next->edgein->next, nextNext->edgein->next->twin);
//
//        Special_vertex p = Special_vertex{};
//        p.edgein = next->edgein->next;
//        p.edgeout = nextNext->edgein->next->twin;
//        p.eti_minus = p.edgein->label;
//        p.eti_plus = p.edgeout->label;
//
//        // Remove the records next and nextNext from l and insert p in their place
//        L.erase(next_it);
//        L.erase(nextNext_it);
//        auto p_it = L.insert(std::next(r_it), &p);
//
//        E--;
//
//        // Update distances and directions
//        Point2D P = circumcenter(r->eti_minus->point, p.eti_minus->point, p.eti_plus->point);
//        if (pointInsideRegion(P, r->edgein->head->point, p.edgeout->head->point)) {
//            r->d_plus = dist(P, r->edgeout->tail->point);
//            p.d_minus = dist(P, p.edgein->head->point);
//        }
//        else {
//            r->d_plus = -dist(P, r->edgeout->tail->point);
//            p.d_minus = -dist(P, p.edgein->head->point);
//        }
//        Point2D P1 = circumcenter(p.eti_minus->point, p.eti_plus->point, (*std::next(p_it))->eti_minus->point);
//        if (pointInsideRegion(P1, p.edgein->head->point, (*std::next(p_it))->edgeout->head->point)) {
//            p.d_plus = dist(P, p.edgeout->tail->point);
//            (*std::next(p_it))->d_minus = dist(P, (*std::next(p_it))->edgein->head->point);
//        }
//        else {
//            p.d_plus = -dist(P, p.edgeout->tail->point);
//            (*std::next(p_it))->d_minus = -dist(P, (*std::next(p_it))->edgein->head->point);
//        }
//        chooseDirection(r);
//        chooseDirection(&p);
//        chooseDirection(*std::next(p_it));
//    }
//    return E;
//}


//void partition(Voronoi::NewDiagram::FacePtr& r_ptr, std::list<Voronoi::NewDiagram::FacePtr>& R) {
//    auto& r = *r_ptr;
//    auto e = std::vector<Voronoi::NewDiagram::HalfEdgePtr>();
//    auto lambda = std::vector<Voronoi::NewDiagram::SitePtr>();
//    auto reg = std::vector<Voronoi::NewDiagram::FacePtr>();
//    
//    e.push_back(r.firstEdge);
//    lambda.push_back(e.at(0)->twin->label);
//    reg.push_back(createRegion(R, r_ptr, lambda.at(0), e.at(0)));
//
//    
//    while (e.at(0)->next != r.firstEdge && e.at(0)->next->twin->label == lambda.at(0)) {
//        e.at(0) = e.at(0)->next;
//    }
//
//    if (e.at(0)->next->twin->label != lambda.at(0)) {
//        lambda.push_back(e.at(0)->next->twin->label);
//        e.push_back(e.at(0)->next);
//        createPendingEdge(e.at(0), e.at(1));
//        reg.push_back(createRegion(R, r_ptr, lambda.at(1), e.at(0)->next->twin));
//        
//        while (e.at(1)->next->twin->label == lambda.at(1)) {
//            e.at(1) = (e.at(1))->next;
//        }
//
//        if ((e.at(1))->next->twin->label == lambda.at(0)) {
//            e.at(0)->next->head = e.at(1)->head;
//            e.at(0)->next->twin->tail = e.at(1)->head;
//            e.at(0)->next->next = e.at(1)->next;
//            e.at(1)->next = e.at(0)->next->twin;
//        }
//        else {
//            int k = 1;
//            lambda.push_back(e.at(1)->next->twin->label);
//            e.push_back(e.at(1)->next);
//            createPendingEdge(e.at(1),e.at(2));
//            bool open = false;
//
//            do {
//                reg.push_back(createRegion(R, r_ptr, lambda.at(k+1), e.at(k)->next->twin));
//                k++;
//                while (e.at(k)->next->twin->label == lambda.at(k)) {
//                    e.at(k) = e.at(k)->next;
//                }
//                lambda.push_back(e.at(k)->next->twin->label);
//                e.push_back(e.at(k)->next);
//
//                if (!(e.at(k)->head->infinite)) {
//                    createPendingEdge(e.at(k), e.back());
//                }
//                else {
//                    open = true;
//                }
//            } while (lambda.back() != lambda.front());
//
//            // Partition in 3 parts
//            std::list<Special_vertex*> L{};
//            int E = 0;
//            for (int h = 0; h < k; h++) {
//                ++E;
//                Special_vertex rv{};
//                rv.eti_minus = lambda.at(h);
//                rv.eti_plus = lambda.at(h%(k+1));
//                rv.reg = reg.at(h % (k + 1));
//                rv.edgein = e.at(h);
//
//                if (rv.edgein->head->infinite) {
//                    rv.edgeout = e.at(h)->next;
//                }
//                else {
//                    rv.edgeout = e.at(h)->next->twin->next;
//                }
//                L.push_back(&rv);
//            }
//            // to make the list circular
//            auto temp = L.end();
//            temp--;
//            temp._Ptr->_Next = L.begin()._Ptr; // Next Node of the last elemen is now first elemen of the List
//            L.begin()._Ptr->_Prev = temp._Ptr; // Prev Node of the first element is now Last node
//
//
//            auto l_iter = L.begin();
//            do{
//                auto l_next = std::next(l_iter);
//                Point2D P = circumcenter((*l_iter)->eti_minus->point, (*l_iter)->eti_plus->point, (*l_next)->eti_plus->point);
//
//                if (pointInsideRegion(P, (*l_iter)->edgein->head->point, (*l_next)->edgein->head->point)) {
//                    if (!(*l_iter)->edgeout->tail->infinite) {
//                        (*l_iter)->d_plus = dist(P, (*l_iter)->edgeout->tail->point);
//                    }
//                    if (!(*l_next)->edgeout->tail) {
//                        (*l_next)->d_minus = dist(P, (*l_next)->edgeout->tail->point);
//                    }
//                }
//                else {
//                    if (!(*l_iter)->edgeout->tail->infinite) {
//                        (*l_iter)->d_plus = -dist(P, (*l_iter)->edgeout->tail->point);
//                    }
//                    if (!(*l_next)->edgeout->tail) {
//                        (*l_next)->d_minus = -dist(P, (*l_next)->edgeout->tail->point);
//                    }
//                }
//                ++l_iter;
//            } while (l_iter==L.begin());
//
//            do{
//                chooseDirection(*l_iter);
//                ++l_iter;
//            } while (l_iter == L.begin());
//
//            // Partition in 3 parts (second part)
//            auto l_last_iter = std::prev(L.end());
//            do {
//                while ((((std::distance(L.begin(), l_last_iter) + 2) % L.size()) != 0) &&
//                    ((*std::next(l_last_iter))->direz != MINUS || (*std::next(l_last_iter, 2))->direz != PLUS || (*std::next(l_last_iter))->d_plus < 0)) {
//                    ++l_last_iter;
//                }
//                // Check if there exists a circumcenter inside
//                if ((*std::next(l_last_iter))->direz != MINUS &&
//                    (*std::next(l_last_iter, 2))->direz != PLUS &&
//                    (*std::next(l_last_iter))->d_plus >= 0 &&
//                    !((*std::next(l_last_iter))->edgein->head->infinite)) {
//
//                    E = circumcenterInside(l_last_iter, L, E, open);
//                } else {
//                    l_iter = L.begin();
//                    do{
//                         (*std::next(l_iter))->edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
//                        //(*std::next(l_iter))->edgein->next->twin->tail->infinite = true;
//                         (*l_iter)->edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
//                        //(*l_iter)->edgein->next->twin->tail->infinite = true;
//                         (*std::next(l_iter))->edgein->next->next = (*l_iter)->edgein->next->twin;
//                         (*l_iter)->reg->firstEdge = (*l_iter)->edgein->next->twin;
//                         ++l_iter;
//                    } while (std::next(l_iter,2) != L.begin());
//                    (*l_iter)->edgein->next->twin->tail->infinite = true;
//                    (*std::next(l_iter))->edgein->next = (*l_iter)->edgein->next->twin;
//                    (*l_iter)->reg->firstEdge = (*l_iter)->edgein->next->twin;
//                    l_iter = std::next(l_iter); //?
//                    (*std::next(l_iter))->edgein->next->head->infinite = true;
//                    (*std::next(l_iter))->edgein->next->next = (*l_iter)->edgeout;
//                    L.clear();
//                    E = 0;
//                }
//            } while (E != 0);
//        }
//    }
//}