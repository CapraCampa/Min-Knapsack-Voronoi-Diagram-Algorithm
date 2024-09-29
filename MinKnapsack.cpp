
#include "MinKnapsack.h"
#include <stdio.h>

long nR = 0;
double total;

enum direction{
    PLUS,
    MINUS,
    EQUAL
};

struct Special_vertex {
    Voronoi::NewDiagram::HalfEdge* edgein;
    Voronoi::NewDiagram::HalfEdge* edgeout;
    Voronoi::NewDiagram::Face* reg; //region at the right of the pending edge
    Voronoi::NewDiagram::Site* eti_minus;
    Voronoi::NewDiagram::Site* eti_plus;
    double d_minus;
    double d_plus;
    direction direz;
};

Voronoi::NewDiagram::Face& createRegion(std::list<Voronoi::NewDiagram::Face>& R,Voronoi::NewDiagram::Face&  r, Voronoi::NewDiagram::Site* lambda, Voronoi::NewDiagram::HalfEdge* e) {
    nR++;
    Voronoi::NewDiagram::Face t{}; //?????????
    t.ID = nR;
    t.flag = 1;
    t.firstEdge = e;
    if ((r.weight + lambda->capacity)>total) {
        t.sites = r.sites;
        t.weight = r.weight;
        t.pivot = lambda;
    }
    else {
        t.sites = r.sites;
        t.sites.push_back(lambda);
        t.weight = r.weight+lambda->capacity;
        t.pivot = nullptr;
    }
    R.push_back(t);
}

void createPendingEdge(Voronoi::NewDiagram::HalfEdge* e, Voronoi::NewDiagram::HalfEdge* f) {
    Voronoi::NewDiagram::HalfEdge a{};
    Voronoi::NewDiagram::HalfEdge b{};
    a.tail = e->head;
    b.head = f->tail;
    a.label = e->twin->label;
    b.label = f->twin->label;
    a.twin = &b;
    b.twin = &a;
    b.next = f;
    e->next = &a;
}


void partition(Voronoi::NewDiagram::Face& r, std::list<Voronoi::NewDiagram::Face>& R) {
    std::vector<Voronoi::NewDiagram::HalfEdge*> e;
    std::vector<Voronoi::NewDiagram::Site*> lambda;
    std::vector<Voronoi::NewDiagram::Face> reg;
    e.push_back(r.firstEdge);
    lambda.push_back(e.at(0)->twin->label);
    reg.push_back(createRegion(R,r, lambda.at(0), e.at(0)));
    while (e.at(0)->next != r.firstEdge && e.back()->next->twin->label == lambda.back()) {
        e.at(0) = e.at(0)->next;
    }
    if (e.at(0)->next->twin->label != lambda.at(0)) {
        lambda.emplace_back(e.at(0)->next->twin->label);
        e.emplace_back(e.at(0)->next);
        createPendingEdge(e.at(0), e.at(1));
        reg.push_back(createRegion(R,r, lambda.at(1), e.at(1)->next->twin));
        while (e.at(1)->next->twin->label == lambda.at(1)) {
            e.at(1) = e.at(1)->next;
        }
        if (e.at(1)->next->twin->label == lambda.at(0)) {
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
            do
            {
                reg.push_back(createRegion(R,r, lambda.at(k + 1), e.at(k)->next->twin));
                k = k + 1;
                while (e.at(k)->next->twin->label == lambda.at(k)) {
                    e.at(k) = e.at(k)->next;
                }
                lambda.push_back(e.at(k)->next->twin->label);
                e.push_back(e.at(k)->next);
                if (!e.at(k)->head->infinite) {
                    createPendingEdge(e.at(k), e.at(k + 1));
                }
                else {
                    open = true;
                }
            } while (lambda.at(k+1)==lambda.at(1));

            //Partition in 3 parts
            std::list<Special_vertex> l{};
            Special_vertex r{};
            
        }
    }
}

void build_minKnapsack(Voronoi::NewDiagram& diagram, std::vector<std::pair<Point2D, double>>& points, double capacity)
{
	auto& R = diagram.getFaces();
    total = capacity;
    nR = R.size() - 1;
    auto it = R.begin();
    while(it != R.end()){
        long last = nR;
        while (it != R.end() && it->ID <= last) {
            if (it->flag == 1 && it->pivot == nullptr && it->weight < capacity) {
                partition(*it, R);
                it->flag = 0;
            }
            ++it;
        }
        Voronoi::NewDiagram::Face firstNewRegion = *it;



        // Delete the element based on a condition
        if (it->flag == 0) { // Example condition: delete even numbers
            it = R.erase(it); // erase returns the next iterator
        }
        else {
            ++it; // only increment if no deletion
        }
        
    }
}
