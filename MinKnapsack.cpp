
#include "MinKnapsack.h"
#include <stdio.h>
#include "Point2D.h"

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
    return t;
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

//to implement
Point2D circumcenter(Point2D p1, Point2D p2, Point2D p3) {

}

//to implement
Point2D circumcenter(std::vector<Voronoi::NewDiagram::Site*>) {

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

int circumcenterInside(Special_vertex* r, Special_vertex* next, Special_vertex* nextNext, std::list<Special_vertex*>& l, int E, bool open) {
    Voronoi::NewDiagram::Vertex v{};
    v.triplet = std::vector{ next->eti_minus, next->eti_plus, nextNext->eti_plus };
    v.point = circumcenter(v.triplet);
    next->edgein->next->head = &v;
    next->edgein->next->twin->tail = &v;
    nextNext->edgein->next->head = &v;
    nextNext->edgein->twin->tail = &v;

    if (E == 3) {
        if (!open) {
            r->edgein->next->head = &v;
            r->edgein->next->twin->tail = &v;
        }
        else {
            createPendingEdge(next->edgein->next, nextNext->edgein->next->twin);
            r->edgein->next = next->edgein->next->next->twin;
            next->edgein->next->next->next = r->edgeout;
            r->edgein->next->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
            r->edgein->next->twin->head = Voronoi::NewDiagram::Vertex::getNullVertex();
        }
        /* Connect the three remaining pending edges*/
        next->edgein->next->next = r->edgein->next->twin;
        nextNext->edgein->next->next = next->edgein->next->twin;
        r->edgein->next->next = nextNext->edgein->next->twin;
        //l.clear() ???????
        E = 0;
    }
    else {
        nextNext->edgein->next->next = next->edgein->next->twin;
        createPendingEdge(next->edgein->next, nextNext->edgein->next->twin);

        Special_vertex p = Special_vertex{};
        p.edgein = next->edgein->next;
        p.edgeout = nextNext->edgein->next->twin;
        p.eti_minus = p.edgein->label;
        p.eti_plus = p.edgeout->label;

        /* Remove the records next and nextNext from l and insert p in their place */
        l.remove(next);
        l.remove(nextNext);
        E--;

        Point2D P = circumcenter(r->eti_minus->point, p.eti_minus->point, p.eti_plus->point);
        if (pointInsideRegion(P, r->edgein->head->point, p.edgeout->head->point)) {
            r->d_plus = dist(P, r->edgeout->tail->point);
            p.d_minus = dist(P, p.edgein->head->point);
        }
        else {
            r->d_plus = -dist(P, r->edgeout->tail->point);
            p.d_minus = -dist(P, p.edgein->head->point);
        }
    }
}


void partition(Voronoi::NewDiagram::Face& r, std::list<Voronoi::NewDiagram::Face>& R) {
    std::list<Voronoi::NewDiagram::HalfEdge*> e;
    std::list<Voronoi::NewDiagram::Site*> lambda;
    std::list<Voronoi::NewDiagram::Face> reg;

    e.push_back(r.firstEdge);
    lambda.push_back(e.front()->twin->label);
    reg.push_back(createRegion(R, r, lambda.front(), e.front()));

    auto e_iter = e.begin();
    auto lambda_iter = lambda.begin();

    while ((*e_iter)->next != r.firstEdge && (*std::prev(e.end()))->next->twin->label == lambda.back()) {
        *e_iter = (*e_iter)->next;
    }

    if ((*e_iter)->next->twin->label != lambda.front()) {
        lambda.emplace_back((*e_iter)->next->twin->label);
        e.emplace_back((*e_iter)->next);
        createPendingEdge(e.back(), e.back());
        reg.push_back(createRegion(R, r, lambda.back(), e.back()->next->twin));

        auto e_back_iter = std::prev(e.end());
        auto lambda_back_iter = std::prev(lambda.end());

        while ((*e_back_iter)->next->twin->label == *lambda_back_iter) {
            *e_back_iter = (*e_back_iter)->next;
        }

        if ((*e_back_iter)->next->twin->label == lambda.front()) {
            (*e_iter)->next->head = (*e_back_iter)->head;
            (*e_iter)->next->twin->tail = (*e_back_iter)->head;
            (*e_iter)->next->next = (*e_back_iter)->next;
            (*e_back_iter)->next = (*e_iter)->next->twin;
        }
        else {
            int k = 1;
            lambda.push_back((*e_back_iter)->next->twin->label);
            e.push_back((*e_back_iter)->next);
            createPendingEdge(e.back(), std::next(e.back()));
            bool open = false;

            do {
                reg.push_back(createRegion(R, r, lambda.back(), (*std::prev(e.end()))->next->twin));
                k++;
                while ((*std::prev(e.end()))->next->twin->label == *std::prev(lambda.end())) {
                    *std::prev(e.end()) = (*std::prev(e.end()))->next;
                }
                lambda.push_back((*std::prev(e.end()))->next->twin->label);
                e.push_back((*std::prev(e.end()))->next);

                if (!(*std::prev(e.end()))->head->infinite) {
                    createPendingEdge(*std::prev(e.end()), *std::next(std::prev(e.end())));
                }
                else {
                    open = true;
                }
            } while (*std::prev(lambda.end()) == *std::next(lambda.begin()));

            // Partition in 3 parts
            std::list<Special_vertex*> l{};
            int E = 0;
            for (auto e_elem = e.begin(); e_elem != e.end(); ++e_elem) {
                Special_vertex rv{};
                rv.eti_minus = *lambda_iter;
                rv.eti_plus = *std::next(lambda_iter);
                rv.reg = &(*std::next(reg.begin()));
                rv.edgein = *e_elem;

                if ((*e_elem)->head->infinite) {
                    rv.edgeout = (*e_elem)->next;
                }
                else {
                    rv.edgeout = (*e_elem)->next->twin->next;
                }
                l.push_back(&rv);
                ++lambda_iter;
                ++E;
            }

            for (auto l_iter = l.begin(); l_iter != l.end(); ++l_iter) {
                auto l_next = std::next(l_iter);
                Point2D P = circumcenter((*l_iter)->eti_minus->point, (*l_iter)->eti_plus->point, (*l_next)->eti_plus->point);

                if (pointInsideRegion(P, (*l_iter)->edgein->head->point, (*l_next)->edgeout->head->point)) {
                    if (!(*l_iter)->edgeout->tail->infinite) {
                        (*l_iter)->d_plus = dist(P, (*l_iter)->edgeout->tail->point);
                    }
                    if (!(*l_next)->edgeout->tail) {
                        (*l_next)->d_minus = dist(P, (*l_next)->edgeout->tail->point);
                    }
                }
                else {
                    if (!(*l_iter)->edgeout->tail->infinite) {
                        (*l_iter)->d_plus = -dist(P, (*l_iter)->edgeout->tail->point);
                    }
                    if (!(*l_next)->edgeout->tail) {
                        (*l_next)->d_minus = -dist(P, (*l_next)->edgeout->tail->point);
                    }
                }
            }

            for (auto l_iter = l.begin(); l_iter != l.end(); ++l_iter) {
                chooseDirection(*l_iter);
            }

            // Partition in 3 parts (second part)
            auto l_last_iter = std::prev(l.end());
            do {
                while ((((std::distance(l.begin(), l_last_iter) + 2) % l.size()) != 0) &&
                    ((*std::next(l_last_iter))->direz != MINUS || (*std::next(l_last_iter, 2))->direz != PLUS || (*std::next(l_last_iter))->d_plus < 0)) {
                    --l_last_iter;
                }

                if ((*std::next(l_last_iter))->direz != MINUS &&
                    (*std::next(l_last_iter, 2))->direz != PLUS &&
                    (*std::next(l_last_iter))->d_plus >= 0 &&
                    !((*std::next(l_last_iter))->edgein->head->infinite)) {

                    E = circumcenterInside(*l_last_iter, *std::next(l_last_iter), *std::next(l_last_iter, 2), l, E, open);
                }
                else {
                    for (auto l_iter = l.begin(); l_iter != std::prev(l.end(), 2); ++l_iter) {
                         (*std::next(l_iter))->edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
                         (*l_iter)->edgein->next->twin->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
                         (*std::next(l_iter))->edgein->next->next = (*l_iter)->edgein->next->twin;
                         (*l_iter)->reg->firstEdge = (*l_iter)->edgein->next->twin;
                    }
                    l.clear();
                    E = 0;
                }
            } while (E != 0);
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
