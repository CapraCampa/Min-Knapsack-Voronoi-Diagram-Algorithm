
#include "MinKnapsack.h"
#include <stdio.h>
#include "Point2D.h"

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
    b->head = f->tail;
    a->label = e->twin->label;
    b->label = f->twin->label;
    a->twin = b;
    b->twin = a;
    b->next = f;
    e->next = a;
}

//to implement
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

//to implement
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

int circumcenterInside(std::list<Special_vertex*>::iterator r_it,
    std::list<Special_vertex*>::iterator next_it,
    std::list<Special_vertex*>::iterator nextNext_it,
    std::list<Special_vertex*>& l, int E, bool open) {
    
    auto v_ptr = std::make_shared<Voronoi::NewDiagram::Vertex>();
    auto& v = *v_ptr;
    auto next = *next_it;
    auto r = *r_it;
    auto nextNext = *nextNext_it;
    v.triplet = std::vector{ next->eti_minus, next->eti_plus, nextNext->eti_plus };
    v.point = circumcenter(v.triplet);
    next->edgein->next->head = v_ptr;
    next->edgein->next->twin->tail = v_ptr;
    nextNext->edgein->next->head = v_ptr;
    nextNext->edgein->twin->tail = v_ptr;

    if (E == 3) {
        if (!open) {
            r->edgein->next->head = v_ptr;
            r->edgein->next->twin->tail = v_ptr;
        }
        else {
            createPendingEdge(next->edgein->next, nextNext->edgein->next->twin);
            r->edgein->next = next->edgein->next->next->twin;
            next->edgein->next->next->next = r->edgeout;
            r->edgein->next->tail = Voronoi::NewDiagram::Vertex::getNullVertex();
            r->edgein->next->twin->head = Voronoi::NewDiagram::Vertex::getNullVertex();
        }
        // Connect the three remaining pending edges
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

        // Remove the records next and nextNext from l and insert p in their place
        auto aux = std::next(nextNext_it);
        l.erase(next_it);
        l.erase(nextNext_it);

        auto p_it = l.insert(std::next(r_it), &p);

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
        Point2D P1 = circumcenter(p.eti_minus->point, p.eti_plus->point, (*std::next(p_it))->eti_minus->point);
        if (pointInsideRegion(P1, p.edgein->head->point, (*std::next(p_it))->edgeout->head->point)) {
            p.d_plus = dist(P, p.edgeout->tail->point);
            (*std::next(p_it))->d_minus = dist(P, (*std::next(p_it))->edgein->head->point);
        }
        else {
            p.d_plus = -dist(P, p.edgeout->tail->point);
            (*std::next(p_it))->d_minus = -dist(P, (*std::next(p_it))->edgein->head->point);
        }
        chooseDirection(r);
        chooseDirection(&p);
        chooseDirection(*std::next(p_it));
    }
    return E;
}


void partition(Voronoi::NewDiagram::FacePtr& r_ptr, std::list<Voronoi::NewDiagram::FacePtr>& R) {
    auto& r = *r_ptr;
    auto e_ptr = std::make_unique<std::vector<Voronoi::NewDiagram::HalfEdgePtr>>();
    auto& e = *e_ptr;
    auto lambda_ptr = std::make_unique<std::vector<Voronoi::NewDiagram::SitePtr>>();
    auto& lambda = *lambda_ptr;
    auto reg_ptr = std::make_unique<std::vector<Voronoi::NewDiagram::FacePtr>>();
    auto& reg = *reg_ptr;

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
            createPendingEdge(e.at(1),e.at(2));
            bool open = false;

            do {
                reg.push_back(createRegion(R, r_ptr, lambda.at(k+1), e.at(k)->next->twin));
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

            // Partition in 3 parts
            std::list<Special_vertex*> l{};
            int E = 0;
            for (int h = 0; h < k; h++) {
                Special_vertex rv{};
                rv.eti_minus = lambda.at(h);
                rv.eti_plus = lambda.at(h%(k+1));
                rv.reg = reg.at(h % (k + 1));
                rv.edgein = e.at(h);

                if (rv.edgein->head->infinite) {
                    rv.edgeout = e.at(h)->next;
                }
                else {
                    rv.edgeout = e.at(h)->next->twin->next;
                }
                l.push_back(&rv);
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

                    E = circumcenterInside(l_last_iter, std::next(l_last_iter), std::next(l_last_iter, 2), l, E, open);
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
        while (it != R.end() && (*it)->ID <= last) {
            std::cout << "I examine the region: " << (*it)->ID << "with pivot:"<< (*it)->pivot << "\n";
            if ((*it)->flag == 1 && (*it)->pivot == nullptr && ((*it)->weight < capacity)) {
                std::cout << "I divide the region: " << (*it)->ID << "\n";
                partition(*it, R);
                (*it)->flag = 0;
            }
            ++it;
        }
        //std::list<Voronoi::NewDiagram::Face>::iterator firstNewRegion = it;



        // Delete the element based on a condition
        /*if ((*it)->flag == 0) { // Example condition: delete even numbers
            it = R.erase(it); // erase returns the next iterator
        }
        else {
            ++it; // only increment if no deletion
        }*/
        
    }
}
