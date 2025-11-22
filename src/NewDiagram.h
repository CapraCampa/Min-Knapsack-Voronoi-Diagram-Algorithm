#pragma once

#ifndef NewDiagram_h
#define NewDiagram_h
 // STL
#include <vector>
#include <memory>
#include <list>
#include <unordered_set>
// My includes
#include "Point2D.h"

namespace Voronoi{
    

    class NewDiagram{
    public:
        struct Face;
        struct HalfEdge;
        struct Vertex;
        struct Site;

        // I define the typename for the shared pointers
        using HalfEdgePtr = std::shared_ptr<HalfEdge> ;
        using VertexPtr = std::shared_ptr<Vertex> ;
        using FacePtr = std::shared_ptr<Face>;
        using SitePtr = std::shared_ptr<Site>;

        /**
         * \brief Point associated with a face of the partitioning:
         * it contains the index of the site, the coordinates and the capacity
         */
        struct Site
        {
            const std::size_t index; /**< Index of the site in mSites */
            const Point2D point; /**< Coordinates of the site */
            const double capacity; /**< Capacity of the site */

            Site(int i, Point2D p, double c)
                : index(i), point(p), capacity(c) {
            }

            friend std::ostream& operator<<(std::ostream& os, const Site& s) {
                os << s.index;
                return os;
            } 
        };

        /**
         * \brief Vertex of a face: 
         * it contains the coordinates and the triplet
         */
        struct Vertex
        {
            Point2D point; /**< Coordinates of the vertex */
            std::vector<SitePtr> triplet; /**< Triplet of sites of which the point is a circumcenter */
            bool infinite = false; /**< If it's true the vertex is at infinity and the point value isn't significant **/

            //Vertex(Point2D p) : point(p) {};

            // Returns the a vertex that is "at infinity", the value of the attribute point and triplet is not significant
            static VertexPtr getNullVertex() {
                static VertexPtr null_vertex = []() {
                    auto v = std::make_shared<Vertex>();
                    v->infinite = true;  // Set infinite once during initialization
                    return v;
                    }();
                return null_vertex;
            }


        private:
            //typename std::list<Vertex>::iterator it;
        };

        

        /**
         * \brief Half-edge of a face:
         * it contains the tail, the head, the twin, the region, the label and next
         */
        struct HalfEdge
        {
            VertexPtr tail = nullptr; /**< Origin vertex of the half-edge */
            VertexPtr head = nullptr; /**< Destination vertex of the half-edge */
            HalfEdgePtr twin = nullptr; /**< Twin half-edge */
            FacePtr region; /**< Face to which this half-edge belongs to */
            SitePtr label = nullptr; /**< Site on the right side of the half-edge */
            HalfEdgePtr next = nullptr; /**< Next half-edge in the face frontier */

            //HalfEdge(SitePtr s) : label(s) {};
        friend std::ostream& operator<<(std::ostream& os, const HalfEdge& e) {

            os << "Halfedge";

            // Safe region access
            if (e.region) {
                os << " belonging to region with ID: " << e.region->ID;
            } else {
                os << " [region: null]";
            }

            // Safe right site (label)
            if (e.label) {
                os << ", with site on the right: " << e.label->index;
            } else {
                os << ", with site on the right: null";
            }

            // Safe left site (through twin)
            if (e.twin && e.twin->label) {
                os << " and with site on the left: " << e.twin->label->index;
            } else {
                os << " and with site on the left: null";
            }

            os << "\n";

            // Safe tail
            if (e.tail) {
                os << ", tail: " << e.tail->point << "\n";
            } else {
                os << ", tail: null\n";
            }

            // Safe head
            if (e.head) {
                os << ", head: " << e.head->point << "\n";
            } else {
                os << ", head: null\n";
            }

            return os;
        }


        private:
            typename std::list<HalfEdge>::iterator it;
        };

        /**
         * \brief Structure representing a cell in the diagram
         *
         * The outer component of the face is represented as a circular doubly linked list.
         */
        struct Face
        {
            std::vector<SitePtr> sites; /**< Vector of sites associated with this face */
            SitePtr pivot = nullptr; /**< Pivot site associated with this face: it can be null */
            double weight = 0; /**< Sum of the weights of all sites associated with this face */
            HalfEdgePtr firstEdge = nullptr; /**< A half-edge of the face. If it open, it should be the first edge */
            long ID; /**< Identifier of this face; it's unique */
            bool flag = true; /**< It's false if the face is to be eliminated */

            friend std::ostream& operator<<(std::ostream& os, const Face& f) {
                os << "Region with ID: " << f.ID << ",weight: " << f.weight << ", sites:\n";
                for (auto& s : f.sites){
                    os << s->index << " ";
                }
                if (f.pivot != nullptr) {
                    os << ", pivot: " << f.pivot->index << "\n";
                }
                else {
                    os << ", pivot: null \n";
                }
                if(f.firstEdge!=nullptr){
                    os << ", firstedge:\n" << *f.firstEdge;
                }else{
                    os << ", firstedge null\n";
                }

                return os;
            } 

        private:
            typename std::list<Face>::iterator it;
        };

        /**
         * \brief Get sites
         *
         * \return Const reference to the vector of sites of the diagram
         */
        const std::vector<SitePtr>& getSites() const
        {
            return mSites;
        }

        /**
         * \brief Get a site
         *
         * \param i Index of the requested site
         *
         * \return Const pointer to the requested site
         */
        const SitePtr getSite(std::size_t i) const
        {
            return mSites[i];
        }

        /**
         * \brief Get the weight of a site
         *
         * \param i Index of the requested site
         *
         * \return Const weight of the requested site
         */
        const double getWeight(std::size_t i) const
        {
            return mSites[i]->capacity;
        }

        /**
         * \brief Get the total capacity to satisfy
         *
         *
         * \return Const total capacity to satisfy
         */
        const double getTotal() const
        {
            return mTotal;
        }

        /**
         * \brief Get the number of sites
         *
         * \return The number of sites
         */
        std::size_t getNbSites() const
        {
            return mSites.size();
        }

        /**
         * \brief Get faces
         *
         * \return Reference to the list of faces of the diagram
         */
        std::list<FacePtr>& getFaces()
        {
            return mFaces;
        }


        /**
         * \brief Get vertices
         *
         * \return Const reference to the list of vertices of the diagram
        std::list<VertexPtr>& getVertices()
        {
            return mVertices;
        }
         */

        /**
         * \brief Get half-edges
         *
         * \return Reference to the list of half-edges of the diagram
         std::list<HalfEdgePtr>& getHalfEdges()
        {
            return mHalfEdges;
        }
         */

       

    private:
        std::vector<SitePtr> mSites; /**< Sites of the diagram */
        std::list<FacePtr> mFaces; /**< Faces of the diagram */
        //std::list<VertexPtr> mVertices; /**< Vertices of the diagram */
        //std::list<HalfEdgePtr> mHalfEdges; /**< Half-edges of the diagram */
        const double mTotal;

        // Diagram construction
        //!!! Qui prima era private
    public:
        NewDiagram(const std::vector<std::pair<Point2D,double>>& points, double total)
            : mTotal(total) {
            mSites.reserve(points.size());
            for (auto i = std::size_t(0); i < points.size(); ++i)
            {
                mSites.push_back(std::make_shared<NewDiagram::Site>(i, points[i].first, points[i].second));
                 
            }
        }

        SitePtr getSite(std::size_t i)
        {
            return mSites[i];
        }

        /*Face* getFace(std::size_t i)
        {
            return &mFaces[i];
        }*/

        VertexPtr createVertex(Point2D point)
        {
            auto temp = std::make_shared<Vertex>();
            temp->point = point;
            return temp;
        }

        

        HalfEdgePtr createHalfEdge(SitePtr site)
        {
            HalfEdgePtr temp = std::make_shared<HalfEdge>();
            temp->label = site;
            return temp;

        }

        FacePtr createFace(long id) {
            mFaces.emplace_back(std::make_shared<Face>());
            mFaces.back()->ID = id;
            mFaces.back()->flag = 1;
            return mFaces.back();
        }

        void deleteHalfEdge(std::shared_ptr<HalfEdge>& he) {
            if (!he) return;

            // Break cycles
            he->twin.reset();
            he->next.reset();
            he->region.reset();
            he->tail.reset();
            he->head.reset();
            he->label.reset();

            he.reset(); // finally delete the shared_ptr itself
        }

        
        /*void removeVertex(Vertex* vertex)
        {
            mVertices.erase(vertex->it);
        }

        */
    };

    

}

#endif
