#pragma once

#ifndef NewDiagram_h
#define NewDiagram_h
 // STL
#include <vector>
#include <list>
#include <unordered_set>
// My includes
#include "Point2D.h"

namespace Voronoi{

    class NewDiagram{
    public:
        struct HalfEdge;
        struct Face;

        /**
         * \brief Point associated with a face of the partitioning:
         * it contains the index of the site, the coordinates and the capacity
         */
        struct Site
        {
            std::size_t index; /**< Index of the site in mSites */
            Point2D point; /**< Coordinates of the site */
            const double capacity; /**< Capacity of the site */
        };

        /**
         * \brief Vertex of a face: 
         * it contains the coordinates and the triplet
         */
        struct Vertex
        {
            Point2D point; /**< Coordinates of the vertex */
            std::vector<Site*> triplet; /**< Triplet of sites of which the point is a circumcenter */
            bool infinite = false; /*< If it's true the vertex is at infinity and the point value isn't significant **/

            static Vertex& getNullVertex() {
                static Vertex null_vertex; // Static instance
                null_vertex.infinite = true;
                return null_vertex;
            }

        private:
            typename std::list<Vertex>::iterator it;
        };

        

        /**
         * \brief Half-edge of a face:
         * it contains the tail, the head, the twin, the region, the label and next
         */
        struct HalfEdge
        {
            Vertex* tail = nullptr; /**< Origin vertex of the half-edge */
            Vertex* head = nullptr; /**< Destination vertex of the half-edge */
            HalfEdge* twin = nullptr; /**< Twin half-edge */
            Face* region; /**< Face to which this half-edge belongs to */
            Site* label = nullptr; /**< Site on the right side of the half-edge */
            HalfEdge* next = nullptr; /**< Next half-edge in the face frontier */

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
            std::list<Site*> sites; /**< List of sites associated with this face */
            Site* pivot = nullptr; /**< Pivot site associated with this face: it can be null */
            double weight = 0; /**< Sum of the weights of all sites associated with this face */
            HalfEdge* firstEdge; /**< A half-edge of the face */
            long ID; /**< Identifier of this face; it's unique */
            bool flag; /**< It's false if the face is to be eliminated */
        private:
            typename std::list<Face>::iterator it;
        };

        /**
         * \brief Get sites
         *
         * \return Const reference to the vector of sites of the diagram
         */
        const std::vector<Site>& getSites() const
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
        const Site* getSite(std::size_t i) const
        {
            return &mSites[i];
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
            return mSites[i].capacity;
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
         * \return Const reference to the list of faces of the diagram
         */
        std::list<Face>& getFaces()
        {
            return mFaces;
        }

        /**
         * \brief Get a face
         *
         * \param i Index of the site associated with the requested face
         *
         * \return Const pointer to the requested face
         
        const Face* getFace(std::size_t i) const
        {
            return &mFaces[i];
        }*/

        /**
         * \brief Get vertices
         *
         * \return Const reference to the list of vertices of the diagram
         */
        std::list<Vertex> getVertices()
        {
            return mVertices;
        }

        /**
         * \brief Get half-edges
         *
         * \return Reference to the list of half-edges of the diagram
         */
         std::list<HalfEdge>& getHalfEdges()
        {
            return mHalfEdges;
        }

       

    private:
        std::vector<Site> mSites; /**< Sites of the diagram */
        std::list<Face> mFaces; /**< Faces of the diagram */
        std::list<Vertex> mVertices; /**< Vertices of the diagram */
        std::list<HalfEdge> mHalfEdges; /**< Half-edges of the diagram */
        double mTotal;

        // Diagram construction
        //!!! Qui prima era private
    public:
        NewDiagram(const std::vector<std::pair<Point2D,double>>& points, double total)
        {
            mSites.reserve(points.size());
            for (auto i = std::size_t(0); i < points.size(); ++i)
            {
                mSites.push_back(NewDiagram::Site{ i, points[i].first, points[i].second});
            }
            mTotal = total;
        }

        Site* getSite(std::size_t i)
        {
            return &mSites[i];
        }

        /*Face* getFace(std::size_t i)
        {
            return &mFaces[i];
        }*/

        Vertex* createVertex(Point2D point)
        {
            mVertices.emplace_back();
            mVertices.back().point = point;
            //mVertices.back().it = std::prev(mVertices.end());
            return &mVertices.back();
        }

        

        HalfEdge* createHalfEdge(Site* site)
        {
            mHalfEdges.emplace_back();
            mHalfEdges.back().label = site;
            //mHalfEdges.back().it = std::prev(mHalfEdges.end());
            // This doesn't work anymore for my code
            /*if (face->firstEdge == nullptr)
                face->firstEdge = &mHalfEdges.back();*/
            return &mHalfEdges.back();
        }

        Face* createFace(long id) {
            mFaces.emplace_back();
            mFaces.back().ID = id;
            return &mFaces.back();
        }

        
        /*void removeVertex(Vertex* vertex)
        {
            mVertices.erase(vertex->it);
        }

        void removeHalfEdge(HalfEdge* halfEdge)
        {
            mHalfEdges.erase(halfEdge->it);
        }*/
    };

}

#endif