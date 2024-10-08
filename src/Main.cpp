#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <MyGAL/FortuneAlgorithm.h>
//#include "NewDiagram.h"
#include "MinKnapsack.h"
#include <SFML/Graphics.hpp>


namespace vor = Voronoi;

long progressiveID = 0;

double readCapacity(const std::string& fileName, int step = 1) {
    std::ifstream in(fileName);
    double capacity = 0.0;
    if (in) {
        in >> capacity;
    }
    return capacity;
}

// Function to read the points and their weights from a given file
std::vector<std::pair<Point2D, double>> readPoints(const std::string& fileName, int step = 1) {

    std::vector<std::pair<Point2D, double>> points;
    std::ifstream in(fileName);

    double x, y, weight;
    int points_n = 0, i = 0;
    double capacity = 0.0;

    if (in) {
        in >> capacity;
        in >> points_n;
        for (int p_i = 0; p_i < points_n; ++p_i) {
            in >> x >> y >> weight;
            if (p_i == i) {
                points.push_back(std::make_pair(Point2D(x, y), weight));
                i += step;
            }
        }
    }
    return points;
}

// Function to draw a grid
void drawGrid(sf::RenderWindow& window, int gridSpacing) {
    sf::VertexArray gridLines(sf::Lines);

    sf::Vector2u windowSize = window.getSize();

    // Vertical lines
    for (int x = 0; x < windowSize.x; x += gridSpacing) {
        gridLines.append(sf::Vertex(sf::Vector2f(x, 0), sf::Color(200, 200, 200)));  // Light gray color
        gridLines.append(sf::Vertex(sf::Vector2f(x, windowSize.y), sf::Color(200, 200, 200)));
    }

    // Horizontal lines
    for (int y = 0; y < windowSize.y; y += gridSpacing) {
        gridLines.append(sf::Vertex(sf::Vector2f(0, y), sf::Color(200, 200, 200)));
        gridLines.append(sf::Vertex(sf::Vector2f(windowSize.x, y), sf::Color(200, 200, 200)));
    }

    window.draw(gridLines);
}

// Function to initialize edge points for SFML
void initEdgePointsVis(Voronoi::NewDiagram::HalfEdgePtr& h, sf::VertexArray& line,
    const std::vector<std::pair<Point2D,double>>& points) {
    if (!h->head->infinite && !h->twin->head->infinite) {
        // Both vertices defined, draw line between them
        line[0].position = sf::Vector2f(h->head->point.x, h->head->point.y);
        line[1].position = sf::Vector2f(h->twin->head->point.x , h->twin->head->point.y);
    }
    else if (!h->head->infinite) {
        // One vertex is defined, project the other
        line[0].position = sf::Vector2f(h->head->point.x, h->head->point.y);
        Point2D norm = (points[h->twin->label->index].first - points[h->label->index].first).normalized().getRotated90CCW();
        line[1].position = sf::Vector2f((line[0].position.x + norm.x * 1000), (line[0].position.y + norm.y * 1000));
    }
    else if (!h->twin->head->infinite) {
        // Only the twin's vertex is defined, project the other
        line[0].position = sf::Vector2f(h->twin->head->point.x , h->twin->head->point.y);
        Point2D norm = (points[h->label->index].first - points[h->twin->label->index].first).normalized().getRotated90CCW();
        line[1].position = sf::Vector2f((line[0].position.x + norm.x * 1000), (line[0].position.y + norm.y * 1000));
    }
    else {
        // No vertices defined, project both points
        Point2D p1 = points[h->twin->label->index].first, p2 = points[h->label->index].first;
        Point2D norm = (p1 - p2).normalized().getRotated90CCW();
        Point2D c = 0.5 * (p1 + p2);
        line[0].position = sf::Vector2f((c.x + norm.x * 1000), (c.y + norm.y * 1000));
        line[1].position = sf::Vector2f((c.x - norm.x * 1000), (c.y - norm.y * 1000));
    }


    // Set the color of the line
    line[0].color = sf::Color::Green;
    line[1].color = sf::Color::Green;
}

// Function to modify the old Voronoi diagram structure to a new one
void modifyStructure(const mygal::Diagram<double>& diagram,
    Voronoi::NewDiagram& newDiagram) {

    // Create a map from old vertices to new vertices
    std::unordered_map<mygal::Diagram<double>::Vertex*, Voronoi::NewDiagram::VertexPtr> mapVertices;
    // Create a map from old halfedges to new halfedges
    std::unordered_map<mygal::Diagram<double>::HalfEdge*, Voronoi::NewDiagram::HalfEdgePtr> mapHalfedges;
    //Create a map from old regions to new regions
    std::unordered_map<mygal::Diagram<double>::Face*, Voronoi::NewDiagram::FacePtr> mapFaces;

    //Connects halfedges & vertices & faces
    for (auto& e : diagram.getHalfEdges()) {
        // I refer to the halfedge based on the pointer of the twin
        if (mapHalfedges.find(e.twin->twin) == mapHalfedges.end()) {
            auto e_new = newDiagram.createHalfEdge(newDiagram.getSite(e.incidentFace->site->index));
            mapHalfedges[e.twin->twin] = e_new;
        }
        

        if (mapHalfedges.find(e.twin) == mapHalfedges.end()) {
            auto e_twin = newDiagram.createHalfEdge(newDiagram.getSite(e.twin->incidentFace->site->index));
            mapHalfedges[e.twin] = e_twin;
        }

        auto e_new = mapHalfedges.at(e.twin->twin);
        e_new->twin = mapHalfedges.at(e.twin); 
        e_new->twin->twin = mapHalfedges.at(e.twin->twin);

        if (e_new->twin->twin != e_new) {
            std::cerr << "Twin consistency issue!" << std::endl;
        }

        if (mapFaces.find(e.incidentFace) == mapFaces.end()) {
            auto f_incident = newDiagram.createFace(progressiveID);
            f_incident->flag = true;
            f_incident->sites = { newDiagram.getSite(e.incidentFace->site->index) };
            f_incident->weight = newDiagram.getWeight(e.incidentFace->site->index);
            if (f_incident->weight <= newDiagram.getTotal()) {
                f_incident->pivot = nullptr;
            }
            else {
                f_incident->pivot = newDiagram.getSite(e.incidentFace->site->index);
            }
            mapFaces[e.incidentFace] = f_incident;
            progressiveID++;
        }
        e_new->region = mapFaces.at(e.incidentFace);

        // In the original structure the halfedges are connected in anti-clockwise sense,
        // but we need a clockwise sense so we assign the prev halfedge as next
        // and we assign the origin vertex as head and the destination vertex as tail

        // Check if there is a previous halfedge of the old halfedge
        if (e.prev != nullptr) {
            if (mapHalfedges.find(e.prev) == mapHalfedges.end()) {
                auto e_prev_new = newDiagram.createHalfEdge(newDiagram.getSite(e.prev->incidentFace->site->index));
                mapHalfedges[e.prev] = e_prev_new;
            }
            e_new->next = mapHalfedges.at(e.prev);
        }

        // Check if there is a origin vertex of the old halfedge
        Voronoi::NewDiagram::VertexPtr v1_new;
        if (e.origin != nullptr) {
            if (mapVertices.find(e.origin)==mapVertices.end()) {
                auto v_origin = newDiagram.createVertex(Point2D(e.origin->point.x, e.origin->point.y));
                v_origin->infinite = false;
                mapVertices[e.origin] = v_origin;
            }
            v1_new = mapVertices.at(e.origin);
        }
        else if (e_new->twin->tail != nullptr) {
            v1_new = e_new->twin->tail;
        }
        else {
            //v1_new = newDiagram.createVertex(NULL);
            //v1_new->infinite = true;
            v1_new = Voronoi::NewDiagram::Vertex::getNullVertex();
        }
        e_new->head = v1_new;


        // Analogous operations with the destination vertex of the old halfedge
        Voronoi::NewDiagram::VertexPtr v2_new;
        if (e.destination != nullptr) {
            if (mapVertices.find(e.destination) == mapVertices.end()) {
                auto v_destination = newDiagram.createVertex(Point2D(e.destination->point.x, e.destination->point.y));
                v_destination->infinite = false;
                mapVertices[e.destination] = v_destination;
            }
            v2_new = mapVertices.at(e.destination);
        }
        else if (e_new->twin->head != nullptr) {
            v2_new = e_new->twin->head;
        }
        else {
            //v2_new = newDiagram.createVertex(NULL);
            //v2_new->infinite = true;
            v2_new = Voronoi::NewDiagram::Vertex::getNullVertex();
        }
        e_new->tail = v2_new;




        // Assign the halfedge of the associated region if there aren't any already
        // or if this the first halfedge
        if ((e_new->region->firstEdge == nullptr) || (e_new->tail->infinite == true)) {
            e_new->region->firstEdge = e_new;
        }

        
    }
    //!!! Per ora la tripletta dei vertici la lascio vuota,
    // se serve la riempio sfruttando i semilati

  

    // Iterate again over all the unbounded faces to connect their last halfedge with their first halfedge,
    // this way the halfedges of each region form a circular list
    for (auto& f : newDiagram.getFaces()) {
        if (f->firstEdge->tail->infinite == true) {
            auto e = f->firstEdge;
            while (e->head->infinite != true) {
                e = e->next;
            }
            e->next = f->firstEdge;
        }
    }

    // Iterate again over all the halfedges to assign the right triplet to all vertices,
    // it's possible only now because we needed to connect all halfedges to each other
    for (auto& e : newDiagram.getHalfEdges()) {
        auto v1_new = e->head;
        // If the vertex to which the new halfedge points has not a triplet yet,
        // we create it based on the old halfedges that had it as their origin
        if (!v1_new->infinite && v1_new->triplet.empty()) {
            std::vector<Voronoi::NewDiagram::SitePtr> temp = std::vector<Voronoi::NewDiagram::SitePtr>(3);
            temp[0] = e->label;
            temp[1] = e->next->twin->label;
            temp[2] = e->next->twin->next->twin->label;
            v1_new->triplet = temp;
        }
    }

}



  int main()
  {
     //double total = 4;

    //std::string fileName = "Data/Lee_2.txt";
    //std::string fileName = "Data/Lee_1.txt"; <- WARNING: maybe there's cocircularity
    std::string fileName = "Data/example_3.txt";
    //std::string fileName = "Data/triangle.txt";
    //std::string fileName = "Data/equilatero_norm.txt";
    //std::string fileName = "Data/example_2.txt";
    //std::string fileName = "Data/example_1.txt";
    //std::string fileName = "Data/degenerate_1.txt";

    // Read points
    std::vector<std::pair<Point2D, double>> points_with_weights = readPoints(fileName);
    double capacity = readCapacity(fileName);

    std::vector<mygal::Vector2<double>> points;
    for (size_t i = 0; i < points_with_weights.size(); ++i) {
        points.push_back(mygal::Vector2<double>(points_with_weights[i].first.x, points_with_weights[i].first.y));
    }
 
    // Initialize an instance of Fortune's algorithm
    auto algorithm = mygal::FortuneAlgorithm<double>(points);
    // Construct the first order Voronoi diagram
    algorithm.construct();
    // Bound the diagram
    //algorithm.bound(Box<double>{-0.05, -0.05, 1.05, 1.05});
    // Get the constructed diagram
    auto diagram = algorithm.getDiagram();
    // Compute the intersection between the diagram and a box
    //diagram.intersect(Box<double>{0.0, 0.0, 1.0, 1.0});

     
    //Create the diagram with the sites and weights
    auto newDiagram = Voronoi::NewDiagram(points_with_weights, capacity);

    // Modify the diagram in order to get the structure I used in the pseudocode
    modifyStructure(diagram, newDiagram);

    /*auto& he = newDiagram.getHalfEdges();
    for (const auto& x : he) {
        std::cout << "Head: ";
        if (x.head->infinite != true) {
            std::cout << x.head->point.x << " " << x.head->point.y << "\n";
        }
        else {
            std::cout << "infinite vertex\n";
        }
        std::cout << "Tail: ";
        if (x.tail->infinite != true) {
            std::cout << x.tail->point.x << " " << x.tail->point.y << "\n";
        }
        else {
            std::cout << "infinite vertex\n";
        }
    }*/



    //Construct the min-knapsack Voronoi diagram
    //build_minKnapsack(newDiagram, points_with_weights, capacity);


    // Visualize the diagram
    // Create a window
    sf::RenderWindow window(sf::VideoMode(800, 600), "Plot Points");

    // Create a view that maps from [0, 1] in both axes to the window's size
    sf::View view;
    view.setCenter(0.5f, 0.5f);   // Center the view at (0.5, 0.5)
    view.setSize(1.0f, -1.0f);    // Set the size to (1, -1) to invert the y-axis

    // Apply this view to the window
    window.setView(view);

    // Define the points to plot
    float radius = 0.01f;
    sf::Color pointColor = sf::Color::Red;
    std::vector<sf::CircleShape> shapes;
    for (const auto& point : points_with_weights) {
        sf::CircleShape shape(radius);
        shape.setFillColor(pointColor);
        shape.setOrigin(radius, radius);
        shape.setPosition(point.first.x, point.first.y);
        shapes.push_back(shape);
    }

    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Handle window resizing
        if (event.type == sf::Event::Resized) {
            // Update the view to match the new window size
            sf::View view(sf::FloatRect(0, 0, 1, 1));
            window.setView(view);
        }

        // Clear window with white background
        window.clear(sf::Color::White); 
        // Draw the grid
        drawGrid(window, 50);

        // Draw all given points
        for (const auto& shape : shapes) {
            window.draw(shape);
        }


        // Draw all edges
        auto& faces = newDiagram.getFaces();
        for (auto& face : faces) {
            auto& edge = face->firstEdge;
            do
            {
                sf::VertexArray line(sf::Lines, 2);  // Create a line with 2 vertices
                initEdgePointsVis(edge, line, points_with_weights);  // Initialize the edge points
                window.draw(line);  // Draw the line
                edge = edge->next;
            } while (edge!=face->firstEdge);
            
        }

        window.display();
    }

        
    return 0;
  }
 
  //old version
  /*Voronoi::NewDiagram* modifyStructure(const Diagram<double>& diagram,
     const std::vector<std::pair<Point2D, double>>* points_with_weights, double total) {

     std::unordered_map<Diagram<double>::Vertex*, NewDiagram::Vertex*> vertexMap;
     std::unordered_map<Diagram<double>::HalfEdge*, NewDiagram::Vertex*> halfEdgeMap;

     //Create the diagram with the sites and weights
     auto newDiagram = Voronoi::NewDiagram(*points_with_weights, total);

     // Create a map from old halfedges to new halfedges
     std::unordered_map<mygal::Diagram<double>::HalfEdge*, Voronoi::NewDiagram::HalfEdge*> mapHalfedges;
     for (auto& e : diagram.getHalfEdges()) {
         auto e_new = newDiagram.createHalfEdge(newDiagram.getSite(e.incidentFace->site->index));
         mapHalfedges[&e] = e_new;
     }

     // Create a map from old vertices to new vertices
     std::unordered_map<mygal::Diagram<double>::Vertex*, Voronoi::NewDiagram::Vertex*> mapVertices;
     for (auto& v : diagram.getVertices()) {
         auto v_new = newDiagram.createVertex(Point2D(v.point.x, v.point.y));
         v_new->infinite = false;
         mapVertices[&v] = v_new;
     }


     //Create a map from old regions to new regions
     std::unordered_map<Diagram<double>::Face*, NewDiagram::Face*> mapFaces;
     for (auto& f : diagram.getFaces()) {
         auto f_new = newDiagram.createFace(progressiveID);
         f_new->flag = true;
         f_new->sites = { newDiagram.getSite(f.site->index) };
         f_new->weight = { newDiagram.getWeight(f.site->index) };
         if (f_new->weight >= newDiagram.getTotal()) {
             f_new->pivot = nullptr;
         }
         else {
             f_new->pivot = newDiagram.getSite(f.site->index);
         }
         mapFaces[&f] = f_new;
         progressiveID++;
     }


     //Connects halfedges & vertices & faces
     for (auto& e : diagram.getHalfEdges()) {
         auto& e_new = mapHalfedges.at(&e);


         e_new->twin = mapHalfedges.at(e.twin); //twin: 0x00000214b2d75850

         e_new->region = mapFaces.at(e.incidentFace);

         // In the original structure the halfedges are connected in anti-clockwise sense,
         // but we need a clockwise sense so we assign the prev halfedge as next
         // and we assign the origin vertex as head and the destination vertex as tail

         // Check if there is a previous halfedge of the old halfedge
         if (e.prev != nullptr) {
             e_new->next = mapHalfedges.at(e.prev);
         }

         // Check if there is a origin vertex of the old halfedge
         NewDiagram::Vertex* v1_new;
         if (e.origin != nullptr) {
             v1_new = mapVertices.at(e.origin);
         }
         else if (e_new->twin->tail != nullptr) {
             v1_new = e_new->twin->tail;
         }
         else {
             v1_new = newDiagram.createVertex(NULL);
             v1_new->infinite = true;
         }
         e_new->head = v1_new;


         // Analogous operations with the destination vertex of the old halfedge
         NewDiagram::Vertex* v2_new;
         if (e.destination != nullptr) {
             v2_new = mapVertices.at(e.destination);
         }
         else if (e_new->twin->head != nullptr) {
             v2_new = e_new->twin->head;
         }
         else {
             v2_new = newDiagram.createVertex(NULL);
             v2_new->infinite = true;
         }
         e_new->tail = v2_new;




         // Assign the halfedge of the associated region if there aren't any already
         // or if this the first halfedge
         if ((e_new->region->firstEdge == nullptr) || (e_new->tail->infinite == true)) {
             e_new->region->firstEdge = e_new;
         }
     }

     // Iterate again over all the unbounded faces to connect their last halfedge with their first halfedge,
     // this way the halfedges of each region form a circular list
     for (auto& f : newDiagram.getFaces()) {
         if (f.firstEdge->tail->infinite == true) {
             auto e = f.firstEdge;
             while (e->head->infinite != true) {
                 e = e->next;
             }
             e->next = f.firstEdge;
         }
     }

     // Iterate again over all the halfedges to assign the right triplet to all vertices,
     // it's possible only now because we needed to connect all halfedges to each other
     for (auto& e : newDiagram.getHalfEdges()) {
         auto v1_new = e.head;
         // If the vertex to which the new halfedge points has not a triplet yet,
         // we create it based on the old halfedges that had it as their origin
         if (v1_new->triplet.empty()) {
             std::vector<NewDiagram::Site*> temp = std::vector<NewDiagram::Site*>(3);
             temp[0] = e.label;
             temp[1] = e.next->twin->label;
             temp[2] = e.next->twin->next->twin->label;
             v1_new->triplet = temp;
         }
     }

     return &newDiagram;
 }*/
