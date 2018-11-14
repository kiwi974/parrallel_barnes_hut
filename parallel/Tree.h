//
// Created by kiwi974 on 01.05.18.
//

#include <math.h>
#include <vector>

#ifndef CODE_TREE_H
#define CODE_TREE_H


/*
 * The gravitational constant.
 */
const double G = 6.67408*pow(10,-11);





/*** Declaration of different structures ***/


/*
 * A structure of point : a point is defined by it's space
 * coordinates.
 */
struct point {

    double x;
    double y;

    point(double c1, double c2) : x(c1), y(c2) {};

    point() : x(0), y(0) {}
};



/*
 * A structure to physically describe a box represented by a node
 * of the tree.
 */
struct box {

    struct point bottom_left;
    struct point bottom_right;
    struct point top_right;
    struct point top_left;


    box() : bottom_left(point(0,0)), bottom_right(point(0,0)), top_right(point(0,0)), top_left(point(0,0)) {}

    box(point p1, point p2, point p3, point p4) :
            bottom_left(p1), bottom_right(p2), top_right(p3), top_left(p4) {}

};



/*
 * A structure to define a particule.
 */
struct particule {

    /* A
     * point with the coordinates
     */
    double x; double y;
    /*
     * A mass.
     */
    double mass;

    /*
     * A velocity.
     */
    double vx; double vy;

    /*
     * Constructors.
     */
    particule() : x(0.), y(0.), mass(0.), vx(0.), vy(0.) {};

    particule(double nx, double ny, double nm, double nvx, double nvy) : x(nx), y(ny), mass(nm), vx(nvx), vy(nvy) {};
};


/*
* A structure of node to build the tree, which is a pointer on
* a node. A node represents a box in the Barnes-Hut (Quad/Octo)tree.
*/
struct node {
    /* The physical box associated to the node. */
    struct box square;

    /* Mass of the set of particules in the node. */
    double totalMass;

    /* The total number of particules in the node. */
    int particulesNumber;

    /* Boolean which is true iff the node is not empty */
    bool isEmpty;

    /* The four sons of one node. */
    struct node *son1;
    struct node *son2;
    struct node *son3;
    struct node *son4;

    /* If the node is a leaf, it contains a at most a unique particule. Else, part
     * is the particule representing the center of mass of the cluster in the node.*/
    struct particule part;



    /*
     * Empty constructor.
     */
    node() : square(box()), totalMass(0), particulesNumber(0), isEmpty(true), son1(0),
                            son2(0), son3(0), son4(0), part(0., 0., 0., 0., 0.) {}

    /*
     * Constructor with a box and a point.
     */
    node(box b, particule p, double mass, bool father, int nbPart) : square(b), totalMass(mass),
         particulesNumber(nbPart), isEmpty(father), son1(0), son2(0), son3(0), son4(0), part(p) {}

};






/*** Beginning of the class Tree ***/

class Tree {


public:

    /* **************************************************** */
    /*                    CONSTRUCTORS                      */
    /* **************************************************** */

    /*
     * Build a tree from a box b and a particule p. All the other fields have
     * a default value (O, true...).
     */
    Tree(box b, particule p);


    /*
     * Build a tree from a node n. The tree will be a pointer on a node
     * with the same attributes than n.
     */
    Tree(node n);










    /* **************************************************** */
    /*                 GETTERS & SETTERS                    */
    /* **************************************************** */

    /*
     * Return the attribute myNode of a tree.
     */
    node* getNode() ;


    /*
     * Return the top box of the tree.
     */
    box getBox();







    /* **************************************************** */
    /*                        METHODS                       */
    /* **************************************************** */



    /*** Methods to build the quadtree associated to a distribution of particules. ***/



    /*
     * Print the features of a node.
     */
    static void printNode(node n);



    /*
     * Print a box.
     */
    static void printBox(box b);



    /*
     * Given a node n, returns the subtree of n in which the particule p lies.
     */
    static node* getSubtree(particule p, node n, int& child);



    /*
     * Iniatilize the four childs of a node n with default values : this way,
     * pointers are not empty and thus memory is allocated.
     */
    static node addChilds(node n);



    /*
     * Move the particule p which is in the node n into the good child.
     */
    static void move(particule p, node n);



    /*
     * Insert a particule in the tree from a node n.
     */
    static void insertNode(node* n, particule p);



    /*
     * Build a tree from a distribution of particules, by inserting particules one
     * after another from the root.
     */
    static void buildTree(node *root, std::vector<particule> distri);



    /*
     * Parallel version of buildTree : create a task per child of the root. Each task
     * is in charge to build it's own local tree, in the shared memory (thanks to OpenMP).
     */
    static void buildTreePara(node *root, std::vector<particule> &d3);









    /*** Methods to build the quadtree associated to a distribution of particules. ***/



    /*
     * Compute the distance between two points.
     */
    static double dist(point p, point q);



    /*
     * Compute the distance between two particules.
     */
    static double distPart(particule p, particule q);



    /*
     * Compute the size of a (squared) box.
     */
    static double sizeBox(box b);



    /*
     * Compute the gravitational force that a (cluster of) particules exert on an other one.
     */
    static point gForce(particule p, particule cm);



    /*
     * Move a particule under the force f which it is subject to, i.e. changes its coordinates
     * and its speed according to f.
     */
    static particule moveParticule(particule p, point force, double t);



    /*
     * Compute the force on a particule p which is in the subtree n, due to all the other
     * particules in the tree.
     */
    static point computeForce(particule p, node *tree, double theta);





private:

    /*
     * Attribute of type node : a tree is a pointer on a node.
     */
    node *myNode;

};

#endif //CODE_TREE_H