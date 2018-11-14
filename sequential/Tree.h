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


/*
 * A structure of point : a point is defined by it's space
 * coordinates.
 */
struct point {

    double x;
    double y;

    point(double c1, double c2) : x(c1), y(c2) {};
};



/*
 * A structure to physically describe a box represented by a node
 * of the tree.
 */
struct box {

    struct point bottom_left;
    struct point bottom_right;
    struct point top_left;
    struct point top_right;


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
    point position;

    /*
     * A mass.
     */
    double mass;

    /*
     * A velocity.
     */
    point velocity;

    /*
     * A constructor.
     */
    particule(point p, double m, point v) : position(p), mass(m), velocity(v) {};
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
                            son2(0), son3(0), son4(0), part(point(0,0), 0., point(0.,0.)) {}

    /*
     * Constructor with a box and a point.
     */
    node(box b, particule p, double mass, bool father, int nbPart) : square(b), totalMass(mass),
         particulesNumber(nbPart), isEmpty(father), son1(0), son2(0), son3(0), son4(0), part(p) {}

};








class Tree {


public:

    /*
     * Empty constructor which build an empty tree.
     */
    Tree(box b, particule p);


    Tree(node n);

    /* **************************************************** */
    /*                 GETTERS & SETTERS                    */
    /* **************************************************** */

    node getNode(Tree tree) ;



    /* **************************************************** */
    /*                        METHODS                       */
    /* **************************************************** */

    /*
     * Print the features of a node.
     */
    static void printNode(node n);


    /*
     * Insert a node in the tree.
     */
    static void insertNode(node* n, particule p);


    /*
     * Given a node n, returns the subtree of n in which the particule p lies.
     */
    static node* getSubtree(point p, node n, int& child);

    /*
     * Iniatilize the four childs of a node n.
     */
    static node addChilds(node n);

    /*
     * Move the particule p which is in n into the good child.
     */
    static void move(particule p, node n);

    /*
     * Build a tree from a distribution of particules.
     */
    static void buildTree(node *root, std::vector<particule> distri);









    /*
     * Compute the distance between two points.
     */
    static double dist(point p, point q);

    /*
     * Compute the size of a box.
     */
    static double sizeBox(box b);

    /*
     * Compute the gravitationnal force that a (cluster of) particules exert on an other one
     */
    static point gForce(particule p, particule cm);

    /*
     * Move a particule under the force which it is subject to.
     */
    static particule moveParticule(particule p, point force, double t);

    /*
     * Compute the force on a particule p which is in the subtree n.
     */
    static point computeForce(particule p, node *tree, double theta);

    /*
     * Print a distribution of particule.
     */
    static void printDistri(std::vector<particule> distri);

//private:

    /*
     * Attribute of type node : a tree is a pointer on a node.
     */
    node *myNode;

    /*
     * Generate a random number between a and b.
     */
    static int randInt(int a, int b);

};

#endif //CODE_TREE_H
