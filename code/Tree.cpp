//
// Created by kiwi974 on 01.05.18.
//

#include "Tree.h"
#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "tools.h"


/* **************************************************** */
/*                    CONSTRUCTOR                       */
/* **************************************************** */

/*
 * Build a tree from a box b and a particule p. All the other fields have
 * a default value (O, true...).
 */

Tree::Tree(box b, particule p) {
    myNode = new node(b, p,0,true,0);
}



/*
 * Build a tree from a node n. The tree will be a pointer on a node
 * with the same attributes than n.
 */
Tree::Tree(node n) {
    myNode = new node(n.square,n.part,n.totalMass,n.isEmpty,n.particulesNumber);
}










/* **************************************************** */
/*                 GETTERS & SETTERS                    */
/* **************************************************** */

/*
 * Return the attribute myNode of a tree.
 */

node* Tree::getNode() {
    return myNode;
}


/*
     * Return the top box of the tree.
     */
box Tree::getBox() {
    return myNode->square;
}







/* **************************************************** */
/*                        METHODS                       */
/* **************************************************** */



/*** Methods to build the quadtree associated to a distribution of particules. ***/



/*
 * Print the features of a node.
 */

void Tree::printNode(node n) {
    std::cout << " " << std::endl;
    std::cout << "*******************************************" << std::endl;
    point p1(0,0);
    point p2(0,0);
    point p3(0,0);
    point p4(0,0);
    p1 = n.square.bottom_left;
    p2 = n.square.bottom_right;
    p3 = n.square.top_right;
    p4 = n.square.top_left;
    std::cout << "In the box : [" << "(" << p1.x << "," << p1.y << ")" <<
              "; " << "(" << p2.x << "," << p2.y << ")" << "; " << "(" << p3.x << "," << p3.y << ")" <<
              "; " << "(" << p4.x << "," << p4.y << ")" << "]."<< std::endl;
    std::cout << "The total mass of this node is : " << n.totalMass << "." << std::endl;
    if (n.isEmpty) {
        std::cout << "This node is an empty leaf." << std::endl;
    } else {
        if (n.particulesNumber == 1) {
            std::cout << "This node is a leaf and contains " << n.particulesNumber << " particule." << std::endl;
            particule p(0. ,0., 0., 0.,0.);
            p = n.part;
            std::cout << "Its particule is located at : " << "(" << p.x << ", " << p.y << ")" << std::endl;
            std::cout << "Its speed is defined by the vector : " << "(" << n.part.vx << ", " << n.part.vy << ")" << std::endl;
        } else {
            std::cout << "This node is not a leaf, and has : " << n.particulesNumber<< " particules." << std::endl;
            std::cout << "Its center of mass is located at : " << "(" << n.part.x << ", " << n.part.y << ")" << std::endl;
            std::cout << "Its speed is defined by the vector : " << "(" << n.part.x << ", " << n.part.y << ")" << std::endl;
        }
    }
    std::cout << "*******************************************" << std::endl;
    std::cout << " " << std::endl;
}




/*
     * Print a box.
     */
void Tree::printBox(box b) {
    point p1(0,0);
    point p2(0,0);
    point p3(0,0);
    point p4(0,0);
    p1 = b.bottom_left;
    p2 = b.bottom_right;
    p3 = b.top_right;
    p4 = b.top_left;
    std::cout << "Box : [" << "(" << p1.x << "," << p1.y << ")" <<
              "; " << "(" << p2.x << "," << p2.y << ")" << "; " << "(" << p3.x << "," << p3.y << ")" <<
              "; " << "(" << p4.x << "," << p4.y << ")" << "]."<< std::endl;
}




/*
 * Given a node n, returns the subtree of n in which the particule p lies.
 */

node* Tree::getSubtree(particule p, node n, int& child) {

    if ((n.son1->square.bottom_left.x <= p.x) && (p.x <= n.son1->square.bottom_right.x) &&
        (n.son1->square.bottom_left.y <= p.y) && (p.y <= n.son1->square.top_left.y)) {
        child = 1;
        return n.son1;
    } else if ((n.son2->square.bottom_left.x <= p.x) && (p.x <= n.son2->square.bottom_right.x) &&
               (n.son2->square.bottom_left.y <= p.y) && (p.y <= n.son2->square.top_left.y)) {
        child = 2;
        return n.son2;
    } else if ((n.son3->square.bottom_left.x <= p.x) && (p.x <= n.son3->square.bottom_right.x) &&
               (n.son3->square.bottom_left.y <= p.y) && (p.y <= n.son3->square.top_left.y)) {
        child = 3;
        return n.son3;
    } else {
        child = 4;
        return n.son4;
    }
}




/*
 * Iniatilize the four childs of a node n with default values : this way,
 * pointers are not empty and thus memory is allocated.
 */

node Tree::addChilds(node n) {

    // Creation of each node with the good square

    double x1(n.square.top_left.x);
    double x2(n.square.top_right.x);
    double midlex(x1 + (x2 - x1)/2.);

    double y1(n.square.bottom_left.y);
    double y2(n.square.top_left.y);
    double midley(y1 + (y2-y1)/2.);

    point p1(x1,midley);
    point p2(midlex,y2);
    point p3(x2,midley);
    point p4(midlex,y1);
    point p5(midlex,midley);

    box b1(n.square.bottom_left,p4,p5,p1);
    box b2(p4,n.square.bottom_right,p3,p5);
    box b3(p5,p3,n.square.top_right,p2);
    box b4(p1,p5,p2,n.square.top_left);

    // Creation of every son

    particule part0(0., 0., 0., 0., 0.);

    Tree childyN(n);

    childyN.myNode->son1 = new node(b1,part0,0,true,0);
    childyN.myNode->son2 = new node(b2,part0,0,true,0);
    childyN.myNode->son3 = new node(b3,part0,0,true,0);
    childyN.myNode->son4 = new node(b4,part0,0,true,0);

    return (*childyN.myNode);
}





/*
 * Move the particule p which is in the node n into the good child.
 */

void Tree::move(particule p, node n) {
    int child;
    node *c = getSubtree(p,n,child);
    insertNode(c,p);
}



/*
 * Add the contribution of a particule to a center of mass.
 */
void addContribution(particule p, particule &cm, double totalmass) {
    double xCmOld = cm.x;
    double yCmOld = cm.y;
    double vxCmOld = cm.vx;
    double vyCmOld = cm.vy;
    double xCm = (totalmass * xCmOld + p.mass * p.x) / (totalmass + p.mass);
    double yCm = (totalmass * yCmOld + p.mass * p.y) / (totalmass + p.mass);
    double vxCm = (totalmass * vxCmOld + p.mass * p.vx) / (totalmass + p.mass);
    double vyCm = (totalmass * vyCmOld + p.mass * p.vy) / (totalmass + p.mass);
    cm.x = xCm;
    cm.y = yCm;
    cm.vx = vxCm;
    cm.vy = vyCm;
}


/*
 * Insert a particule in the tree from a node n.
 */

void Tree::insertNode(node* n, particule p) {
    if ((*n).particulesNumber > 1) {
        //Modificaton of the center of mass
        addContribution(p,(*n).part,(*n).totalMass);
        (*n).particulesNumber += 1;
        (*n).totalMass += p.mass;
        //--------------------------------
        int child;
        node *c = getSubtree(p,*n,child);
        insertNode(c,p);
    } else if ((*n).particulesNumber == 1) {
        if ((p. x != n->part.x) && (p.y != n->part.y)) {
            *n = addChilds(*n);
            move((*n).part, *n);
            //Modificaton of the center of mass
            addContribution(p,(*n).part,(*n).totalMass);
            //--------------------------------
            (*n).particulesNumber += 1;
            (*n).totalMass += p.mass;
            int child;
            node *c = getSubtree(p, *n, child);
            insertNode(c, p);
        }
    } else {
        (*n).particulesNumber = 1;
        (*n).part = p;
        (*n).totalMass += p.mass;
        (*n).isEmpty = false;
    }
}



/*
 * Build a tree from a distribution of particules, by inserting particules one
 * after another from the root.
 */

void Tree::buildTree(node *root, std::vector<particule> distri) {
    for (int i(0) ; i < distri.size() ; i++) {
        insertNode(root,distri[i]);
    }
}



/*
 * Parallel version of buildTree : create a task per child of the root. Each task
 * is in charge to build it's own local tree, in the shared memory (thanks to OpenMP).
 */

void Tree::buildTreePara(node *root, std::vector<particule> &d) {


    /*** Evenly split the initial distribution in four distributions ***/

    std::vector<particule> dleft, dright;
    std::vector<particule> d1, d2, d3, d4;
    double middle;

    // Cut along x
    fusionSortingx(d,0,d.size()-1);
    point cm1(0,0), cm2(0,0);
    splitDistri(d, dleft,dright, middle, 1, cm1, cm2);

    // Cut along y
    fusionSortingy(dleft,0,dleft.size()-1);
    fusionSortingy(dright,0,dright.size()-1);

    splitDistri(dleft,d1,d4, middle, 0, cm1, cm2);
    splitDistri(dright,d2,d3, middle, 0, cm1, cm2);

    /*** Insertion in the tree ***/

    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task
            {
                for (int i(0) ; i < d1.size() ; i++) {
                    insertNode(root,d1[i]);
                }
            }
            #pragma omp task
            {
                for (int i(0) ; i < d2.size() ; i++) {
                    insertNode(root,d2[i]);
                }
            }
            #pragma omp task
            {
                for (int i(0) ; i < d3.size() ; i++) {
                    insertNode(root,d3[i]);
                }
            }
            #pragma omp task
            {
                for (int i(0) ; i < d4.size() ; i++) {
                    insertNode(root,d4[i]);
                }
            }
        }
    };
}





/*** Methods to build the quadtree associated to a distribution of particules. ***/



/*
 * Compute the distance between two points.
 */
double Tree::dist(point p, point q) {
    double xdiff = p.x - q.x;
    double ydiff = p.y - q.y;
    double d = sqrt(pow(xdiff,2) + pow(ydiff,2));
    return d;
}


/*
 * Compute the distance between two particules.
 */
double Tree::distPart(particule p, particule q) {
    double xdiff = p.x - q.x;
    double ydiff = p.y - q.y;
    double d = sqrt(pow(xdiff,2) + pow(ydiff,2));
    return d;
}



/*
 * Compute the size of a (squared) box.
 */
double Tree::sizeBox(box b) {
    return (abs(b.bottom_left.x - b.bottom_right.x));
}



/*
 * Compute the gravitational force that a (cluster of) particules exert on an other one
 */
point Tree::gForce(particule p, particule pcm) {
    point force(0.,0.);
    if ((p.x != pcm.x) && (p.y != pcm.y)) {
        double xdiff = pcm.x - p.x;
        double ydiff = pcm.y - p.y;
        double r = sqrt(pow(xdiff, 2) + pow(ydiff, 2));
        double G_m_mcm = G * p.mass * pcm.mass;
        double forcex = G_m_mcm * (xdiff) / (pow(r, 3));
        double forcey = G_m_mcm * (ydiff) / (pow(r, 3));
        force.x = forcex;
        force.y = forcey;
    }
    return force;
}







/*
 * Move a particule under the force f which it is subject to, i.e. changes its coordinates
 * and its speed according to f.
 */

particule Tree::moveParticule(particule p, point force, double dt) {
    double fx = force.x;
    double fy = force.y;
    double xt = p.x;
    double yt = p.y;
    double vxt = p.vx;
    double vyt = p.vy;
    double m = p.mass;
    // Computation of the new position
    double xdt = (fx/m)*(pow(dt,2)/2.) + vxt*dt + xt;
    double ydt = (fy/m)*(pow(dt,2)/2.) + vyt*dt + yt;
    // Computation of the new velocity
    double vxdt = (fx/m)*dt + vxt;
    double vydt = (fy/m)*dt + vyt;
    // Return the corresponding new particule
    particule movedParticule(xdt,ydt,p.mass,vxdt,vydt);
    return movedParticule;
}



/*
 * Compute the force on a particule p which is in the subtree n, due to all the other
 * particules in the tree.
 */

point Tree::computeForce(particule p, node *tree, double theta) {
    point f(0,0);
    if (tree->particulesNumber == 1) {
        f = gForce(p,tree->part);
    } else if (!tree->isEmpty){
        double r = distPart(p,tree->part);
        double D = sizeBox(tree->square);
        if (D / r < theta) {
            f = gForce(p,tree->part);
        } else {
            point f1(0,0), f2(0,0), f3(0,0), f4(0,0);
            if (&tree->son1 != 0) {
                f1 = computeForce(p, tree->son1, theta);
            }
            if (&tree->son2 != 0) {
                f2 = computeForce(p,tree->son2,theta);
            }
            if (&tree->son3 != 0) {
                f3 = computeForce(p,tree->son3,theta);
            }
            if (&tree->son4 != 0) {
                f4 = computeForce(p,tree->son4,theta);
            }
            double fx = f1.x + f2.x + f3.x + f4.x;
            double fy = f1.y + f2.y + f3.y + f4.y;
            f = point(fx,fy);
        }
    }
    return f;
}