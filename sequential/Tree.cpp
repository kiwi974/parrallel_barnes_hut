//
// Created by kiwi974 on 01.05.18.
//

#include "Tree.h"
#include <iostream>
#include <cstdlib>


/* **************************************************** */
/*                    CONSTRUCTOR                       */
/* **************************************************** */

/*
 * Empty constructor which build an empty tree.
 */

Tree::Tree(box b, particule p) {
    myNode = new node(b, p,0,true,0);
}


Tree::Tree(node n) {
    myNode = new node(n.square,n.part,n.totalMass,n.isEmpty,n.particulesNumber);
}




/* **************************************************** */
/*                 GETTERS & SETTERS                    */
/* **************************************************** */




/* **************************************************** */
/*                        METHODS                       */
/* **************************************************** */

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
            particule p(point(0,0),0.,point(0.,0.));
            p = n.part;
            std::cout << "Its particule is located at : " << "(" << p.position.x << ", " << p.position.y << ")" << std::endl;
            std::cout << "Its speed is defined by the vector : " << "(" << n.part.velocity.x << ", " << n.part.velocity.y << ")" << std::endl;
        } else {
            std::cout << "This node is not a leaf, and has : " << n.particulesNumber<< " particules." << std::endl;
            std::cout << "Its center of mass is located at : " << "(" << n.part.position.x << ", " << n.part.position.y << ")" << std::endl;
            std::cout << "Its speed is defined by the vector : " << "(" << n.part.velocity.x << ", " << n.part.velocity.y << ")" << std::endl;
        }
    }
    std::cout << "*******************************************" << std::endl;
    std::cout << " " << std::endl;
}


/*
 * Generate a random number between a and b.
 */
int Tree::randInt(int a, int b){
    return rand()%(b-a) +a;
}



/*
 * Given a node n, returns the subtree of n in which the particule p lies.
 */

node* Tree::getSubtree(point p, node n, int& child) {

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
 * Iniatilize the four childs of a node n.
 */
node Tree::addChilds(node n) {

    // Creation of each node with he good square

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

    particule part0(point(0,0),0.,point(0.,0.));

    Tree childyN(n);

    childyN.myNode->son1 = new node(b1,part0,0,true,0);
    childyN.myNode->son2 = new node(b2,part0,0,true,0);
    childyN.myNode->son3 = new node(b3,part0,0,true,0);
    childyN.myNode->son4 = new node(b4,part0,0,true,0);

    return (*childyN.myNode);

}





/*
 * Move the particule p which is in n into the good child.
 */
void Tree::move(particule p, node n) {
    int child;
    node *c = getSubtree(p.position,n,child);
    insertNode(c,p);
}



/*
* Insert a node in the tree.
*/

void Tree::insertNode(node* n, particule p) {
    if ((*n).particulesNumber > 1) {
        //std::cout << "-----Stricly more than 1 particule in the current node.-----" << std::endl;
        //Modificaton of the center of mass
        double mass = (*n).totalMass;
        double xCmOld = (*n).part.position.x;
        double yCmOld = (*n).part.position.y;
        double xCm = (mass*xCmOld + p.mass*p.position.x)/(mass+p.mass);
        double yCm = (mass*yCmOld + p.mass*p.position.y)/(mass+p.mass);
        (*n).part.position.x = xCm;
        (*n).part.position.y = yCm;
        //--------------------------------
        int child;
        node *c = getSubtree(p.position,*n,child);
        (*n).particulesNumber += 1;
        (*n).totalMass += p.mass;
        insertNode(c,p);
    } else if ((*n).particulesNumber == 1) {
        //std::cout << "-----Exactly 1 particule in the current node.-----" << std::endl;
        if ((p.position. x != n->part.position.x) && (p.position.y != n->part.position.y)) {
            *n = addChilds(*n);
            move((*n).part, *n);
            //Modificaton of the center of mass
            double mass = (*n).totalMass;
            double xCmOld = (*n).part.position.x;
            double yCmOld = (*n).part.position.y;
            double xCm = (mass * xCmOld + p.mass * p.position.x) / (mass + p.mass);
            double yCm = (mass * yCmOld + p.mass * p.position.y) / (mass + p.mass);
            (*n).part.position.x = xCm;
            (*n).part.position.y = yCm;
            //--------------------------------
            (*n).particulesNumber += 1;
            (*n).totalMass += p.mass;
            int child;
            node *c = getSubtree(p.position, *n, child);
            insertNode(c, p);
        }
    } else {
        //std::cout << "We insert particule : (" << p.position.x << " , " << p.position.y << ")" << std::endl;
        (*n).particulesNumber = 1;
        (*n).part = p;
        (*n).totalMass += p.mass;
        (*n).isEmpty = false;
    }
}




/*
 * Build a tree from a distribution of particules.
 */
void Tree::buildTree(node *root, std::vector<particule> distri) {

    for (int i(0) ; i < distri.size() ; i++) {
        insertNode(root,distri[i]);
    }
}









/*
 * Compute the gravitationnal force that a (cluster of) particules exert on an other one
 */
point Tree::gForce(particule p, particule pcm) {
    point force(0.,0.);
    if ((p.position.x != pcm.position.x) && (p.position.y != pcm.position.y)) {
        double xdiff = pcm.position.x - p.position.x;
        double ydiff = pcm.position.y - p.position.y;
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
 * Compute the distance between two points.
 */
double Tree::dist(point p, point q) {
    double xdiff = p.x - p.y;
    double ydiff = p.y - q.y;
    double d = sqrt(pow(xdiff,2) + pow(ydiff,2));
    return d;
}

/*
 * Compute the size of a box.
 */
double Tree::sizeBox(box b) {
    return (abs(b.bottom_left.x - b.bottom_right.x));
}


/*
 * Move a particule under the force which it is subject to.
 */
particule Tree::moveParticule(particule p, point force, double t) {
    double fx = force.x;
    double fy = force.y;
    double vxt = p.velocity.x;
    double vyt = p.velocity.y;
    double xt = p.position.x;
    double yt = p.position.y;
    double m = p.mass;
    // Computation of the new position
    double xdt = xt + (fx*(pow(t,2)/2.) + vxt*t + xt)/m;
    double ydt = yt + (fy*(pow(t,2)/2.) + vyt*t + yt)/m;
    // Computation of the new velocity
    double vxdt = (fx*t + vxt)/m;
    double vydt = (fy*t + vyt)/m;
    //std::cout << " " << std::endl;
    //std::cout << "*****************************" << std::endl;
    //std::cout << "t = " << t << std::endl;
    //std::cout << "m = " << m << std::endl;
    //std::cout << "fx = " << fx << std::endl;
    //std::cout << "fy = " << fy << std::endl;
    //std::cout << "xt =" << xt << ", " << "xdt = " << xdt << std::endl;
    //std::cout << "yt =" << yt << ", " << "ydt = " << ydt << std::endl;
    //std::cout << "vxt =" << vxt << ", " << "vxdt = " << vxdt << std::endl;
    //std::cout << "vyt =" << vyt << ", " << "vydt = " << vydt << std::endl;
    //std::cout << "*****************************" << std::endl;
    //std::cout << " " << std::endl;
    // Return the corresponding new particule
    particule movedParticule(point(xdt,ydt),p.mass,point(vxdt,vydt));
    return movedParticule;
}



/*
 * Compute the force on a particule p which is in the subtree n.
 */

point Tree::computeForce(particule p, node *tree, double theta) {
    point f(0,0);
    if (tree->particulesNumber == 1) {
        //std::cout << "   Just one particule in that node" << std::endl;
        f = gForce(p,tree->part);
        return f;
    } else if (!tree->isEmpty){
        double r = dist(p.position,tree->part.position);
        double D = sizeBox(tree->square);
        if (D / r < theta) {
            //std::cout << "   We consider the action of the cluster of particules of this node." << std::endl;
            f = gForce(p,tree->part);
            return f;
        } else {
            //std::cout << " Can't appriximate this cluster by its center of mass." << std::endl;
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
            //std::cout << " " << std::endl;
            //std::cout << "*****************************" << std::endl;
            //std::cout << "f1 = " << f1.x << " , " << f1.y << std::endl;
            //std::cout << "f2 = " << f2.x << " , " << f2.y << std::endl;
            //std::cout << "f3 = " << f3.x << " , " << f3.y << std::endl;
            //std::cout << "f4 = " << f4.x << " , " << f4.y << std::endl;
            //std::cout << "*****************************" << std::endl;
            //std::cout << " " << std::endl;
            double fx = f1.x + f2.x + f3.x + f4.x;
            double fy = f1.y + f2.y + f3.y + f4.y;
            f = point(fx,fy);
            return f;
        }
    }
}

/*
 * Print a distribution of particule.
 */
void Tree::printDistri(std::vector<particule> distri) {
    std::cout << "****************************************************" << std::endl;
    for (int i(0) ; i < distri.size() ; i++) {
        std::cout << "     (" << distri[i].position.x << ", " << distri[i].position.y << ")" << std::endl;
    }
    std::cout << "****************************************************" << std::endl;
}



