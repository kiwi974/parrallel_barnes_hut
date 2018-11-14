//
// Created by kiwi974 on 14.05.18.
//


#include <iostream>
#include <typeinfo>
#include <vector>

#include "Tree.cpp"

int main() {

    /********* Building of the main tree *********/
    int limSquare = 1500;
    int limDistri = 10;

    point p1(0, 0);
    point p2(limSquare, 0);
    point p3(limSquare, limSquare);
    point p4(0, limSquare);

    box b(p1, p2, p3, p4);

    /********** Generate random particules *********/


    int numbPart(5);

    std::vector<particule> distribution;

    double bornSup = pow(10,4);
    double beta = (limDistri-bornSup)/(1-bornSup);
    double alpha = 1 - beta;

    srand(time(NULL));

    for (int i(0) ; i < numbPart ; i++) {
        /* Generating a random initial position */
        double x = Tree::randInt(1,bornSup);
        double y = Tree::randInt(1,bornSup);
        point pos(limSquare/4. + alpha*x + beta, limSquare/4. + alpha*y + beta);
        double m = Tree::randInt(1,9);
        /* Generating a random initial velocity */
        double vx = Tree::randInt(1,bornSup)/(bornSup);
        double vy = Tree::randInt(1,bornSup)/(bornSup);
        double vxThreshold = rand()/float(RAND_MAX);
        double vyThreshold = rand()/float(RAND_MAX);
        if (vxThreshold > 0.5) {
            vx = -vx;
        }
        if (vyThreshold > 0.5) {
            vy = -vy;
        }
        particule p(pos, m, point(vx,vy));
        distribution.push_back(p);
    }

    std::cout << "  " << std::endl;
    std::cout << " Initiale distribution. " << std::endl;
    Tree::printDistri(distribution);

    /********* Barnes-Hut *********/

    double dt = pow(10,-6);
    double t = dt;
    double theta = 0.2;
    int iteration = 1;

    std::vector<particule> *distridt;

    while (t <= 10*dt) {

        if (iteration == 1) {
            distridt = &distribution;
        }

        /* Building of the matching tree */

        Tree tree(b, particule(p1, 0., point(0., 0.)));
        Tree::buildTree(tree.myNode, *distridt);

        distridt->clear();

        /* Computing of the new positions after time dt */
        for (int i(0); i < numbPart; i++) {
            particule p = distribution[i];
            point f = Tree::computeForce(p, tree.myNode, theta);
            //std::cout << "force : (" << f.x << " , " << f.y << ")" << std::endl;
            particule movedParticule = Tree::moveParticule(p, f, dt);
            bool lim1 = movedParticule.position.x >= 0;
            bool lim2 = movedParticule.position.y >= 0;
            bool lim3 = movedParticule.position.x <= limSquare;
            bool lim4 = movedParticule.position.y <= limSquare;
            bool lim = lim1 && lim2 && lim3 && lim4;
            if (lim) {
                (*distridt).push_back(movedParticule);
            }
        }

        //std::cout << "Next distribution is : " << std::endl;
        //Tree::printDistri(*distridt);

        free(tree.myNode);

        t += dt;
        iteration += 1;

        //std::cout << "iteration : " << iteration << std::endl;

        //std::cin.get();
    }

    std::cout << "  " << std::endl;
    std::cout << "Final distribution." << std::endl;
    Tree::printDistri(*distridt);





    return 0;
}
