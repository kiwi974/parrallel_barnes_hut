#include <iostream>
#include <typeinfo>
#include <vector>

#include "Tree.cpp"

int main() {


    /* Building of the main tree */

    point p1(0, 0);
    point p2(15, 0);
    point p3(15, 15);
    point p4(0, 15);

    box b(p1, p2, p3, p4);

    Tree tree(b, particule(p1, 0., point(0., 0.)));

    point initVelocity(0., 0.);

    /* Declaration of the particules that we want to insert in that tree */

    particule part1(point(2, 5), 8.9, initVelocity);
    particule part2(point(12, 10), 3.8, initVelocity);
    particule part3(point(13, 13), 7.9, initVelocity);
    particule part4(point(3, 13), 10.8, initVelocity);
    particule part5(point(14, 1), 7.5, initVelocity);
    particule part6(point(12, 12), 8.9, initVelocity);
    particule part7(point(13, 14), 6.6, initVelocity);



    /* Building of the corresponding distribution */
    std::vector<particule> distribution;
    distribution.push_back(part1);
    distribution.push_back(part2);
    distribution.push_back(part3);
    distribution.push_back(part4);
    distribution.push_back(part5);
    distribution.push_back(part6);
    distribution.push_back(part7);


    Tree::printDistri(distribution);

    float dt = 0.01;
    float theta = 0.5;

    std::vector<particule> *distridt;

    while (dt <= 0.1) {

        if (dt == 0.01) {
            distridt = &distribution;
        }

        /* Building of the matching tree */

        Tree::buildTree(tree.myNode, *distridt);

        distridt->clear();

        /* Computing of the new positions after time dt */
        for (int i(0); i < 7; i++) {
            particule p = distribution[i];
            point f = Tree::computeForce(p, tree.myNode, theta);
            particule movedParticule = Tree::moveParticule(p, f, dt);
            (*distridt).push_back(movedParticule);
        }

        Tree::printDistri(*distridt);

        dt += dt;

        std::cin.get();
    }








    /* Printing of the non empty nodes of the tree */

    //std::cout << " " << std::endl;
    //std::cout << " --> Root tree" << std::endl;
    //Tree::printNode(*tree.myNode);
    //std::cout << " --> Its son 1." << std::endl;
    //Tree::printNode(*(tree.myNode->son1));
    //std::cout << " --> Its son 2." << std::endl;
    //Tree::printNode(*(tree.myNode->son2));
    //std::cout << " --> Its son 3.2." << std::endl;
    //Tree::printNode(*(tree.myNode->son3->son2));
    //std::cout << " --> Its son 3.3.1.3." << std::endl;
    //Tree::printNode(*(tree.myNode->son3->son3->son1->son3));
    //std::cout << " --> Its son 3.3.1.1." << std::endl;
    //Tree::printNode(*(tree.myNode->son3->son3->son1->son1));
    //std::cout << " --> Its son 3.3.4." << std::endl;
    //Tree::printNode(*(tree.myNode->son3->son3->son4));
    //std::cout << " --> Its son 4." << std::endl;
    //Tree::printNode(*(tree.myNode->son4));

    return 0;
}