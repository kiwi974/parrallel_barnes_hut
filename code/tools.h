//
// Created by kiwi974 on 17.05.18.
//


#include <math.h>
#include "Tree.h"
#include <vector>
#include <mpi.h>


#ifndef CODE_SORTTOOLS_H
#define CODE_SORTTOOLS_H


/*** Declaration of structures to represent distributions of particules ***/



/*
 * Structure distribution which contains the limits of the domain.
 */
struct distribution {
    int id;
    std::vector<particule> distri;
    struct box area;
    struct point cm;

    distribution(int id, std::vector<particule> d, box b, point cm) :
            id(id), distri(d), area(b), cm(cm) {}


    distribution() : id(0), distri({particule(0.,0.,0.,0.,0.)}), area(box()), cm(point(0,0)) {}
};



/*
 * Structure which is a list of distributions.
 */
struct distriList{
    struct distribution distri;
    struct distriList *next;

    distriList(distribution d) : distri(d), next(0) {}
};





/*** Methods to print distributions of particules ***/



/*
 * Print a distribution of particule.
 */
static void printDistri(std::vector<particule> distri);



/*
 * Print a structure distribution.
 */
void printDistriStruct(distribution distri);





/*** Merge sort methods for particules vectors : sort according to x or y, merge according to x or y ***/



/*
 * Split the distribution d in two and merge the two parts w.r.t the croissant order
 * of the abscissas.
 */
void fusionx(std::vector<particule> & T,int bg,int k,int bd);



/*
 * Sort the particule in the distribution according to the abscissa component.
 */
void fusionSortingx(std::vector<particule> & T,int bg =-1, int bd =-1);



/*
 * Split the distribution d in two and merge the two parts w.r.t the croissant order
 * of the ordinate.
 */
void fusiony(std::vector<particule> & T,int bg,int k,int bd);



/*
 * Sort the particule in the distribution according to the ordinate component.
 */
void fusionSortingy(std::vector<particule> & T,int bg =-1, int bd =-1);





/*** Methods usefull to split the distribution in the root process ***/



/*
 * Add a distribution to a list of distributions.
 */
distriList* addDistri(distriList *d, distribution dp);



/*
 * Reshape the box of a distribution after splitting : we cut the box in order to have to boxes
 * on one side and on the other of middle, according to a direction (1 for x, 0 for y).
 */
void reshapeBox(box d, box &d1, box &d2, double middle, int direction);



/*
 * Evenly split the distribution d among d1 and d2 (number of particules in d1 = number
 * of particules in d2 +- 1).
 * x and y coordinates of the center of mass of both distributions are computed in the same
 * time and stored into cm1 and cm2.
 */
void splitDistri(std::vector<particule> &d, std::vector<particule> &d1, std::vector<particule> &d2, double &middle,
                 int direction, point &cm1, point &cm2);



/*
 * Split a distribution according to a number of processors : this method call the previous one, as
 * many time as necessary to split the initial distribution in dlist in p subdivisions, all stored
 * in dlist.
 */
distriList* splitDistri(distriList *dlist, int p);




/*** Other methods used in order to perform the computation ***/



/*
 * Compute the longer side of a box.
 */
double dn(box b);



/*
 * Fill a vector paticules2Send with coordinates of particules that belongs to the Local
 * Essential Tree of proc otherproc.
 */
void findlocalLET(node nodeMyProc, distribution dOtherProc, std::vector<particule> &particules2Send, double threshold);



/*
 * Write the distribution d in a file named fileName (sequentil I/O).
 */
void writeDistri(char* fileName, std::vector<particule> d);



/* Write the disstribution d in a file named fileName with MPI */
void writeDistriMPI(char* fileName, std::vector<particule> d, int prank, int step, int nbTotalPart, int psize);



/*
 * Generate a random distribution of points.
 */
void generateDistri(int limDistri, std::vector<particule> &distri, int numbPart);



/*
 * Fill the vector toFill with data of type dtype got from process sender via communicator.
 */
int receiveDistriComm(int sender, MPI_Datatype dtype, MPI_Comm communicator, std::vector<particule> &toFill);



/*
 * Compute the new distribution after timestep dt.
 */
int moveParticules(std::vector<particule> &localDistri, double threshold, Tree localTree, double dt, int limSquare,
                   std::vector<particule> &partsOut, box localBox);



/*
 * Find the processor in which lie the particule p thanks to the routing table rootDistri.
 */
int findProcLie(std::vector<distribution> rootDistri, particule p);



#endif //CODE_SORTTOOLS_H
