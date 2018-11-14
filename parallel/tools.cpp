//
// Created by kiwi974 on 17.05.18.
//

/* File: trifusion.cc -- Francois -- Last modified on 11 Oct 2001
 *
 *      Algorithme de Tri Fusion - C++
 *
 */


#include "tools.h"
#include "Tree.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>


/*** Methods to print distributions of particules ***/



/*
 * Print a distribution of particule.
 */
void printDistri(std::vector<particule> distri) {
    std::cout << "****************************************************" << std::endl;
    std::cout << "size distri = " << distri.size() << std::endl;
    for (int i(0) ; i < distri.size() ; i++) {
        std::cout << "     (" << distri[i].x << ", " << distri[i].y << ", " << distri[i].mass
                  << ", " << distri[i].vx << ", " << distri[i].vy << ")" << std::endl;
    }
    std::cout << "****************************************************" << std::endl;
}

/*
 * Print a structure distribution.
 */
void printDistriStruct(distribution d) {
    std::vector<particule> distri(d.distri);
    point p1(d.area.bottom_left), p2(d.area.bottom_right), p3(d.area.top_right), p4(d.area.top_left);
    std::cout << "****************************************************" << std::endl;
    std::cout << "The distribution contains = " << distri.size() << " particules." << std::endl;
    std::cout << "The id of this distributions is : " << d.id << std::endl;
    std::cout << "Its center of mass is located in : (" << d.cm.x << " , " << d.cm.y << ")" << std::endl;
    std::cout << "It is in the box :[" << "(" << p1.x << "," << p1.y << ")" <<
              "; " << "(" << p2.x << "," << p2.y << ")" << "; " << "(" << p3.x << "," << p3.y << ")" <<
              "; " << "(" << p4.x << "," << p4.y << ")" << "]."<< std::endl;
    std::cout << "Its particules are located : " << std::endl;
    for (int i(0) ; i < distri.size() ; i++) {
        std::cout << "     (" << distri[i].x << ", " << distri[i].y << ")" << std::endl;
    }
    std::cout << "****************************************************" << std::endl;
}




/*** Merge sort methods for particules vectors : sort according to x or y, merge according to x or y ***/



/*
 * Split the distribution d in two and merge the two parts w.r.t the croissant order
 * of the abscissas.
 */
void fusionx(std::vector<particule> & d,int bg,int k,int bd) {
    std::vector<particule> temp;

    // Merging of the two subdistributions until the minium of their respective sizes
    int i1=bg, i2=k+1, i=0;
    while ((i1 <= k) && (i2 <= bd)) {
        if (d[i1].x < d[i2].x) {
            temp.push_back(d[i1]);
            i1++;
        }
        else {
            temp.push_back(d[i2]);
            i2++;
        }
        i++;
    }

    // Copy of the eventual queue if both subdistibution didn't had the same size
    if (i1 > k) {
        while (i2<=bd) {
            temp.push_back(d[i2]);
            i2++;
            i++;
        }
    }
    else {
        while (i1<=k) {
            temp.push_back(d[i1]);
            i1++;
            i++;
        }
    }

    // Copy of temp in the original vector
    for (i=bg;i<=bd;i++)
        d[i]=temp[i-bg];
}


/*
 * Sort the particule in the distribution according to the abscissa component.
 */

void fusionSortingx(std::vector<particule> & d,int bg, int bd) {
    int k;
    if (bg < bd) {
        k = (bg+bd)/2;
        fusionSortingx(d,bg,k);
        fusionSortingx(d,k+1,bd);
        fusionx(d,bg,k,bd);
    }
}



/*
 * Split the distribution d in two and merge the two parts w.r.t the croissant order
 * of the ordinate.
 */

void fusiony(std::vector<particule> & d,int bg,int k,int bd)
{
    std::vector<particule> temp;  //

    // Merging of the two subdistributions until the minium of their respective sizes
    int i1=bg, i2=k+1, i=0;
    while ((i1 <= k) && (i2 <= bd)) {
        if (d[i1].y < d[i2].y) {
            temp.push_back(d[i1]);
            i1++;
        }
        else {
            temp.push_back(d[i2]);
            i2++;
        }
        i++;
    }

    // Copy of the eventual queue if both subdistibution didn't had the same size
    if (i1 > k) {
        while (i2<=bd) {
            temp.push_back(d[i2]);
            i2++;
            i++;
        }
    }
    else {
        while (i1<=k) {
            temp.push_back(d[i1]);
            i1++;
            i++;
        }
    }

    // Copy of temp in the original vector
    for (i=bg;i<=bd;i++)
        d[i]=temp[i-bg];
}



/*
 * Sort the particule in the distribution according to the ordinate component.
 */

void fusionSortingy(std::vector<particule> & d,int bg, int bd) {
    int k;
    if (bg < bd) {
        k = (bg+bd)/2;
        fusionSortingy(d,bg,k);
        fusionSortingy(d,k+1,bd);
        fusiony(d,bg,k,bd);
    }
}





/*** Methods usefull to split the distribution in the root process ***/



/*
 * Add a distribution to a list of distributions.
 */
distriList* addDistri(distriList *d, distribution dp) {
    distriList *dtemp = new distriList(dp);
    dtemp->next = d;
    return dtemp;
}



/*
 * Reshape the box of a distribution after splitting : we cut the box in order to have to boxes
 * on one side and on the other of middle, according to a direction (1 for x, 0 for y).
 */
void reshapeBox(box d, box &d1, box &d2, double middle, int direction) {
    point bl(d.bottom_left), br(d.bottom_right), tr(d.top_right), tl(d.top_left);
    point bl1(0,0), br1(0,0), tr1(0,0), tl1(0,0);
    point bl2(0,0), br2(0,0), tr2(0,0), tl2(0,0);
    if (direction) {
        /* Coordinates of the left box */
        bl1 = bl;
        tl1 = tl;
        br1.x = middle; br1.y = br.y;
        tr1.x = middle; tr1.y = tr.y;
        /* Coordinates of the right box */
        bl2.x = middle; bl2.y = bl.y;
        tl2.x = middle; tl2.y = tl.y;
        br2 = br;
        tr2 = tr;
    } else {
        /* Coordinates of the under box */
        br1 = br;
        bl1 = bl;
        tl1.x = tl.x; tl1.y = middle;
        tr1.x = tr.x; tr1.y = middle;
        /* Coordinates of the above box */
        tr2 = tr;
        tl2 = tl;
        bl2.x = bl.x; bl2.y = middle;
        br2.x = br.x; br2.y = middle;
    }
    d1.bottom_left = bl1; d1.bottom_right = br1; d1.top_right = tr1; d1.top_left = tl1;
    d2.bottom_left = bl2; d2.bottom_right = br2; d2.top_right = tr2; d2.top_left = tl2;
}



/*
 * Evenly split the distribution d among d1 and d2 (number of particules in d1 = number
 * of particules in d2 +- 1).
 * x and y coordinates of the center of mass of both distributions are computed in the same
 * time and stored into cm1 and cm2.
 */
void splitDistri(std::vector<particule> &d, std::vector<particule> &d1, std::vector<particule> &d2,
                 double &middle, int direction, point &cm1, point &cm2) {
    int k(0), l(d.size()), l2(d.size() / 2);
    double tm1(0), tm2(0);
    while (k < l2) {
        cm1.x += d[k].mass*d[k].x;
        cm1.y += d[k].mass*d[k].y;
        tm1 += d[k].mass;
        d1.push_back(d[k]);
        k += 1;
    }
    cm1.x = cm1.x/tm1;
    cm1.y = cm1.y/tm1;

    if (direction) {
        middle = d[k-1].x + (d1[k].x - d[k-1].x) / 2.;
    } else {
        middle = d[k-1].y + (d1[k].y - d[k-1].y) / 2.;
    }

    while (k < l) {
        cm2.x += d[k].mass*d[k].x;
        cm2.y += d[k].mass*d[k].y;
        tm2 += d[k].mass;
        d2.push_back(d[k]);
        k += 1;
    }
    cm2.x = cm2.x/tm2;
    cm2.y = cm2.y/tm2;
}



/*
 * Split a distribution according to a number of processors : this method call the previous one, as
 * many time as necessary to split the initial distribution in dlist in p subdivisions, all stored
 * in dlist.
 */
distriList* splitDistri(distriList *dlist, int p) {

    int nbCut(p-1), direction(1);
    distriList *dtemp = NULL;
    distriList *dlistTemp = dlist;
    int id(1);

    while (nbCut != 0) {
        if (dlistTemp != NULL) {
            distribution currentDistri(dlistTemp->distri);
            std::vector<particule> d(currentDistri.distri);
            std::vector<particule> d1, d2;
            box b1, b2;
            point cm1(0,0), cm2(0,0);

            double middle(0);
            if (direction) {
                fusionSortingx(d, 0, d.size() - 1);
                splitDistri(d, d1, d2, middle, 1, cm1, cm2);
            } else {
                fusionSortingy(d, 0, d.size() - 1);
                splitDistri(d, d1, d2, middle, 0, cm1, cm2);
            }
            reshapeBox(currentDistri.area, b1, b2, middle, direction);
            distribution dist1(currentDistri.id,d1, b1, cm1), dist2(id,d2, b2, cm2);
            id += 1;
            dtemp = addDistri(dtemp, dist1);
            dtemp = addDistri(dtemp, dist2);
            dlistTemp = dlistTemp->next;
            nbCut -= 1;
            direction = (direction+1)%2;
        } else {
            dlistTemp = dtemp;
            dtemp = NULL;
        }
    }

    /* Merge of dtemp and dlist */
    dlist = NULL;
    while (dlistTemp != NULL) {
        dlist = addDistri(dlist,dlistTemp->distri);
        dlist = dlist->next;
    }
    while (dtemp != NULL) {
        dlist = addDistri(dlist,dtemp->distri);
        dtemp = dtemp->next;
    }
    return dlist;
}



/*** Other methods used in order to perform the computation ***/



/*
 * Compute the longer side of a box.
 */
double dn(box b) {
    double d1 = Tree::dist(b.bottom_left, b.bottom_right);
    double d2 = Tree::dist(b.bottom_left, b.top_left);
    double max(0);
    if (d1 < d2) {
        max = d1;
    } else {
        max = d2;
    }
    return max;
}



/*
 * Fill a vector paticules2Send with coordinates of particules that belongs to the Local
 * Essential Tree of proc otherproc.
 */
void findlocalLET(node nodeMyProc, distribution dOtherProc, std::vector<particule> &particules2Send, double threshold) {

    if (not(nodeMyProc.isEmpty)) {
        if (nodeMyProc.particulesNumber == 1) {
            particule p(nodeMyProc.part);
            particules2Send.push_back(p);
        } else {
            double dMyProc(dn(nodeMyProc.square));
            particule cmMyProc(nodeMyProc.part);
            point cmOtherProc(dOtherProc.cm);
            double distBetweenCm = Tree::dist(point(cmMyProc.x, cmMyProc.y), cmOtherProc);
            double ratio = dMyProc/distBetweenCm;
            if ((dMyProc/distBetweenCm) < threshold) {
                particule p(nodeMyProc.part);
                particules2Send.push_back(p);
            } else {
                // Recursive call on the four childs of nodeMyProc (given that
                // nodeMyProc is not empty, these childs exist by construction)
                findlocalLET(*nodeMyProc.son1,dOtherProc,particules2Send,threshold);
                findlocalLET(*nodeMyProc.son2,dOtherProc,particules2Send,threshold);
                findlocalLET(*nodeMyProc.son3,dOtherProc,particules2Send,threshold);
                findlocalLET(*nodeMyProc.son4,dOtherProc,particules2Send,threshold);
            }
        }
    }
}



/*
 * Write the distribution d in a file named fileName (sequentil I/O).
 */
void writeDistri(char* fileName, std::vector<particule> d){
    std::ofstream myfile;
    myfile.open(fileName,std::ios::out|std::ios::app);
    for (int i(0) ; i < d.size() ; i++) {
        if (i == d.size()-1) {
            myfile << d[i].x << "," << d[i].y << "<-->";
        } else {
            myfile << d[i].x << "," << d[i].y << "<>";
        }
    }
    myfile.close();
}



inline int digit(double x, int n) {
    int ndigit(std::trunc(x * std::pow(10., n)) - std::trunc(x * std::pow(10., n - 1)) *10.);
    if ((ndigit != 0) and (ndigit != 1) and (ndigit != 2) and (ndigit != 3) and (ndigit != 3) and (ndigit != 4)
        and (ndigit != 5) and (ndigit != 6) and (ndigit != 7) and (ndigit != 8) and (ndigit != 9)) {
        //std::cout << "ndigit = " << ndigit << std::endl;
        //std::cout << "converted = " << zero + ndigit << std::endl;
        ndigit = 0;
    }
    if ((ndigit != 0) and (ndigit != 1) and (ndigit != 2) and (ndigit != 3) and (ndigit != 3) and (ndigit != 4)
        and (ndigit != 5) and (ndigit != 6) and (ndigit != 7) and (ndigit != 8) and (ndigit != 9)) {
        std::cout << "ndigit = " << ndigit << std::endl;
        //std::cout << "converted = " << zero + ndigit << std::endl;
    }
    return ndigit;
}



/* Write the disstribution d in a file named fileName with MPI */
void writeDistriMPI(char* fileName, std::vector<particule> d, int prank, int step, int nbTotalPart, int psize) {

    int lenDistri(d.size());
    int lenDistriChar(17*lenDistri + 2*(lenDistri-1) + 2);
    int dstart(lenDistriChar*prank);
    int nbCharPreviousStep(lenDistriChar*(psize-1) + (lenDistriChar+2));
    int departure((step-1)*nbCharPreviousStep);

    /* Transform d in its string represenation */
    std::vector<char> toWrite;
    char zero = '0';
    for (int k(0) ; k < d.size() ; k++) {

        // Write coordinate x
        for (int i(0) ; i < 7 ; i++) {
            int ndigit(digit(d[k].x*pow(10,-2),i));
            toWrite.push_back(zero + ndigit);
            if (i==2){
                toWrite.push_back('.');
            }
        }

        toWrite.push_back(',');

        // Write coordinate y
        for (int i(0) ; i < 7 ; i++) {
            int ndigit(digit(d[k].y*pow(10,-2),i));
            toWrite.push_back(zero + ndigit);
            if (i==2){
                toWrite.push_back('.');
            }
        }

        if ((k == d.size()-1) and (prank == (psize-1))) {
            toWrite.push_back('<');
            toWrite.push_back('-');
            toWrite.push_back('-');
            toWrite.push_back('>');
        } else {
            toWrite.push_back('<');
            toWrite.push_back('>');
        }
    }

    /* open a file */
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    //MPI_File_set_size(file, 16);

    /* write the vector with MPI_File_write_at */
    MPI_File_write_at(file, departure + dstart, toWrite.data(), toWrite.size(), MPI_CHAR, MPI_STATUS_IGNORE);

    /* close the file */
    MPI_File_close(&file);

}


/*
 * Generate a random number between -1 and 1.
 */
double genNormal() {
    int limSup(pow(10,8));
    double gen((rand()%limSup)/(double)limSup);
    int threshold(rand()%2);
    if (threshold == 0) {
        gen = -gen;
    }
    return gen;
}



/*
 * Generate a random number between a and b.
 */
static int randInt(int a, int b);



/*
 * Generate a random boltzman distribution of points.
 */
void generateDistri(int limDistri, std::vector<particule> &distri, int numbPart) {
    double bornSup = pow(10, 4);
    double beta = (limDistri - bornSup) / (1 - bornSup);
    double alpha = 1 - beta;

    for (int i(0); i < numbPart; i++) {
        /* Generating a random initial position */
        double x = randInt(150, 850);
        double y = randInt(150, 850);
        double m = randInt(1, 9)*pow(10,-26);
        /* Generating a random initial velocity */
        double vx = genNormal()*pow(10,1); //sqrt(kb*T/m)*genNormal();
        double vy = genNormal()*pow(10,1); //sqrt(kb*T/m)*genNormal();
        particule p(x, y, m, vx, vy);
        distri.push_back(p);
    }
}



/*
 * Fill the vector toFill with data of type dtype got from process sender via communicator.
 */
int receiveDistriComm(int sender, MPI_Datatype dtype, MPI_Comm communicator, std::vector<particule> &toFill) {

    MPI_Status status;
    int nbPartReceived(0);

    // Probe for an incoming message from process otherproc */
    MPI_Probe(sender, 0, communicator, &status);

    // Get the message size when probe returns.
    MPI_Get_count(&status, dtype, &nbPartReceived);

    // Declare a buffer of length msgSize and receive the message.
    particule partReceive[nbPartReceived];
    MPI_Recv(&partReceive, nbPartReceived, dtype, sender, 0, communicator, MPI_STATUS_IGNORE);

    for (int i(0) ; i < nbPartReceived ; i++) {
        toFill.push_back(partReceive[i]);
    }

    return nbPartReceived;
}




/*
 * Compute the new distribution after timestep dt.
 */
int moveParticules(std::vector<particule> &localDistri, double threshold,
                   Tree localTree, double dt, int limSquare, std::vector<particule> &partsOut, box localBox) {

    int nbStayed(0);
    int n(localDistri.size());
    std::vector<particule> newLocalDistri;

    for (int i(0); i < n ; i++) {
        particule p = localDistri[i];
        point f = Tree::computeForce(p, localTree.getNode(), threshold);
        particule pMoved = Tree::moveParticule(p, f, dt);
        //std::cout << "pMoved = ( " << pMoved.x << ", " << pMoved.y << ")" << std::endl;
        bool lim1 = pMoved.x > 0.;
        bool lim2 = pMoved.y > 0.;
        bool lim3 = pMoved.x < limSquare;
        bool lim4 = pMoved.y < limSquare;
        bool lim = lim1 and lim2 and lim3 and lim4;
        if (lim) {
            /* Test wether moved particule is still in the right box or not */
            //Tree::printBox(localBox);
            bool limxinf = (localBox.bottom_left.x < pMoved.x);
            bool limxsup = (pMoved.x < localBox.bottom_right.x);
            bool limyinf = (localBox.bottom_left.y < pMoved.y);
            bool limysup = (pMoved.y < localBox.top_left.y);
            bool inLocalBox = limxinf and limxsup and limyinf and limysup;
            if (inLocalBox) {
                newLocalDistri.push_back(pMoved);
                nbStayed += 1;
            } else {
                partsOut.push_back(pMoved);
            }
        }
    }

    localDistri.clear();

    /* Root never enter that loop so it works */
    for (int i(0) ; i < nbStayed ; i++) {
        localDistri.push_back(newLocalDistri[i]);
    }

    return nbStayed;
}


/*
 * Find the processor in which lie the particule p thanks to the routing table
 */
int findProcLie(std::vector<distribution> rootDistri, particule p) {
    int proc(-1), n(rootDistri.size()), i(0);
    while ((i < n) and (proc == -1)) {
        box localBox = rootDistri[i].area;
        bool lim1 = localBox.bottom_left.x < p.x;
        bool lim2 = p.x < localBox.bottom_right.x;
        bool lim3 = localBox.bottom_left.y < p.y;
        bool lim4 = p.y < localBox.top_left.y;
        bool lim = lim1 and lim2 and lim3 and lim4;
        if (lim) {
            proc = i;
        }
        i += 1;
    }
    return proc;
}







/*
 * Generate a random number between a and b.
 */
int randInt(int a, int b){
    return rand()%(b-a) +a;
}









