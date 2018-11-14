#include <iostream>
#include <typeinfo>
#include <vector>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <stddef.h>
#include <sstream>


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "tools.cpp"
#include "Tree.cpp"


int main(int argc, char *argv[]) {

    /*** Declaration ***/

    // Total number of iterations performed
    int nbIt;

    // Number of particules we generate for the initial distribution
    int numbPart;


    std::istringstream iss1( argv[1] );
    std::istringstream iss2( argv[2] );
    if (iss1 >> nbIt) {}
    if (iss2 >> numbPart) {}


    // Number of procesors and rank of eahc one of th
    int psize, prank;  // psize = 2^m

    //theta : threshold to know whether we consider a center of mass to compute forces or not
    double threshold(0.5);

    // Timestamp and current time during the simulation
    float dt(0.01), t(0.01);

    // Limit of the simulation domain (0,limSquare)x(0,limSquare)
    int limSquare = 1500;

    // Side of the square in which all particules of the initial distribution lie
    int limDistri = 10;

    // Global box representing the computation area
    box b0(point(0,0), point(limSquare,0), point(limSquare,limSquare), point(0,limSquare));

    // Sunchronization of the clock with the module random.
    srand(time(NULL));

    // Name of the text file in which we will write the successive distributions.
    char* fileName("distriTest");

    // Know wether we need to build all again or not
    bool buildAll(t==dt);

    // Number of iteration we wait before to rebuild the whole tree
    int load = 10;

    // Average time pet step
    double averageTimeStep(0);

    // MPI initialization
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);


    //Vector of struct distribution which is the rooting table owned by each processor
    std::vector<distribution> rootDistri(psize);

    // Buffer to share the routing table among processors
    double boxes[11*psize];

    // Vector of particules that will contain all the particules
    std::vector<particule> dAllParticules;

    // Vector of particules reprensenting the local distribution of each processor
    std::vector<particule> localDistri;

    MPI_Request sendLET;
    MPI_Status status;



    /************************************* Decalration of MPIType particule *******************************************/


    MPI_Datatype mpiParticule;
    MPI_Type_contiguous(5,MPI_DOUBLE,&mpiParticule);
    MPI_Type_commit(&mpiParticule);

    /******************************************************************************************************************/



    /*** Global chrono ***/
    using clk = std::chrono::high_resolution_clock;
    using second = std::chrono::duration<double>;
    auto t1 = clk::now();
    auto tPerStep = clk::now();



    /*** Delete the file distriTest previously wrote */
    if (prank == 0) {
        if (remove(fileName) != 0) {
            perror("Error deleting file");
        } else {
            puts("File successfully deleted");
        }
    }


    while (t <= nbIt*dt) {

        /*** Start scaling chrono ***/
        if (prank ==0) {
            tPerStep = clk::now();
        }

        /********************* Generation of an intiale distribution or reconstruction of *****************************/
        /*********************************** the whole tree if its necessary ******************************************/
        /************************************************ or **********************************************************/
        /********************* Simply receive the new local distribution of moved particules **************************/
        /*********************************** without rebuilding all the tree ******************************************/

        if (buildAll) { /*** In that case we (re)build all the distribution to (re)build all the tree among procs ***/

            if (prank == 0) {

                if ((int(t/dt)%10)==0) {
                    //std::cout << " " << std::endl;
                    std::cout << "******************** STEP " << t / dt << " ********************" << std::endl;
                    //std::cout << "t = " << t << std::endl;
                }

                //std::cout << "Rebuilding the whole tree." << std::endl;

                if (t == dt) {
                    /* Generation of the initial distribution */
                    std::cout << "Beginning of the simulation for " << numbPart << " particules..." << std::endl;
                    generateDistri(limDistri, dAllParticules, numbPart);
                } else {
                    dAllParticules.clear();
                    /* Receive the moved particules from every process which is not the root */
                    int nbCoordPreviousIt(0);
                    for (int i(1); i < psize; i++) {
                        nbCoordPreviousIt = receiveDistriComm(i, mpiParticule, MPI_COMM_WORLD, dAllParticules);
                    }
                }

                for (int i(0) ; i < localDistri.size() ; i++) {
                    dAllParticules.push_back(localDistri[i]);
                }
                std::cout << "We deal with " << dAllParticules.size() << " particules" << std::endl;

                /* Initialization of the global distribution */
                distribution *distri = new distribution(0, dAllParticules, b0, point());

                /* Write the distribution in the designed text file to plot the simulation later. */
                writeDistri(fileName, dAllParticules);


                /* Split distribution to evenly deal particules with all processors. */
                distriList *distriLProc = new distriList(*distri);
                distriLProc = splitDistri(distriLProc, psize);

                /* Built the routing table. */
                int k(0);
                distriList *distriLProcTemp = distriLProc;
                while (distriLProcTemp != NULL) {
                    distribution d = distriLProcTemp->distri;
                    boxes[k] = d.id;
                    k += 1;
                    boxes[k] = d.cm.x;
                    boxes[k + 1] = d.cm.y;
                    k += 2;
                    boxes[k] = d.area.bottom_left.x;
                    boxes[k + 1] = d.area.bottom_left.y;
                    k += 2;
                    boxes[k] = d.area.bottom_right.x;
                    boxes[k + 1] = d.area.bottom_right.y;
                    k += 2;
                    boxes[k] = d.area.top_right.x;
                    boxes[k + 1] = d.area.top_right.y;
                    k += 2;
                    boxes[k] = d.area.top_left.x;
                    boxes[k + 1] = d.area.top_left.y;
                    k += 2;
                    rootDistri[d.id] = d;
                    distriLProcTemp = distriLProcTemp->next;
                }


                /* Broadcast the routing table among all processors. */
                MPI_Bcast(&boxes, 11 * psize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


                /* Send its local distribution to every process. */
                distriLProcTemp = distriLProc;
                while (distriLProcTemp != NULL) {
                    distribution d = distriLProcTemp->distri;
                    int receiver = d.id;
                    if (receiver != 0) {
                        std::vector<particule> partsDistri = d.distri;
                        MPI_Send(&partsDistri[0], (partsDistri.size()), mpiParticule, receiver, 0, MPI_COMM_WORLD);
                    } else {
                        // Directly record the local distribution of the root
                        localDistri = d.distri;
                    }
                    distriLProcTemp = distriLProcTemp->next;
                }

                /* Process data that need to be/ */
                delete distriLProc;
                delete distri;
                dAllParticules.clear();
                /**********************************************************************************************************/
            } else {

                /****************** Reception of the routing table and of the local distribution **************************/
                /* Get the rooting table coming from the root */
                MPI_Bcast(&boxes, 11 * psize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                /* Build it. */
                int k(0);
                std::vector<particule> v0 = {particule()};
                while (k < 11 * psize) {
                    int id = boxes[k];
                    k += 1;
                    point cm(boxes[k], boxes[k + 1]);
                    k += 2;
                    point bl(boxes[k], boxes[k + 1]);
                    k += 2;
                    point br(boxes[k], boxes[k + 1]);
                    k += 2;
                    point tr(boxes[k], boxes[k + 1]);
                    k += 2;
                    point tl(boxes[k], boxes[k + 1]);
                    k += 2;
                    distribution *d = new distribution(id, v0, box(bl, br, tr, tl), cm);
                    rootDistri[id] = *d;
                }
                /* Receive local distribution. */
                int nbCoordPartsReceived(receiveDistriComm(0, mpiParticule, MPI_COMM_WORLD, localDistri));

                /******************************************************************************************************/
            }
        } else { /*** In this case we just receive local distributions of the end of the last loop ***/
            for (int p(0) ; p < psize ; p++) {
                if (p != prank) {
                    int nb(receiveDistriComm(p, mpiParticule, MPI_COMM_WORLD, localDistri));
                }
            }
        }

        /* Write the local distribution */
        //writeDistriMPI(fileName,localDistri,prank,int(t/dt),numbPart,psize);

        /************************** Build the local tree which will become the LET ************************************/
        /* Note that here, every local tre has the same top box : it permits to insert particules in it
         * during the building of the LET */
        Tree localTree(b0, particule());
        Tree::buildTree(localTree.getNode(), localDistri);

        /**************************************************************************************************************/



        /**************** Build Local Essential Tree (LET) which will be usefull to compute the forces ****************/

        /* Send to every processor the particules in the localTree that it needs to build its LET */

        std::vector<particule> parts2Send;
        int myproc(prank);

        for (int k(0); k < psize; k++) {
            int otherproc(k);
            if (myproc != otherproc) {
                /* Find the particules in the localTree that belongs to the LET of otherproc. */
                findlocalLET(*localTree.getNode(), rootDistri[otherproc], parts2Send, threshold);

                /* Send these particules to processor otherproc. */
                MPI_Isend(&parts2Send[0], parts2Send.size(), mpiParticule, otherproc, 0, MPI_COMM_WORLD, &sendLET);
            }
            parts2Send.clear();
        }

        /* Receive the particules sent by all the other processors to build LET. */

        int nbPartReceived(0);
        std::vector<particule> distriLET;

        for (int k(0); k < psize; k++) {
            int otherproc(k);
            if (myproc != otherproc) {
                int partsForLET(receiveDistriComm(otherproc,mpiParticule,MPI_COMM_WORLD,distriLET));
                nbPartReceived += partsForLET;
            }
        }

        /* Insert these particules into the localTree to built the LET. */
        for (int d(0); d < nbPartReceived; d++) {
            Tree::insertNode(localTree.getNode(), distriLET[d]);
        }

        dAllParticules.clear();

        /**************************************************************************************************************/



        /******************** Computation of the forces on every particule of localDistri *****************************/

        // Vector which will contain particules that are supposed to go into an other proc
        std::vector<particule> partsOut;
        // Local box in which lie particules in local distri
        box localBox = rootDistri[prank].area;

        /* Compute the new positions after time dt. */
        int nbStayed(moveParticules(localDistri,threshold,localTree,dt,limSquare,partsOut,localBox));

        buildAll = ((int((t+dt)/dt)%load) == 0);

        /* If it is not the last iteration, perform the following communication step */
        if (t + dt <= nbIt * dt) {
            if (buildAll) {
                /*** We'll rebuild all the quadtree during next iteration, so we send all particules to the root */
                // This could be done inside previous function
                for (int i(0); i < partsOut.size(); i++) {
                    localDistri.push_back(partsOut[i]);
                }
                if (prank != 0) {
                    /* Send moved particules to processors 0 to do the loop again. */
                    MPI_Send(&localDistri[0], localDistri.size(), mpiParticule, 0, 0, MPI_COMM_WORLD);
                    localDistri.clear();
                }
            } else {
                /*** We just deal with particules which left their node and went to an other one */
                std::vector<std::vector<particule>> table2Send(psize);
                /* Send every particule in partsOut to the corresponding processor */
                for (int k(0); k < partsOut.size(); k++) {
                    int idReceiver(findProcLie(rootDistri, partsOut[k]));
                    if (idReceiver != -1) {
                        table2Send[idReceiver].push_back(partsOut[k]);
                    }
                }

                for (int k(0); k < psize; k++) {
                    if (k != prank) {
                        MPI_Send(&table2Send[k][0], table2Send[k].size(), mpiParticule, k, 0, MPI_COMM_WORLD);
                    }
                }
            }
        }

        partsOut.clear();

        /**************************************************************************************************************/

        /*** Stop scaling chrono ***/
        if (prank ==0) {
            second elapsedStep = clk::now() - tPerStep;
            //std::printf("step clock time (chrono)        = %.4gs\n", elapsedStep.count());
            averageTimeStep = (elapsedStep.count() + averageTimeStep*(int(t/dt)-1))/(int(t/dt));
        }

        t += dt;

        //MPI_Barrier(MPI_COMM_WORLD);

        //if (prank == 0) {
            //std::cout << "***********************************************" << std::endl;
            //std::cout << " " << std::endl;
        //}

    }


    if (prank ==0) {
        second elapsed = clk::now() - t1;
        std::printf("average step time (chrono)      = %f\n", averageTimeStep);
        std::printf("wall clock time (chrono)        = %.4gs\n", elapsed.count());
    }

    MPI_Type_free(&mpiParticule);
    MPI_Finalize();

}

