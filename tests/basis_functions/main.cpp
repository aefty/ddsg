#include <iostream>
#include <ostream>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip> // setprecision

#include "../../core/HDMR/lib/SGwrite.cpp"
#include "../../core/HDMR/lib/SGread.cpp"

using namespace std;

void funNull(double* x ,  int xDim, double* val ) {
    val[0] = x[0];
}

// Global vars
string outputFolder = "";
int mpirank ;
int mpisize;


int main(int argc, char* argv[]) {

    int SGgridType;
    int SGmaxLevel;

    // General predefined parameters not so important
    int dim = 2;
    int dof = 1;
    double SGcutOff = 0.0;
    double dataPoints = 11.0;

    MPI_Init (&argc, &argv);
    if (argc < 3) {
        cout << "./main <SGgridType><SGmaxLevel><folderName>" << endl;
        MPI_Finalize ();
        return 0;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpisize);

    SGgridType      = atoi(argv[1]);
    SGmaxLevel      = atoi(argv[2]);
    dataPoints      = atoi(argv[3]);
    outputFolder    = argv[4];

    if (mpirank == 0) {
        fstream index_output;
        index_output.open(outputFolder + "/index.txt", ios::out | ios::trunc);

        index_output << "==================================" << endl;
        index_output << "SGmaxLevel: "     << SGmaxLevel     << endl;
        index_output << "SGgridType: "     << SGgridType     << endl;
        index_output << "==================================" << endl;
        cout << "Output index.txt : Successful" << endl;
        index_output.close();


        // Pre processing tests
        SGwrite sgwrite = SGwrite(funNull, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);



        SGread sgread = SGread(0);
        sgread.mpiCOMM = MPI_COMM_WORLD;
        sgread.rank = mpirank;
        sgread.size = mpisize;


        double h = 1.0 / dataPoints;
        int l = SGmaxLevel;
        double temp ;
        vector<double> x(dataPoints + 1);


        fstream xLine_output;
        xLine_output.open(outputFolder + "/X.csv", ios::out | ios::trunc);
        for (int k = 0 ; k <= dataPoints; k++) {
            x[k] = k / dataPoints;
            xLine_output << x[k] << "\t";
        }
        xLine_output.close();


        fstream data_output;
        fstream data_output_post;
        data_output.open(outputFolder + "/basisFunction.csv", ios::out | ios::trunc);
        data_output.open(outputFolder + "/basisFunctionPost.csv", ios::out | ios::trunc);

        if (l == 1) {

            cout <<  "Level :" << l << endl;
            cout << "--------" << endl;
            for (int k = 0 ; k < x.size(); k++) {
                temp = sgwrite.BasisFunction(x[k], l, 0 );
                //cout << "f(" << x[k] << ")="  << temp << endl;
                data_output <<  temp << "\t";
                data_output_post << sgread.BasisFunction(x[k], l, 0 ) << endl;
            }
            data_output << "\n";
        } else {

            cout << "Level :" << l << endl;
            cout << "--------" << endl;
            for (int i = 2; i <= (1 << l); i = i + 2) { // index

                cout <<  endl;
                cout <<  "Index :" << i << endl;
                cout << "--------" << endl;

                for (int k = 0 ; k < x.size(); k++) {
                    temp = sgwrite.BasisFunction(x[k], l, i);
                    cout << "f(" << x[k] << ")="  << temp << endl;
                    data_output <<  temp << "\t";
                }

                data_output << "\n";
            }
        }
        data_output.close();


        MPI_Barrier(MPI_COMM_WORLD);
    }



    MPI_Finalize ();


    return 0;
}

