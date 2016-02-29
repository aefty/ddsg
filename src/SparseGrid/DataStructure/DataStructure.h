#ifndef DataStructure_h_IS_INCLUDED
#define DataStructure_h_IS_INCLUDED

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
using std::cout;
using std::endl;
using std::string;

#ifdef _MSC_VER // if using Microsoft C++ compiler
    #pragma warning( disable : 4786 ) // disable warning C4786
#endif


//! Class of array data structure
template<class Data>
class Array {
  public:
    Array() {
        n = 0;
        pData = NULL;
    };
    virtual ~Array() {
        cleanup();
    };
    virtual void cleanup() {
        if (pData) {
            delete[] pData;

            pData = NULL;
            n = 0;
        }
    }
    Array& operator=(const Array& A) {
        redim(A.n );
        int i;
        for (i = 0; i < n; i++) {
            pData[i] = A.pData[i];
        }
        return *this;
    }

    virtual void redim(int n_) {
        cleanup();
        this->n = n_;
        pData = new Data[n];
    }
    Data& operator () (int i) {
        return pData[i - 1];
    }
    virtual void fill(Data data) {
        int i;
        for (i = 0; i < n; i++) {
            pData[i] = data;
        }
    }
    void add(Array<Data>& deta, Data ratio_) {
        int i;
        for (i = 0; i < n; i++) {
            pData[i] += deta.pData[i] * ratio_;
        }
    }
    Data norm() {
        Data sum = 0;
        int i;
        for (i = 0; i < n; i++) {
            sum += pData[i] * pData[i];
        }
        return sqrt(sum);
    }
    void print(int nrow = 2) {
        int i;
        for (i = 0; i < n; i++) {
            //printf("%e ", pData[i]);
            cout << pData[i] << " ";
            if ( i % nrow == (nrow - 1) ) {
                printf("\n");
            }
        }
    }
    Data* pData;
    int n;
};



//! Class of matrix data structure
template<class T>
class MATRIX1 {
  public:
    MATRIX1() {
        nx = 0; ny = 0;
        pData = NULL;
    };
    MATRIX1(int nx_, int ny_) {
        nx = 0; ny = 0;
        pData = NULL;
        redim(nx_, ny_);
    };
    virtual ~MATRIX1() {
        cleanup();
    };
    void cleanup() {
        if (pData) {
            delete[] pData;
            pData = NULL;
        }
        pData = NULL;
    }
    virtual void redim(int nx_, int ny_) {
        cleanup();
        this->nx = nx_;
        this->ny = ny_;
        pData = new T[nx * ny];

    }
    virtual T& operator () (int i, int j) {
        return pData[ (i - 1) * ny + j - 1 ];
    }

    T* pData;
    int nx;
    int ny;
};



//! Class of matrix data structure (including operations)
template <class T>
class Matrix1: public MATRIX1<T> {

  public:
    using MATRIX1<T>::nx;
    using MATRIX1<T>::ny;
    using MATRIX1<T>::pData;
    using MATRIX1<T>::cleanup;
    using MATRIX1<T>::redim;
    using MATRIX1<T>::operator();
    Matrix1(): MATRIX1<T>() {
    };

    Matrix1(int nx, int ny): MATRIX1<T>(nx, ny) {
    };
    ~Matrix1() {
    };

    Matrix1& operator=(Matrix1& A) {
        redim(A.nx, A.ny );
        for (int i = 0; i < nx * ny; i++) {
            pData[i] = A.pData[i];
        }
        return *this;
    }

    virtual void print() {
        int i, j;
        for (i = 1; i <= nx; i++) {
            for (j = 1; j <= ny; j++) {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }

    virtual void print(char* filename) {
        FILE* fp = fopen( filename, "a");
        int i, j;
        for (i = 1; i <= nx; i++) {
            for (j = 1; j <= ny; j++) {
                double val = (*this)(i, j);
                fprintf(fp, "%e ", val );
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
        fclose(fp);
    }

    virtual void fill(T data) {
        int i;
        for (i = 0; i < nx * ny; i++) {
            pData[i] = data;
        }
    }
};
#endif
/*Class:DataStructure

NAME:  DataStructure - Basic data structure used in this tool box

DESCRIPTION:

  This class defines two simple data structures : Matrix and Array.


CONSTRUCTORS AND INITIALIZATION:

  Before using the elements of the Marix and Array, you must allocate
  memory by calling
        Array.redim(dim);
        Matrix.redim(dimx,dimy);


MEMBER FUNCTIONS:

   Most member functions are self-explanatory.


End:
*/
