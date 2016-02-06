#ifndef AdaptiveDataStructure_h_is_included
#define AdaptiveDataStructure_h_is_included

#include "../DataStructure/DataStructure.h"
#include <set>


template<class Data>
class AdaptiveARRAY : public Array<Data> {
 public:
  AdaptiveARRAY(): Array<Data>() {
  };
  virtual ~AdaptiveARRAY() {Delete();};

  //! Init for BlockAllocator, called in BlockAllocator<AdaptiveARRAY>::NewItem
  virtual void Init(int dim) {this->redim(dim); this->fill(0.0);};
  //! Clear for BlockAllocator, called in the DeleteItem of BlockAllocator<AdaptiveARRAY>
  virtual void Clear() {this->fill(0.0);} ;
  //! Delete for BlockAllocator, called in the destructor BlockAllocator<AdaptiveARRAY>
  virtual void Delete() {this->cleanup();} ;

};

class AdaptiveNodeData {
 public:
  //! default constructor
  AdaptiveNodeData() {index = NULL; surplus = NULL;};
  virtual ~AdaptiveNodeData() {Delete();}

  AdaptiveARRAY<int> *index;
  double *surplus ;

  //! Init for BlockAllocator, called in BlockAllocator<AdaptiveData>::NewItem
  virtual void Init() {};
  //! Clear for BlockAllocator, called in the DeleteItem of BlockAllocator<AdaptiveData>
  virtual void Clear() {index = NULL; surplus = NULL;};
  //! Delete for BlockAllocator, called in the destructor BlockAllocator<AdaptiveData>
  virtual void Delete() {
    if (surplus) {delete[] surplus; surplus = NULL;}
    index = NULL;
  } ;
};

struct AdaptiveNodeDataCompare { bool operator()( AdaptiveNodeData*, AdaptiveNodeData* )  ; } ;


class AdaptiveData {
 public:
  //! default constructor
  AdaptiveData();
  ~AdaptiveData() {Delete();}

  AdaptiveARRAY<int> *index;
  set< AdaptiveNodeData*, AdaptiveNodeDataCompare> NodeData;

  //! Init for BlockAllocator, called in BlockAllocator<AdaptiveData>::NewItem
  virtual void Init();
  //! Clear for BlockAllocator, called in the DeleteItem of BlockAllocator<AdaptiveData>
  virtual void Clear() ;
  //! Delete for BlockAllocator, called in the destructor BlockAllocator<AdaptiveData>
  virtual void Delete() ;
};

struct AdaptiveDataCompare { bool operator()( AdaptiveData*, AdaptiveData* )  ; } ;
#endif
/*Class:AdaptiveDataStructure

NAME:  AdaptiveDataStructure - This is the class derived from the datastructure class.

DESCRIPTION:

  <b>This class defines the datastructure for the implmentation of BlockAllocator.</b>

  - AdatpvieARRAY defines an D-dimensional array to store the coorinates for each interpolation point.

  - AdaptiveNodeData stores the information associated with each multi-index j.

  - AdaptiveData stores the information associated with each mult-index i.

  - AdaptiveNodeDataCompare  a struct used in the standard template librarty <set> to compare two multi-index j

  - AdaptiveDataCompare      a struct used in the standard template librarty <set> to compare two multi-index i

CONSTRUCTORS AND INITIALIZATION:

    The constructor takes no arguments. No further initialization is
    necessary.


MEMBER FUNCTIONS:


  Most member functions are self-explanatory.



*/
