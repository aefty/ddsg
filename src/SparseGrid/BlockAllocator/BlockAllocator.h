#ifndef BlockAllocator_h_is_included
#define BlockAllocator_h_is_included

#include <list>
using namespace std;

template<class I>
class BlockAllocator {
public:
    //! default constructor
    BlockAllocator() {;}
    //! default destructor, deletes the pending, previously allocated items
    ~BlockAllocator();
    //! allocate an item, I::Init() has been called at least once in the lifetime of the object, but not necessarily in this
    //! call of NewItem
    I*    NewItem(int dim) ;
    I*    NewItem();
    //! release item and call I::Clear()
    void  DeleteItem(I* p) ;
    void  ClearAll();

    list<I*> List   ;
    list<I*> Blocks ;

    //! define the blocksize of the object. you can modify it according to your size of memory.
    enum { BLOCKSIZE=10000 } ;
};

template<class I>
BlockAllocator<I>::~BlockAllocator() {

    ClearAll();
}

template<class I>
void BlockAllocator<I>::ClearAll() {
    I *p ;
    while (!Blocks.empty()) {
        p=Blocks.front() ;
        for (int i=0; i<BLOCKSIZE; i++) p[i].Delete() ;
        delete[] p ;
        Blocks.pop_front() ;
    }
    List.clear();
}


template<class I>
I* BlockAllocator<I>::NewItem(int dim) {
    I *p ;

    if (List.empty()) {
        p=new I[BLOCKSIZE];
        Blocks.push_front(p) ;
        for (int i=0; i<BLOCKSIZE; i++) {
            p[i].Init( dim) ;
            List.push_front(&p[i]) ;
        }
    }

    // return an init. instance
    p=List.front() ;
    List.pop_front() ;
    return p ;
}

template<class I>
I* BlockAllocator<I>::NewItem() {
    I *p ;

    if (List.empty()) {
        p=new I[BLOCKSIZE];
        Blocks.push_front(p) ;
        for (int i=0; i<BLOCKSIZE; i++) {
            p[i].Init() ;
            List.push_front(&p[i]) ;
        }
    }

    // return an init. instance
    p=List.front() ;
    List.pop_front() ;
    return p ;
}

template<class I>
void BlockAllocator<I>::DeleteItem(I *p) {
    p->Clear() ;
    List.push_front(p) ;
}
#endif
/*Class:BlockAllocator

NAME:  BlockAllocator

DESCRIPTION:

  //! BlockAllocator provides a basic mechanism for the heap friendly dynamic allocation and de-allocation
  //! of many instances of I. Requirements on I: <br>
  //! * must have a member function void Init () which mimics a default constructor <br>
  //! * must have a member function void Clear() which clears proprietary members of I<br>
  //! * must have a member function void Delete() which mimics a destructor <br>

CONSTRUCTORS AND INITIALIZATION:

    The constructor takes no arguments. No further initialization is
    necessary.


MEMBER FUNCTIONS:


  Most member functions are self-explanatory.


End:
*/ 