#include "AdaptiveDataStructure.h"

// Implementations
//
// functions for BlockAllocator

AdaptiveData::AdaptiveData(){
	
  index = NULL;
}

void AdaptiveData::Init(){
	
}

void AdaptiveData::Clear(){
     index = NULL;
	 NodeData.clear();
}

void AdaptiveData::Delete()
{  
	 index = NULL;
	 NodeData.clear();
	 
}

bool AdaptiveNodeDataCompare::operator()( AdaptiveNodeData* a,  AdaptiveNodeData* b)
{
	
 for ( int i = 1; i<= a->index->n ; i++){
	 if ( (*(a->index))(i) < (*(b->index))(i) ) return true;
	 if ( (*(a->index))(i) > (*(b->index))(i) ) return false;
 }
 return false;
}

bool AdaptiveDataCompare::operator()( AdaptiveData* a,  AdaptiveData* b)
{
	
 for ( int i = 1; i<= a->index->n ; i++){
	 if ( (*(a->index))(i) < (*(b->index))(i) ) return true;
	 if ( (*(a->index))(i) > (*(b->index))(i) ) return false;
 }
 return false;
}