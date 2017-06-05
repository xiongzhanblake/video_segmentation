#ifndef __CLASSUTILITIES__
#define __CLASSUTILITIES__

#include "CUtilities.hpp"
#include <stdio.h>
#include <stdlib.h>


//-----------------------------------------------------------------------------------------------//
//								Basic structures & classes					  				 	 //
//-----------------------------------------------------------------------------------------------//

#ifndef __BASICSTRUCTURES__
#define __BASICSTRUCTURES__

template<typename T>
struct UnaryElement
{
public:
	T d;

	UnaryElement(T d) : d(d) {}

	bool operator==(UnaryElement &A)
	{
		if (this->d == A.d) return true;
		else return false;
	}

	bool operator<=(UnaryElement &A)
	{
		if (this->d <= A.d) return true;
		else return false;
	}

	bool operator<(UnaryElement &A)
	{
		if (this->d < A.d) return true;
		else return false;
	}

	bool operator>=(UnaryElement &A)
	{
		if (this->d >= A.d) return true;
		else return false;
	}

	bool operator>(UnaryElement &A)
	{
		if (this->d > A.d) return true;
		else return false;
	}

	// Memory deep copy;
	void operator=(UnaryElement &A)
	{
		this->d = A.d;
	}

	void operator=(T A)
	{
		this->d = A;
	}

	UnaryElement& operator+(UnaryElement& A)
	{
		return (new UnaryElement(this->d + A.d)); 
	}

	UnaryElement& operator-(UnaryElement& A)
	{
		return (new UnaryElement(this->d - A.d)); 
	}

	UnaryElement& operator*(UnaryElement& A)
	{
		return (new UnaryElement(this->d * A.d)); 
	}

	UnaryElement& operator/(UnaryElement& A)
	{
		if (A.d != 0)
			return (new UnaryElement(this->d / A.d)); 
		else
		{
			printf("Divided by zero, exit!\n");
			exit(-1);
		}
	}
};

template<typename T1,typename T2>
struct PairElement
{
public:
	T1 d1;
	T2 d2;

	PairElement(){};

	PairElement(T1 d1, T2 d2) : d1(d1), d2(d2){};

	PairElement(T1 d) : d1(d), d2((T2)d){};

//	PairElement(T2 d) : d1((T1)d), d2(d){};

	bool operator!=(PairElement &A)
	{
		if (this->d1 != A.d1 && this->d2 != A.d2) return true;
		else return false;
	}

	bool operator==(PairElement &A)
	{
		if (this->d1 == A.d1 && this->d2 == A.d2) return true;
		else return false;
	}

	bool operator<=(PairElement &A)
	{
		if (this->d1 <= A.d1 && this->d2 <= A.d2) return true;
		else return false;
	}

	bool operator<(PairElement &A)
	{
		if (this->d1 < A.d1 && this->d2 < A.d2) return true;
		else return false;
	}

	bool operator>=(PairElement &A)
	{
		if (this->d1 >= A.d1 && this->d2 >= A.d2) return true;
		else return false;
	}

	bool operator>(PairElement &A)
	{
		if (this->d1 > A.d1 && this->d2 > A.d2) return true;
		else return false;
	}

	// Memory deep copy;
	void operator=(PairElement &A)
	{
		this->d1 = A.d1;
		this->d2 = A.d2;
	}

//	void operator=(T1 A)
//	{
//		this->d1 = A;
//		this->d2 = (T2)A;
//	}

	void operator=(T2 A)
	{
		this->d1 = (T1)A;
		this->d2 = A;
	}

	PairElement* operator+(PairElement& A)
	{
		return (new PairElement(this->d1 + A.d1, this->d2 + A.d2)); 
	}

	PairElement* operator-(PairElement& A)
	{
		return (new PairElement(this->d1 - A.d1, this->d2 - A.d2)); 
	}

	PairElement* operator*(PairElement& A)
	{
		return (new PairElement(this->d1 * A.d1, this->d2 * A.d2)); 
	}

	PairElement* operator/(PairElement& A)
	{
		if (A.d1 != 0 && A.d2 != 0)
			return (new PairElement(this->d1 / A.d1, this->d2 / A.d2)); 
		else
		{
			printf("Divided by zero, exit!\n");
			exit(-1);
		}
	}
};

template<typename T1,typename T2>
struct MapPair
{
public:
	T1 key;
	T2 value;

	MapPair(){};

	MapPair(T1 d1, T2 d2) : key(d1), value(d2){};

	MapPair(T1 d) : key(d), value((T2)d){};

//	PairElement(T2 d) : d1((T1)d), d2(d){};

	bool operator!=(MapPair &A)
	{
		if (this->key != A.key) return true;
		else return false;
	}

	bool operator==(MapPair &A)
	{
		if (this->key == A.key) return true;
		else return false;
	}

	bool operator<=(MapPair &A)
	{
		if (this->key <= A.key) return true;
		else return false;
	}

	bool operator<(MapPair &A)
	{
		if (this->key < A.key) return true;
		else return false;
	}

	bool operator>=(MapPair &A)
	{
		if (this->key >= A.key) return true;
		else return false;
	}

	bool operator>(MapPair &A)
	{
		if (this->key > A.key) return true;
		else return false;
	}

	// Memory deep copy;
	void operator=(MapPair &A)
	{
		this->key = A.key;
		this->value = A.value;
	}

//	void operator=(T1 A)
//	{
//		this->d1 = A;
//		this->d2 = (T2)A;
//	}

//	void operator=(T2 A)
//	{
//		this->d1 = (T1)A;
//		this->d2 = A;
//	}

	MapPair* operator+(MapPair& A)
	{
		return (new MapPair(this->key, this->value + A.value));
	}

	MapPair* operator-(MapPair& A)
	{
		return (new MapPair(this->key, this->value - A.value));
	}

	MapPair* operator*(MapPair& A)
	{
		return (new MapPair(this->key, this->value * A.value));
	}

	MapPair* operator/(MapPair& A)
	{
		if (A.value != 0)
			return (new MapPair(this->key, this->value / A.value));
		else
		{
			printf("Divided by zero, exit!\n");
			exit(-1);
		}
	}
};

#endif //__BASICSTRUCTURES__

//-----------------------------------------------------------------------------------------------//
//									Class: Memory Manager					  				 	 //
//-----------------------------------------------------------------------------------------------//
#ifndef __MEMORYMANAGER__
#define __MEMORYMANAGER__
/*
 * Memory management for linked block list.
 * Implement adding and deleting elements of the same type in blocks.
 */
template<typename ElementType>
class MemoryManager
{
private:

	typedef struct MemoryBlockItem
	{
		ElementType item;
		MemoryBlockItem *next_free;
	}BlockItem;

	typedef struct MemoryBlock
	{
		MemoryBlockItem *blockitem;
		MemoryBlock *next;
	}Block;

protected:
	Block *headBlock;
	BlockItem *idleList;

	uint64_t	BLOCKSIZE;
public:
	MemoryManager(uint64_t size) {BLOCKSIZE = size; headBlock = NULL; idleList = NULL;};
	~MemoryManager(){while (headBlock){	MemoryBlock *next = headBlock -> next; delete headBlock->blockitem;	delete headBlock; headBlock = next;	}}
	ElementType *New()
	{
		BlockItem *item;
		if (!idleList)
		{
			Block *next = headBlock;
			headBlock = new Block();
			if (headBlock == NULL) {printf("Memory allocationi failed!\n");exit(-1);}
			headBlock->next = next;
			headBlock->blockitem = new BlockItem[BLOCKSIZE];
			if (headBlock->blockitem == NULL) {printf("Memory allocationi failed!\n");exit(-1);}
			for (uint64_t i = 0; i < BLOCKSIZE - 1; i++ )
				headBlock->blockitem[i].next_free = &headBlock->blockitem[i + 1];
			headBlock->blockitem[BLOCKSIZE - 1].next_free = idleList;
			idleList = headBlock->blockitem;
		}
		item = idleList;
		idleList = item->next_free;
		return (ElementType *)item;
	}
	void Delete(ElementType *item)
	{
		((BlockItem *) item) -> next_free = idleList;
		idleList = ((BlockItem *) item);
	}
	void Clear()
	{
		MemoryBlock *traversal = headBlock;
		while (traversal != NULL)
		{
			for (uint64_t i = 0; i < BLOCKSIZE - 1; i++ )
				traversal->blockitem[i].next_free = traversal->blockitem[i + 1];
			if (traversal->next != NULL) traversal->blockitem[BLOCKSIZE - 1].next_free = traversal->next->blockitem;
			traversal = traversal -> next;
		}
	}
};
#endif//__MEMORYMANAGER__

//-----------------------------------------------------------------------------------------------//
//									Class: linked block list				  				 	 //
//-----------------------------------------------------------------------------------------------//
#ifndef __LINKEDBLOCKLIST__
#define __LINKEDBLOCKLIST__

template<typename T>
class LinkedBlockList
{

//The type of data stored in the linked list
private:
	typedef uint64_t SPINT;
	#define BLOCKSIZE 512
	typedef struct LLBlockNodeStruct
	{
		T *item;
		LLBlockNodeStruct *next;
		LLBlockNodeStruct *prev;
	}LLBlockNode;
//typedef MemoryManager<LLBlockNode> MemoryManagerE;
//MemoryManagerE *MemSource;
int conflict_flag;
void (LinkedBlockList::*conflictPtr)(LLBlockNode *conflict, T item);
// Logical structures of Linked node list;
// For block traversal, pointers to logical
// first, last and current node in the list.
	// LLBlockNode *middle;
public:
	LLBlockNode *cursor;
	LLBlockNode *begin, *end;
	SPINT length;
	// SPINT Indmiddle;

//-----------------------------------------------------------------------//
//	       			  Implementation of member functions 			 	 //
//-----------------------------------------------------------------------//
protected:
	// inline void FindMiddle()
	// {
	// 	SPINT tmp = (SPINT)ceil(length/2);
	// 	if ((tmp - Indmiddle) == 1) { middle = middle->next; Indmiddle ++; }
	// 	else if ((Indmiddle - tmp) == 1) { middle = middle->prev; Indmiddle --; }
	// }
	inline LLBlockNode *New(T item)
	{
//		LLBlockNode *node = MemSource->New();
		LLBlockNode *node = new LLBlockNode();
		node->item = new T();
		*node->item = item;
		return node;
	}

	inline void conflict_remove(LLBlockNode *conflict, T item)
	{
		conflict -> prev -> next = conflict -> next;
		conflict -> next -> prev = conflict -> prev;
//		MemSource -> Delete(conflict);
		delete conflict;
	}

	inline void conflict_default(LLBlockNode *conflict, T item)
	{ 
		printf("The item has existed, do nothing and return!\n"); 
	}

	inline void conflict_addition(LLBlockNode *conflict, T item)
	{ 
		T *tmp = *conflict->item + item;
		delete conflict->item;
		conflict->item = tmp;
	}

	inline void conflict_substract(LLBlockNode *conflict, T item)
	{ 
		T *tmp = *conflict->item - item;
		T *zero = new T(0);
		if (*tmp != *zero)
		{
			delete conflict->item;
			conflict->item = tmp;
		}
		else
			conflict_remove(conflict, item);
		delete zero;
	}

	inline void conflict_multiply(LLBlockNode *conflict, T item)
	{
		T *tmp = *conflict->item * item;
		delete conflict->item;
		conflict->item = tmp;
	}

	inline void conflict_divide(LLBlockNode *conflict, T item)
	{
		T *tmp = *conflict->item / item;
		delete conflict -> item;
		conflict -> item = tmp;
	}

	LLBlockNode *popFront()
		{
			if (!isEmpty())
			{
				setCursorFront();
				LLBlockNode *pos = begin;
				if (hasNext()) begin->next->prev = NULL;
				begin = begin->next;
				length--;
				pos->next = pos->prev = NULL;
				// if (middle == pos) { middle = pos -> prev; Indmiddle --;}
				// FindMiddle();
				return pos;
			}
			else
				return NULL;
		}

	LLBlockNode *popEnd()
	{
		if (!isEmpty())
		{
			setCursorEnd();
			LLBlockNode *pos = end;
			if (hasPrev()) end->prev->next = NULL;
			end = end->prev;
			length--;
			pos->next = pos->prev = NULL;
			// if (middle == pos) { middle = pos -> prev; Indmiddle --;}
			// FindMiddle();
			return pos;
		}
		else return NULL;
	}

private:
	
	void addFront(LLBlockNode* item)
	{
		if (!isEmpty())
		{
			item->prev = NULL;
			item->next = begin;
			begin->prev = item;
			begin = item;
		}
		else
		{
			end = begin = item;
			item->prev = item->next = NULL;
		}
		length++;
		// FindMiddle();
	}

	void addEnd(LLBlockNode* item)
	{
		if (!isEmpty())
		{
			item->prev = end;
			item->next = NULL;
			end->next = item;
			end = item;
		}
		else
		{
			begin = end = item;
			item->prev = item->next = NULL;
		}
		length++;
		// FindMiddle();
	}

	void insert(LLBlockNode *pos, LLBlockNode *hit)
	{
		hit->next = pos->next;
		hit->prev = pos;
		pos->next->prev = hit;
		pos->next = hit;
		length++;
		// FindMiddle();
	}

	// When current link is empty, add all elements of L at end;
	// Otherwise, find position and insert;
	// When position is hit, conflict operation will be invoked;
	void copyLinkdedBlockList(LinkedBlockList* L)
	{
		if (this->isEmpty())
		{
			for (LLBlockNode *ptr = L->begin; ptr != NULL; ptr = ptr -> next)
				this->addEnd(*ptr->item);
		}
		else
		{
			for (LLBlockNode *ptr = L->begin; ptr != NULL; ptr = ptr -> next)
				this->insert(*ptr->item);
		}
	}

public:

	T* MoveTo(SPINT step)
	{
		SPINT dist = 0;
		while (cursor != NULL && dist < step)
		{
			cursor = cursor -> next;
			dist++;
		}
		if (dist == step && cursor != NULL)
			return cursor->item;
		else
		{
			return NULL;
		}
	}

	T* MoveBack(SPINT step)
	{
		SPINT dist = 0;
		while (cursor != NULL && dist < step)
		{
			cursor = cursor -> prev;
			dist++;
		}
		if (dist == step && cursor != NULL)
			return cursor->item;
		else
		{
			return NULL;
		}
	}

	// addFront and addEnd will not cheak repeat elements in the link;
	void addFront(T data)
	{
		LLBlockNode *item = this->New(data);
		addFront(item);
	}

	void addEnd(T data)
	{
		LLBlockNode *item = this->New(data);
		addEnd(item);
	}
	// find repeat elements meanwhile insert new elements
	void insert(T data)
	{
		if (!isEmpty())
		{
			LLBlockNode *pos = NULL;
			SPINT step = 0;
			if(!search(data, &pos, &step))
			{
				LLBlockNode *hit = this->New(data);
				if (pos == NULL && step == 0) addFront(hit);
				else if (pos == NULL && step == this->length)	addEnd(hit);
				else insert(pos,hit);
			}
			else (this->*conflictPtr)(pos, data);
		}
		else
			this->addEnd(data);
	}

	LinkedBlockList *operator+(LinkedBlockList *L)
	{
		LinkedBlockList *tmp = new LinkedBlockList();
		tmp->setConfilicFlag(1);
		tmp->copyLinkdedBlockList(this);
		tmp->copyLinkdedBlockList(L);
		return tmp;
	}

	LinkedBlockList *operator-(LinkedBlockList *L)
	{
		LinkedBlockList *tmp = new LinkedBlockList();
		tmp->setConfilicFlag(2);
		tmp->copyLinkdedBlockList(this);
		tmp->copyLinkdedBlockList(L);
		return tmp;
	}

	LinkedBlockList *operator*(LinkedBlockList *L)
	{
		LinkedBlockList *tmp = new LinkedBlockList();
		tmp->setConfilicFlag(3);
		tmp->copyLinkdedBlockList(this);
		tmp->copyLinkdedBlockList(L);
		return tmp;
	}

	LinkedBlockList *operator/(LinkedBlockList *L)
	{
		LinkedBlockList *tmp = new LinkedBlockList();
		tmp->setConfilicFlag(4);
		tmp->copyLinkdedBlockList(this);
		tmp->copyLinkdedBlockList(L);
		return tmp;
	}

	void operator=(LinkedBlockList *L)
	{
		if (!this->isEmpty()) this->clear();
		this->copyLinkdedBlockList(L);
	}

	T *popHead()
	{
		LLBlockNode *ret = popFront();
		if (ret == NULL) return NULL;
		else return ret->item;
	}

	T *popRear()
	{
		LLBlockNode *ret = popEnd();
		if (ret == NULL) return NULL;
		else return ret->item;
	}

	bool search(T item, LLBlockNode **result, SPINT *step)
	{
		if (item < *begin->item) { *result = NULL; *step = 0; return false;}
		if (item > *end->item) { *result = NULL; *step = this->length; return false;}
		for(cursor = begin; cursor != NULL; cursor = cursor -> next, step++)
		{
				if (*cursor->item <  item) continue;
				if (*cursor->item == item) {*result = cursor; return true;}
				if (*cursor->item >  item) break;
		}
		*result = cursor -> prev;
		return false;
	}

	inline bool isEmpty()
	{
		if (begin == NULL)
		{
			cursor = end = begin;
			return true;
		}
		else
			return false;
	};

	inline LinkedBlockList(int conflict_flag = 0)
	{
		cursor = begin = end = NULL;
		length = 0;
		this->conflict_flag = conflict_flag;
		conflictPtr = NULL;
		setConfilicFlag(conflict_flag);
	};

	inline void setConfilicFlag(int conflict_flag)
	{
		this->conflict_flag = conflict_flag;
		switch (conflict_flag)
		{
			case 0: conflictPtr = &LinkedBlockList<T>::conflict_default; break;
			case 1: conflictPtr = &LinkedBlockList<T>::conflict_addition; break;
			case 2: conflictPtr = &LinkedBlockList<T>::conflict_substract; break;
			case 3: conflictPtr = &LinkedBlockList<T>::conflict_multiply; break;
			case 4: conflictPtr = &LinkedBlockList<T>::conflict_divide; break;
			case 5: conflictPtr = &LinkedBlockList<T>::conflict_remove; break;
			default: conflictPtr = &LinkedBlockList<T>::conflict_default; break;
		}
	}

	inline void setCursorFront()
	{
		cursor = begin;
	};

	inline void setCursorEnd()
	{
		cursor = end;
	};

	// inline void setCursorMiddle()
	// {
	// 	cursor = middle;
	// };

	~LinkedBlockList()
	{
		this->clear();
		cursor = begin = end = NULL;
		length = 0;
	};

	inline LLBlockNode* next()
	{
		return(cursor = cursor -> next);
	}

	inline LLBlockNode* prev()
	{
		return(cursor = cursor -> prev);
	}

	inline bool hasNext()
	{
		if ( cursor->next != NULL ) return true;
		else return false;
	}

	inline bool hasPrev()
	{
		if ( cursor->prev != NULL ) return true;
		else return false;
	}

	inline void clear()
	{
		LLBlockNode* tmp = NULL;
		while((tmp = popFront()))
		{
			delete tmp->item;
			tmp->item = NULL;
			tmp->prev = tmp->next = NULL;
			delete tmp;
		}
	}
};

#endif //__LINKEDBLOCKLIST__

//-----------------------------------------------------------------------------------------------//
//										Class: sparse matrix 	 				 	 			 //
//-----------------------------------------------------------------------------------------------//

#ifndef __SPARSEMATRIX__
#define __SPARSEMATRIX__

template <typename T, typename P>
class SparseMatrix
{
private:
#define BLOCKSIZE 512
#define SafeDelete(ptr) if(ptr != NULL){delete ptr; ptr = NULL;}

typedef uint64_t SPINT;
typedef MapPair<T,P> PairT;
typedef LinkedBlockList<PairT> LBLP;
typedef struct SparseMatrixElement
{
	P item;
	LBLP *link;
	SparseMatrixElement(P item, LBLP *link = NULL) : item(item), link(link){};
	SparseMatrixElement() : link(NULL){};
	void setElement(P item, LBLP *link) {this->item = item; this->link = link;}
	SparseMatrixElement* operator+(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item + A.item, *this->link + A.link));}

	SparseMatrixElement* operator-(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item - A.item, *this->link - A.link));}

	SparseMatrixElement* operator*(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item * A.item, *this->link * A.link));}

	SparseMatrixElement* operator/(SparseMatrixElement& A)
	{
		if (A.item != 0)
			return (new SparseMatrixElement(this->item / A.item, *this->link / A.link));
		else
		{
			printf("Divided by zero, exit!\n");
			exit(-1);
		}
	}
	void operator=(SparseMatrixElement& A)
	{
		this->item = A.item;
		this->link = A.link;
	}
	bool operator!=(SparseMatrixElement& A)
	{
		if (this->item != A.item) return true;
		else return false;
	}

}SPMElement;

//MemoryManager<LBLP> *LLIcPv;
//MemoryManager<SPMElement> *LLIr;
typedef LinkedBlockList<SPMElement> LBLE;

public:
	SPINT junk;
	SPINT *counter;
	SPINT total;
	SPINT rRows;
	SPINT rows;
	SPINT cols;
	SPMElement *matrix;
// For data interface to other functions;
	P *Ir;
	P *Ic;
	T *Pv;
//-----------------------------------------------------------------------//
//	         		Implementation of member functions 			 		 //
//-----------------------------------------------------------------------//

protected:

void findValidRows(P *Ir)
{
	SPINT *tmp = CMemoryInit<SPINT,SPINT>(0, this->total);
	SPINT index = 0;
	P record = Ir[0];
	for (SPINT i = 0; i < this->total; i++)
	{
		if (record == Ir[i]) tmp[index]++;
		else
		{
			record = Ir[i];
			index++;
			tmp[index]++;
		}
	}
	this->rRows = index+1;
	this->counter = new SPINT[this->rRows];
	for (SPINT i = 0; i < this->rRows; i++)
		this->counter[i] = tmp[i];
	CMemoryRelease<SPINT>(&tmp);
}

void CreateSparseMatrix(P *Ir, P *Ic, T *Pv)
{
	SPINT indexGlobe = 0;
	this->matrix = (SPMElement *)malloc(sizeof(SPMElement)*this->rRows);
	if (this->matrix == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (SPINT i = 0; i < this->rRows; i++)
	{
		this->matrix[i].item = Ir[indexGlobe];
		this->matrix[i].link = new LBLP();
		for (SPINT j = 0; j < this->counter[i]; j++, indexGlobe++)
		{
			PairT *tmp = new PairT(Ic[indexGlobe],Pv[indexGlobe]);
			this->matrix[i].link->insert(*tmp);
			delete tmp;
		}
	}
}

void CopySparseMatrix(SPMElement *matrix)
{
	LBLP *tmp = NULL;
	if (this->matrix != NULL) clearMatrix();
	this->matrix = (SPMElement *)malloc(sizeof(SPMElement)*this->rRows);
	if (this->matrix == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (SPINT i = 0; i < this->rRows; i++)
	{
		this->matrix[i].item = matrix[i].item;
		tmp = new LBLP();
		this->matrix[i].link = *tmp + matrix[i].link;
		this->total += this->matrix[i].link->length;
		delete tmp; tmp = NULL;
	}
}

SPMElement *RestoreMatrix(LBLE *TMatrix, SPINT length)
{
	SPMElement *ret = (SPMElement *)malloc(sizeof(SPMElement)*length);
	SPMElement *e = NULL;
	SPINT index = 0;
	while (!TMatrix->isEmpty())
	{
		e = TMatrix->popHead();
		ret[index].item = e->item;
		ret[index].link = e->link;
		e->link = NULL; delete e; e = NULL;
		++index;
	}
	return ret;
}

public:

void readSparseMatrix()
{
	Ir = new P[total];
	Ic = new P[total];
	Pv = new T[total];
	if (Ir == NULL || Ic == NULL || Pv == NULL){printf("Memory allocation failed!\n");exit(-1);}
	SPINT index = 0;
//	LBLP *traversal = NULL;
	PairT *e = NULL;
	for (SPINT i = 0; i < rRows; i++)
	{
		this->matrix[i].link -> setCursorFront();
		while((e = this->matrix[i].link -> MoveTo(1)))
		{
			Ir[index] = this->matrix[i].item;
			Ic[index] = e -> key;
			Pv[index] = e -> value;
			index ++;
		}
	}
}

P* getIr()
{
	if (Ir == NULL) readSparseMatrix();
	return Ir;
}

P* getIc()
{
	if (Ic == NULL) readSparseMatrix();
	return Ic;
}


T* getPv()
{
	if (Pv == NULL) readSparseMatrix();
	return Pv;
}

SparseMatrix(SPMElement *matrix, P length, P rows, P cols)
{
	this->rows = (SPINT)rows;
	this->cols = (SPINT)cols;
	this->rRows = (SPINT)length;
	this->total = 0;
	this->junk = 0;
//	LLIcPv = new MemoryManager<LBLP>(BLOCKSIZE);
//	LLIr = new MemoryManager<SPMElement>(BLOCKSIZE);
	counter = NULL;
	CopySparseMatrix(matrix);
	Ic = NULL;
	Ir = NULL;
	Pv = NULL;
}


SparseMatrix(P *sp_i, P *sp_j, T *sp_v, P length, P rows, P cols)
{
	Ir = NULL;
	Ic = NULL;
	Pv = NULL;
	this->rows = (SPINT)rows;
	this->cols = (SPINT)cols;
	this->rRows = 0;
	this->total = (SPINT)length;
	this->junk = 0;
	matrix = NULL;
	counter = NULL;
	findValidRows(sp_i);
	CreateSparseMatrix(sp_i, sp_j, sp_v);
}

~SparseMatrix()
{
	clearMatrix();
	matrix = NULL;
//	delete LLIr;
//	delete LLIcPv;
}

void clearMatrix()
{
	for (SPINT i = 0; i < this->rRows; i++)
		this->matrix->link->clear();
	free(this->matrix);
	this->matrix = NULL;
	SafeDelete(Ir);
	SafeDelete(Ic);
	SafeDelete(Pv);
}

void operator=(SparseMatrix &A)
{
	if (A.rows != this->rows || A.cols != this->cols) {printf("Dimensions are mismatched!\n");exit(-1);}
	this->rRows = A.rRows;
	this->clearMatrix();
	CopySparseMatrix(A.matrix,A.rRows,A.rows,A.cols);
}

SparseMatrix* operator+(const SparseMatrix& A)
{
	if (A.rows != this->rows || A.cols != this->cols) {printf("Dimensions are mismatched!\n");exit(-1);}

	LBLE *TMatrix = new LBLE();
	SPMElement *e = NULL;
	LBLP *link = NULL;
	LBLP *trash = NULL;
	SPINT ia = 0, ib = 0;
	SPINT TotalRows = 0;
	while (ia < this->rRows && ib < A.rRows)
	{
		if (this->matrix[ia].item < A.matrix[ib].item)
		{
			trash = new LBLP();
			link = *trash + this->matrix[ia].link;
			delete trash; trash = NULL;
			e = new SPMElement();
			e ->setElement(this->matrix[ia].item, link);
			link = NULL;
			TMatrix->addEnd(*e);
			e->link = NULL; delete e; e = NULL;
			TotalRows ++;
			ia ++;
		}
		else if (this->matrix[ia].item == A.matrix[ib].item)
		{
			link = *this->matrix[ia].link + A.matrix[ib].link;
			e = new SPMElement();
			e->setElement(this->matrix[ia].item, link);
			link = NULL;
			TMatrix->addEnd(*e);
			e->link = NULL; delete e; e = NULL;
			TotalRows ++;
			ia ++; ib ++;
		}
		else if (this->matrix[ia].item > A.matrix[ib].item)
		{
			trash = new LBLP();
			link = *trash + A.matrix[ib].link;
			delete trash; trash = NULL;
			e = new SPMElement();
			e->setElement(A.matrix[ib].item, link);
			link = NULL;
			TMatrix->addEnd(*e);
			e->link = NULL; delete e; e = NULL;
			TotalRows ++;
			ib ++;
		}
	}

	while (ia < this->rRows)
	{
		trash = new LBLP();
		link = *trash + this->matrix[ia].link;
		delete trash; trash = NULL;
		e = new SPMElement();
		e->setElement(this->matrix[ia].item, link);
		link = NULL;
		TMatrix->addEnd(*e);
		e->link = NULL; delete e; e = NULL;
		TotalRows ++;
		ia ++;
	}

	while (ib < A.rRows)
	{
		trash = new LBLP();
		link = *trash + A.matrix[ib].link;
		delete trash; trash = NULL;
		e = new SPMElement();
		e->setElement(A.matrix[ib].item, link);
		link = NULL;
		TMatrix->addEnd(*e);
		e->link = NULL; delete e; e = NULL;
		TotalRows ++;
		ib ++;
	}
	SPMElement *tmp = RestoreMatrix(TMatrix, TotalRows);
	SparseMatrix *ret = new SparseMatrix(tmp, (P)TotalRows, (P)this->rows, (P)this->cols);
	delete TMatrix; TMatrix = NULL;
	delete tmp; tmp = NULL;
	return ret;
}

};
#endif //__SPARSEMATRIX__

#endif //__CLASSUTILITIES__

