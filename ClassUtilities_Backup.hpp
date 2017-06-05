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
		T item;
		LLBlockNodeStruct *next;
		LLBlockNodeStruct *prev;
	}LLBlockNode;
typedef MemoryManager<LLBlockNode> MemoryManagerE;
MemoryManagerE *MemSource;
int conflict_flag;
void (LinkedBlockList::*conflictPtr)(LLBlockNode *conflict, T item);
// Logical structures of Linked node list;
// For block traversal, pointers to logical
// first, last and current node in the list.
	// LLBlockNode *middle;
public:
	LLBlockNode *cursor;
	LLBlockNode begin, end;
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
		LLBlockNode *node = MemSource->New();
		node->item = item;
		return node;
	}

	inline void conflict_remove(LLBlockNode *conflict, T item)
	{
		conflict -> prev -> next = conflict -> next;
		conflict -> next -> prev = conflict -> prev;
		MemSource -> Delete(conflict);
	}

	inline void conflict_default(LLBlockNode *conflict, T item)
	{ 
		printf("The item has existed, do nothing and return!\n"); 
	}

	inline void conflict_addition(LLBlockNode *conflict, T item)
	{ 
		T *tmp = conflict -> item + item;
		conflict -> item = *tmp;
		delete tmp;
	}

	inline void conflict_substract(LLBlockNode *conflict, T item)
	{ 
		T *tmp = conflict -> item - item;
		T *zero = new T(0);
		if (*tmp != *zero)
			conflict -> item = *tmp;
		else
			conflict_remove(conflict, item);
		delete tmp;
		delete zero;
	}

	inline void conflict_multiply(LLBlockNode *conflict, T item)
	{
		T *tmp = conflict -> item * item;
		conflict -> item = *tmp;
		delete tmp;
	}

	inline void conflict_divide(LLBlockNode *conflict, T item)
	{
		T *tmp = conflict -> item / item;
		conflict -> item = *tmp;
		delete tmp;
	}

private:
	
	void addFront(LLBlockNode* item)
	{
		item->next = begin.next;
		item->prev = &begin;
		begin.next->prev = item;
		begin.next = item;
		length++;
		// FindMiddle();
	}

	void addEnd(LLBlockNode* item)
	{
		item->prev = end.prev;
		item->next = &end;
		end.prev->next = item;
		end.prev = item;
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
			for (LLBlockNode *ptr = L->begin.next; ptr != &(L->end); ptr = ptr -> next)
				this->addEnd(ptr->item);
		}
		else
		{
			for (LLBlockNode *ptr = L->begin.next; ptr != &(L->end); ptr = ptr -> next)
				this->insert(ptr->item);
		}
	}

public:

	T* MoveTo(SPINT step)
	{
		SPINT dist = 0;
		while (cursor != & end && dist < step)
		{
			cursor = cursor -> next;
			dist++;
		}
		if (dist == step && cursor != & end)
			return &cursor->item;
		else
		{
//			printf("Access beyond the range of this link, return NULL!\n");
			return NULL;
		}
	}

	T* MoveBack(SPINT step)
	{
		SPINT dist = 0;
		while (cursor != & begin && dist < step)
		{
			cursor = cursor -> prev;
			dist++;
		}
		if (dist == step)
			return &cursor->item;
		else
		{
			printf("Access beyond the range of this link, return NULL!\n");
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
			if(!search(data, &pos))
			{
				LLBlockNode *hit = this->New(data);
				if (pos == &begin) addFront(hit);
				else if (pos == &end)	addEnd(hit);
				else insert(pos,hit);
			}
			else (this->*conflictPtr)(pos, data);
		}
		else
			this->addFront(data);
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

	LLBlockNode *popFront()
	{
		LLBlockNode *pos = begin.next;
		begin.next = pos -> next;
		pos -> next -> prev = &begin;
		length--;
		// if (middle == pos) { middle = pos -> prev; Indmiddle --;}
		// FindMiddle();
		return pos;
	}

	LLBlockNode *popEnd()
	{
		LLBlockNode *pos = end.prev;
		end.prev = pos -> prev;
		pos -> prev -> next = &end;
		length--;
		// if (middle == pos) { middle = pos -> prev; Indmiddle --;}
		// FindMiddle();
		return next;
	}

	bool search(T item, LLBlockNode **result)
	{
		if (item < begin.next->item) { *result = &begin; return false;}
		if (item > end.prev->item) { *result = &end; return false;}
		for(cursor = begin.next; cursor != &end; cursor = cursor -> next)
		{
				if (cursor->item <  item) continue;
				if (cursor->item == item) {*result = cursor; return true;}
				if (cursor->item >  item) break;
		}
		*result = cursor -> prev;
		return false;
	}

	inline bool isEmpty()
	{
		if (begin.next == &end)
			return(true);
		else
			return(false);
	};

	inline LinkedBlockList(int conflict_flag = 0)
	{
		cursor = &begin;
		begin.next = &end;
		end.prev = &begin;
		length = 0;
		MemSource = new MemoryManagerE(BLOCKSIZE);
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
		cursor = &begin;
	};

	inline void setCursorEnd()
	{
		cursor = &end;
	};

	// inline void setCursorMiddle()
	// {
	// 	cursor = middle;
	// };

	~LinkedBlockList()
	{
		cursor = &begin;
		begin.next = &end;
		end.prev = &begin;
		delete MemSource;
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
		if ( cursor->next != 0 ) return (true);
		else return(false);
	}

	inline bool hasPrev()
	{
		if ( cursor->prev != 0 ) return (true);
		else return(false);
	}

	inline void clear(){ MemSource -> clear(); begin.next = end.prev = NULL;}
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
//typedef struct SparseMatrixElement
//{
//	P item;
//	LBLP *link;
//	SparseMatrixElement(P item, LBLP *link = NULL) : item(item), link(link){};
//	SparseMatrixElement() : link(NULL){};
//	void setElement(P item, LBLP *link) {this->item = item; this->link = link;}
//	SparseMatrixElement* operator+(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item + A.item, *this->link + A.link));}
//
//	SparseMatrixElement* operator-(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item - A.item, *this->link - A.link));}
//
//	SparseMatrixElement* operator*(SparseMatrixElement& A) { return (new SparseMatrixElement(this->item * A.item, *this->link * A.link));}
//
//	SparseMatrixElement* operator/(SparseMatrixElement& A)
//	{
//		if (A.item != 0)
//			return (new SparseMatrixElement(this->item / A.item, *this->link / A.link));
//		else
//		{
//			printf("Divided by zero, exit!\n");
//			exit(-1);
//		}
//	}
//	void operator=(SparseMatrixElement& A)
//	{
//		this->item = A.item;
//		this->link = A.link;
//	}
//	bool operator!=(SparseMatrixElement& A)
//	{
//		if (this->item != A.item) return true;
//		else return false;
//	}
//
//}SPMElement;

MemoryManager<LBLP> *LLIcPv;
//MemoryManager<SPMElement> *LLIr;
//typedef LinkedBlockList<SPMElement> LBLE;

public:
	SPINT junk;
	SPINT *counter;
	SPINT total;
	SPINT rRows;
	SPINT rows;
	SPINT cols;
	LBLP *matrix;
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
	this->matrix = new SPMElement*[this->rRows];
	if (this->matrix == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (SPINT i = 0; i < this->rRows; i++)
	{
		matrix[i].item = Ir[indexGlobe];
		matrix[i].link = new LBLP();
		for (SPINT j = 0; j < this->counter[i]; j++, indexGlobe++)
			matrix[i].link->insert(*(new PairT(Ic[indexGlobe],Pv[indexGlobe])));
	}
}

void CopySparseMatrix(SPMElement *matrix)
{
	if (this->matrix != NULL) delete this->matrix;
	this->matrix = new SPMElement*[this->rRows];
	if (this->matrix == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (SPINT i = 0; i < this->rRows; i++)
	{
		this->matrix[i].item = matrix[i].item;
		this->matrix[i].link = new LBLP();
		this->matrix[i].link = *this->matrix[i].link + matrix[i].link;
		this->total += this->matrix[i].link->length;
	}
}

SPMElement *RestoreMatrix(LBLE *TMatrix, SPINT length)
{
	SPMElement *ret = new SPMElement[length];
	SPMElement *e = NULL;
	TMatrix -> setCursorFront();
	for (SPINT i = 0; i < length; i++)
	{
		e = TMatrix -> MoveTo(1);
		ret[i].item = e->item;
		ret[i].link = LLIcPv->New();
		ret[i].link = *ret[i].link + e->link;
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
	LLIcPv = new MemoryManager<LBLP>(BLOCKSIZE);
	LLIr = new MemoryManager<SPMElement>(BLOCKSIZE);
	counter = NULL;
	CopySparseMatrix(matrix);
	Ic = NULL;
	Ir = NULL;
	Pv = NULL;
}


SparseMatrix(P *sp_i, P *sp_j, T *sp_v, P length, P rows, P cols)
{
	// Ir = new P[length];
	// Ic = new P[length];
	// Pv = new T[length];
	// if (Ir == NULL || Ic == NULL || Pv == NULL){printf("Memory allocation failed!\n");exit(-1);}
	// CMemoryCopy<P>(Ir,sp_i);
	// CMemoryCopy<P>(Ic,sp_j);
	// CMemoryCopy<P>(Pv,sp_v);
	Ir = NULL;
	Ic = NULL;
	Pv = NULL;
	this->rows = (SPINT)rows;
	this->cols = (SPINT)cols;
	this->rRows = 0;
	this->total = (SPINT)length;
	this->junk = 0;
	LLIcPv = new MemoryManager<LBLP>(BLOCKSIZE);
	LLIr = new MemoryManager<SPMElement>(BLOCKSIZE);
	matrix = NULL;
	counter = NULL;
	findValidRows(sp_i);
	CreateSparseMatrix(sp_i, sp_j, sp_v);
}

~SparseMatrix()
{
	SafeDelete(Ir);
	SafeDelete(Ic);
	SafeDelete(Pv);
	delete []matrix;
//	delete LLIr;
//	delete LLIcPv;
}

void operator=(SparseMatrix &A)
{
	if (A.rows != this->rows || A.cols != this->cols) {printf("Dimensions are mismatched!\n");exit(-1);}
	this->rRows = A.rRows;
	LLIr->clear();
	LLIcPv->clear();
	CopySparseMatrix(A.matrix,A.rRows,A.rows,A.cols);
	SafeDelete(Ir);
	SafeDelete(Ic);
	SafeDelete(Pv);
}

SparseMatrix* operator+(const SparseMatrix& A)
{
	if (A.rows != this->rows || A.cols != this->cols) {printf("Dimensions are mismatched!\n");exit(-1);}
	SPINT minLength = min(this->rRows, A.rRows);
	SPINT maxLength = max(this->rRows, A.rRows);
	
	LBLE *TMatrix = new LBLE();
	SPMElement *e = NULL;
	LBLP *link = NULL;
	SPINT ia = 0, ib = 0;
	SPINT TotalRows = 0;
	for (P i = 0; i < minLength; i++)
	{
		if (this->matrix[ia].item < A.matrix[ib].item)
		{
			link = LLIcPv->New();
			link = *link + this->matrix[ia].link;
			e = LLIr->New();
			e ->setElement(this->matrix[ia].item, link);
			TMatrix->addEnd(*e);
			TotalRows ++;
			ia ++;
		}
		else if (this->matrix[ia].item == A.matrix[ib].item)
		{
			link = LLIcPv->New();
			link = *this->matrix[ia].link + A.matrix[ib].link;
			e = LLIr->New();
			e->setElement(this->matrix[ia].item, link);
			TMatrix->addEnd(*e);
			TotalRows ++;
			ia ++; ib ++;
		}
		else if (this->matrix[ia].item > A.matrix[ib].item)
		{
			link = LLIcPv->New();
			link = *link + A.matrix[ib].link;
			e = LLIr->New();
			e->setElement(A.matrix[ib].item, link);
			TMatrix->addEnd(*e);
			TotalRows ++;
			ib ++;
		}
	}

	while (ia < this->rRows)
	{
		link = LLIcPv->New();
		link = *link + this->matrix[ia].link;
		e = LLIr->New();
		e->setElement(this->matrix[ia].item, link);
		TMatrix->addEnd(*e);
		TotalRows ++;
		ia ++;
	}

	while (ib < A.rRows)
	{
		link = LLIcPv->New();
		link = *link + A.matrix[ib].link;
		e = LLIr->New();
		e->setElement(A.matrix[ib].item, link);
		TMatrix->addEnd(*e);
		TotalRows ++;
		ib ++;
	}
	SPMElement *tmp = RestoreMatrix(TMatrix, TotalRows);
	SparseMatrix *ret = new SparseMatrix(tmp, (P)TotalRows, (P)this->rows, (P)this->cols);
	delete TMatrix;
	delete tmp;
	return ret;
}

};
#endif //__SPARSEMATRIX__

#endif //__CLASSUTILITIES__

