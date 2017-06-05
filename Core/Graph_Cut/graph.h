// One combination version of GCMEX 2.0 and GCMEX 3.0
// In this header file, use block to restore nodes and arcs

#ifndef __GRAPH_H__
#define __GRAPH_H__

// #include <tmwtypes.h>
#include "block.h"
#include <iostream>
#include <assert.h>
#include "tmwtypes.h"
// NOTE: in UNIX you need to use -DNDEBUG preprocessor option to supress assert's!!!

// captype: type of edge capacities (excluding t-links)
// tcaptype: type of t-links (edges between nodes and terminals)
// flowtype: type of total flow
//
// Current instantiations are in instances.inc
template <typename captype, typename tcaptype, typename flowtype> class Graph
{
protected:

	struct node;
	struct arc;

	typedef struct node* node_id;
	typedef struct arc*  arc_id;
	
	arc_id TERMINAL; /*	special constants for node->parent */
	arc_id ORPHAN;

	typedef enum
	{
		SOURCE	= 0,
		SINK	= 1,
	} termtype; /* terminals */

	typedef enum
	{
		S = 0,
		T = 1,
		F = 2
	} treetype; /* Trees */

	typedef enum arctype
	{
    	FOR = 0,
    	REV = 1,
  	} arctype;

public:

	/* interface functions */

	/* Constructor. Optional argument is the pointer to the
	   function which will be called if an error occurs;
	   an error message is passed to this function. If this
	   argument is omitted, exit(1) will be called. */
	Graph(void (*err_function)(const char *) = NULL);

	/* Destructor */
	~Graph();

	/* Adds a node to the graph */
	node_id add_node();

	/* Adds a bidirectional edge between 'from' and 'to'
	   with the weights 'cap' and 'rev_cap' */
	void add_edge(node_id from, node_id to, captype cap, captype rev_cap);

	// Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights.
	// Can be called multiple times for each node.
	// Weights can be negative.
	// NOTE: the number of such edges is not counted in edge_num_max.
	//       No internal memory is allocated by this call.
	void add_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink);
	void set_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink);
	// Computes the maxflow. Can be called several times.
	// FOR DESCRIPTION OF reuse_trees, SEE mark_node().
	// FOR DESCRIPTION OF changed_list, SEE remove_from_changed_list().
	flowtype maxflow();
	// After the maxflow is computed, this function returns to which
	// segment the node 'i' belongs (Graph<captype,tcaptype,flowtype>::SOURCE or Graph<captype,tcaptype,flowtype>::SINK).
	//
	// Occasionally there may be several minimum cuts. If a node can be assigned
	// to both the source and the sink, then default_segm is returned.
	termtype what_segment(node_id i, termtype default_segm = SOURCE);

	//////////////////////////////////////////////
	//       ADVANCED INTERFACE FUNCTIONS       //
	//      (provide access to the graph)       //
	//////////////////////////////////////////////

public:

	////////////////////////////////////////////////////////////////////////////////
	// 1. Functions for getting pointers to arcs and for reading graph structure. //
	//    NOTE: adding new arcs may invalidate these pointers (if reallocation    //
	//    happens). So it's best not to add arcs while reading graph structure.   //
	////////////////////////////////////////////////////////////////////////////////

	// The following two functions return arcs in the same order that they
	// were added to the graph. NOTE: for each call add_edge(i,j,cap,cap_rev)
	// the first arc returned will be i->j, and the second j->i.
	// If there are no more arcs, then the function can still be called, but
	// the returned arc_id is undetermined.

	// arc_id get_first_arc();
	// arc_id get_next_arc(arc_id a);

	// other functions for reading graph structure
	// int get_node_num() { return node_num; }
	// int get_arc_num() { return arc_num; }
	void get_arc_ends(arc_id a, node_id& i, node_id& j); // returns i,j to that a = i->j

	///////////////////////////////////////////////////
	// 2. Functions for reading residual capacities. //
	///////////////////////////////////////////////////

	// returns residual capacity of SOURCE->i minus residual capacity of i->SINK
	tcaptype get_trcap(node_id i); 
	// returns residual capacity of arc a
	captype get_rcap(arc_id a);

	/////////////////////////////////////////////////////////////////
	// 3. Functions for setting residual capacities.               //
	//    NOTE: If these functions are used, the value of the flow //
	//    returned by maxflow() will not be valid!                 //
	/////////////////////////////////////////////////////////////////

	void set_trcap(node_id i, tcaptype trcap); 
	void set_rcap(arc_id a, captype rcap);

	/////////////////////////////////////////////////////////////////
	// 4. Functions for init node and arc.                         //
	/////////////////////////////////////////////////////////////////
	void readNode(node_id n);
	void setNodeID(node_id n, int id) {n->siteid = id;};
	void node_init(node_id n);
	void arc_init(arc_id a, arc_id b, node_id c, node_id d, captype cap, captype rev_cap);
/////////////////////////////////////////////////////////////////////////
///                  Internal variables and functions                  //
/////////////////////////////////////////////////////////////////////////

private:
//#define int64_T long long
#define NODE_BLOCK_SIZE 512
#define ARC_BLOCK_SIZE 1024
#define NODEPTR_BLOCK_SIZE 128

	typedef struct node
	{
		int			siteid;
		arc			*first;		// first outcoming arc

		arc			*parent;	// node's parent
		node		*next;		// pointer to the next active node
								//   (or to itself if it is the last node in the list)
		int			TS;			// timestamp showing when DIST was computed
		int			DIST;		// distance to the terminal
		treetype	tree;		// flag showing whether the node is in the source or in the sink tree (if parent!=NULL)

		tcaptype	tr_cap;		// if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
								// otherwise         -tr_cap is residual capacity of the arc node->SINK 
	} node;

	typedef struct arc
	{
		node		*head;		// node the arc points to
		arc			*next;		// next arc with the same originating node
		arc			*sister;	// reverse arcs
		arctype     direction;  // arc direction
		captype		r_cap;		// residual capacity
	} arc;

	typedef struct nodeptr_st
	{
		node    	*ptr;
		nodeptr_st	*next;
	} nodeptr;

	typedef struct node_block_st
	{
		node					*current;
		struct node_block_st	*next;
		node					nodes[NODE_BLOCK_SIZE];
	} node_block;

	typedef struct arc_block_st
	{
		arc 					*current;
		struct arc_block_st     *next;
		arc						arcs[ARC_BLOCK_SIZE];
	} arc_block;

	int64_T				node_num;
	int64_T				arc_num;
	int					maxflow_iteration;
	node_block			*node_block_first;
	arc_block			*arc_for_block_first;
	arc_block			*arc_rev_block_first;
	DBlock<nodeptr>		*nodeptr_block;
	void	(*error_function)(const char *);	// this function is called if a error occurs, with a corresponding error message (or exit(1) is called if it's NULL)
	flowtype			flow;					// total flow
	/////////////////////////////////////////////////////////////////////////

	node				*queue_first[2], *queue_last[2];	// list of active nodes
	nodeptr				*orphan_first, *orphan_last;		// list of pointers to orphans
	int64_T				TIME;								// monotonically increasing global counter

	/////////////////////////////////////////////////////////////////////////
protected:
	// functions for processing active list
	void set_active(node_id i);
	node_id next_active();

	// functions for processing orphans list
	void set_orphan_front(node_id i); // add to the beginning of the list
	void set_orphan_rear(node_id i);  // add to the end of the list

	void maxflow_init();             // called if reuse_trees == false
	void augment(arc_id middle_arc);
	void process_source_orphan(node_id i);
	void process_sink_orphan(node_id i);
};


///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::node_init(node_id n)
{
	n -> first = NULL;
	n -> parent = NULL;
	n -> next = NULL;
	n -> TS = 0;
	n -> DIST = 0;
	n -> tr_cap = 0;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::arc_init(arc_id a_for, arc_id a_rev, node_id from, node_id to, captype cap, captype rev_cap)
{
	a_for -> head = to;
	a_rev -> head = from;
	a_for -> r_cap = cap;
	a_rev -> r_cap = rev_cap;
	a_for -> sister = a_rev;
	a_rev -> sister = a_for;
	a_for -> direction = FOR;
	a_rev -> direction = REV;
	a_for -> next = from -> first;
	from -> first = a_for;
	a_rev -> next = to -> first;
	to -> first = a_rev;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline typename Graph<captype,tcaptype,flowtype>::node_id Graph<captype,tcaptype,flowtype>::add_node()
{
	node_id n;

	if (!node_block_first || node_block_first->current+1 > &node_block_first->nodes[NODE_BLOCK_SIZE-1])
	{
		node_block *next = node_block_first;
		node_block_first = new node_block();
		if (!node_block_first) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
		node_block_first -> current = & ( node_block_first -> nodes[0] );
		node_block_first -> next = next;
	}

	n = node_block_first -> current++;
	node_init(n);
	node_num++;
	return n;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::add_edge(node_id from, node_id to, captype cap, captype rev_cap)
{
	arc_id a_for, a_rev;

	if (!arc_for_block_first || arc_for_block_first->current+1 > &arc_for_block_first->arcs[ARC_BLOCK_SIZE])
	{
		arc_block *next = arc_for_block_first;
		arc_for_block_first = new arc_block();
		if (!arc_for_block_first) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
		arc_for_block_first -> current = &( arc_for_block_first -> arcs[0] );
		arc_for_block_first -> next = next;
	}

	if (!arc_rev_block_first || arc_rev_block_first->current+1 > &arc_rev_block_first->arcs[ARC_BLOCK_SIZE])
	{
		arc_block *next = arc_rev_block_first;
		arc_rev_block_first = new arc_block();
		if (!arc_rev_block_first) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
		arc_rev_block_first -> current = &( arc_rev_block_first -> arcs[0] );
		arc_rev_block_first -> next = next;
	}

	a_for = arc_for_block_first -> current++;
	a_rev = arc_rev_block_first -> current++;

	arc_init(a_for, a_rev, from, to, cap, rev_cap);
	arc_num++;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::set_tweights(node_id n, tcaptype cap_source, tcaptype cap_sink)
{
	flow += (cap_source < cap_sink) ? cap_source : cap_sink;
	n -> tr_cap = cap_source - cap_sink;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::add_tweights(node_id n, tcaptype cap_source, tcaptype cap_sink)
{
	tcaptype delta = n -> tr_cap;
	if (delta > 0) cap_source += delta;
	else           cap_sink   -= delta;
	flow += (cap_source < cap_sink) ? cap_source : cap_sink;
	n -> tr_cap = cap_source - cap_sink;	
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::get_arc_ends(arc_id a, node_id& from, node_id& to)
{
	from = a->sister->head;
	to = a->head;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline tcaptype Graph<captype,tcaptype,flowtype>::get_trcap(node_id i)
{
	return i->tr_cap;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline captype Graph<captype,tcaptype,flowtype>::get_rcap(arc_id a)
{
	return a->r_cap;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::set_trcap(node_id i, tcaptype trcap)
{
	i->tr_cap = trcap;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::set_rcap(arc_id a, captype rcap)
{
	a->r_cap = rcap;
}

template <typename captype, typename tcaptype, typename flowtype> 
inline typename Graph<captype,tcaptype,flowtype>::termtype Graph<captype,tcaptype,flowtype>::what_segment(node_id i, termtype default_segm)
{
	if (i->parent)
	{
		if (i->tree == S) return SOURCE;
		else if(i->parent && i->tree == T) return SINK;
	}
	return default_segm;
}

template <typename captype, typename tcaptype, typename flowtype> 
Graph<captype,tcaptype,flowtype>::Graph(void (*err_function)(const char *))
{
	error_function = err_function;
	node_block_first = NULL;
	arc_for_block_first = NULL;
	arc_rev_block_first = NULL;
	flow = 0;
	TIME = 0;
	orphan_first = NULL;
	orphan_last = NULL;
	nodeptr_block = NULL;
	TERMINAL = new arc();
	ORPHAN = new arc();
	arc_num = 0;
	maxflow_iteration = 0;
	node_num = 0;
}

template <typename captype, typename tcaptype, typename flowtype> 
Graph<captype,tcaptype,flowtype>::~Graph()
{
	while (node_block_first)
	{
		node_block *next = node_block_first -> next;
		delete node_block_first;
		node_block_first = next;
	}

	while (arc_for_block_first)
	{
		arc_block *next = arc_for_block_first -> next;
		delete arc_for_block_first;
		arc_for_block_first = next;
	}

	while (arc_rev_block_first)
	{
		arc_block *next = arc_rev_block_first -> next;
		delete arc_rev_block_first;
		arc_rev_block_first = next;
	}

	if (nodeptr_block)
	{
		delete nodeptr_block;
		nodeptr_block = NULL;
	}
	delete TERMINAL;
	delete ORPHAN;
}

template <typename captype, typename tcaptype, typename flowtype>
void Graph<captype,tcaptype,flowtype>::readNode(node_id n)
{
	std::cout<<"-------------------------------------------------------------------------\n";
	std::cout<<"siteID: "<<n->siteid+1<<"\n";
	std::cout<<"tcapacity: "<<n->tr_cap<<"\n";
	std::cout<<"parent -> ";
	if (n -> parent == NULL)
		std::cout<<"NULL\n";
	else if(n -> parent == TERMINAL)
		{
			if(n->tree == S)
				std::cout<<"Source Terminal\n";
			else if (n->tree == T)
				std::cout<<"Sink Terminal\n";
		}
	else if(n -> parent == ORPHAN)
		std::cout<<"Orphan node\n";
	else std::cout<<"siteID: "<<n -> parent -> head -> siteid+1<<"\n";
	std::cout<<"edges(e<i,j>):capacity\n";
	arc_id a = n->first;
	while(a)
	{
		std::cout<<"e("<<n->siteid+1<<", "<< a->head->siteid+1 <<"): "<<a->r_cap<<"\t";
		a = a->next;
	}
	std::cout<<"\n-------------------------------------------------------------------------\n";
}

#endif
