#include <stdio.h>
#include "graph.h"

// #define INFINITE_D 1000000000		/* infinite distance to the terminal */

//#define TERMINAL ( (arc *) 1 )		/* terminal */
//#define ORPHAN   ( (arc *) 2 )		/* orphan */

#define INFINITE_D 1000000000		/* infinite distance to the terminal */

template <typename captype, typename tcaptype, typename flowtype> 
inline void Graph<captype,tcaptype,flowtype>::set_active(node_id i)
{
	if (!i->next)
	{
		/* it's not in the list yet */
		if (queue_last[1]) queue_last[1] -> next = i;
		else               queue_first[1]        = i;
		queue_last[1] = i;
		i -> next = i;
	}
}


template <typename captype, typename tcaptype, typename flowtype> 
void Graph<captype,tcaptype,flowtype>::augment(arc_id middle_arc)
{
	node_id i = NULL;
	arc_id a = NULL;
	tcaptype bottleneck;
	nodeptr *np = NULL;

	/* 1. Finding bottleneck capacity */
	/* 1a - the source tree */
	bottleneck = middle_arc -> r_cap;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > a->sister->r_cap) bottleneck = a -> sister -> r_cap;
	}
	if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
	/* 1b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
	}
	if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


	/* 2. Augmenting */
	/* 2a - the source tree */
	middle_arc -> sister -> r_cap += bottleneck;
	middle_arc -> r_cap -= bottleneck;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		a -> r_cap += bottleneck;
		a -> sister -> r_cap -= bottleneck;
		if (!a->sister->r_cap)
		{
			i -> parent = ORPHAN;
			np = nodeptr_block -> New();
			np -> ptr = i;
			np -> next = orphan_first;
			orphan_first = np;
		}
	}
	i -> tr_cap -= bottleneck;
	if (!i->tr_cap)
	{
		i -> parent = ORPHAN;
		np = nodeptr_block -> New();
		np -> ptr = i;
		np -> next = orphan_first;
		orphan_first = np;
	}
	/* 2b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		a -> sister -> r_cap += bottleneck;
		a -> r_cap -= bottleneck;
		if (!a->r_cap)
		{
			i -> parent = ORPHAN;
			np = nodeptr_block -> New();
			np -> ptr = i;
			np -> next = orphan_first;
			orphan_first = np;
		}
	}
	i -> tr_cap += bottleneck;
	if (!i->tr_cap)
	{
		i -> parent = ORPHAN;
		np = nodeptr_block -> New();
		np -> ptr = i;
		np -> next = orphan_first;
		orphan_first = np;
	}


	flow += bottleneck;
}

template <typename captype, typename tcaptype, typename flowtype> 
void Graph<captype,tcaptype,flowtype>::process_source_orphan(node_id i)
{
	node_id j;
	arc_id a0, a0_min = NULL, a;
	int d, d_min = INFINITE_D;
	nodeptr *np;

	/* trying to find a new parent */
	for (a0 = i->first; a0; a0 = a0->next)
	if (a0->sister->r_cap)
	{
		j = a0 -> head;
		if (j->tree == S && (a=j->parent))
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==ORPHAN) { d = INFINITE_D; break; }
				j = a -> head;
			}
			if (d<INFINITE_D) /* j originates from the source - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=a0->head; j->TS!=TIME; j=j->parent->head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if ((i->parent = a0_min))
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		i -> TS = 0;
		i -> tree = F;
		/* process neighbors */
		for (a0=i->first; a0; a0=a0->next)
		{
			j = a0 -> head;
			if (j->tree == S && (a=j->parent))
			{
				if (a0->sister->r_cap) set_active(j);
				if (a!=TERMINAL && a!=ORPHAN && a->head==i)
				{
					j -> parent = ORPHAN;
					np = nodeptr_block -> New();
					np -> ptr = j;
					if (orphan_last) orphan_last -> next = np;
					else             orphan_first        = np;
					orphan_last = np;
					np -> next = NULL;
				}
			}
		}
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
	void Graph<captype,tcaptype,flowtype>::process_sink_orphan(node_id i)
{
	node_id j;
	arc_id a0, a0_min = NULL, a;
	int d, d_min = INFINITE_D;
	nodeptr *np;

	/* trying to find a new parent */
	for (a0=i->first; a0; a0=a0->next)
	if (a0->r_cap)
	{
		j = a0 -> head;
		if (j->tree == T && (a=j->parent))
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==ORPHAN) { d = INFINITE_D; break; }
				j = a -> head;
			}
			if (d<INFINITE_D) /* j originates from the sink - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=a0->head; j->TS!=TIME; j=j->parent->head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if ((i->parent = a0_min))
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		i -> TS = 0;
		i -> tree = F;
		/* process neighbors */
		for (a0=i->first; a0; a0=a0->next)
		{
			j = a0 -> head;
			if (j->tree == T && (a=j->parent))
			{
				if (a0->r_cap) set_active(j);
				if (a!=TERMINAL && a!=ORPHAN && a->head==i)
				{
					j -> parent = ORPHAN;
					np = nodeptr_block -> New();
					np -> ptr = j;
					if (orphan_last) orphan_last -> next = np;
					else             orphan_first        = np;
					orphan_last = np;
					np -> next = NULL;
				}
			}
		}
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
inline typename Graph<captype,tcaptype,flowtype>::node_id Graph<captype,tcaptype,flowtype>::next_active()
{
	node_id i;

	while ( 1 )
	{
		if (!(i=queue_first[0]))
		{
			queue_first[0] = i = queue_first[1];
			queue_last[0]  = queue_last[1];
			queue_first[1] = NULL;
			queue_last[1]  = NULL;
			if (!i) return NULL;
		}

		/* remove it from the active list */
		if (i->next == i) queue_first[0] = queue_last[0] = NULL;
		else              queue_first[0] = i -> next;
		i -> next = NULL;

		/* a node in the list is active iff it has a parent */
		if (i->parent) return i;
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
void Graph<captype,tcaptype,flowtype>::maxflow_init()
{
	node_id i;
	node_block *nb;

	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = NULL;

	TIME = 0;

	for (nb = node_block_first; nb; nb = nb -> next)
		for (i = &nb->nodes[0]; i < nb -> current; i++)
		{
			i -> next = NULL;
			i -> TS = TIME;
			if (i->tr_cap > 0)
			{
				/* i is connected to the source */
				i -> tree = S;
				i -> parent = TERMINAL;
				set_active(i);
				i -> DIST = 1;
			}
			else if (i->tr_cap < 0)
			{
				/* i is connected to the sink */
				i -> tree = T;
				i -> parent = TERMINAL;
				set_active(i);
				i -> DIST = 1;
			}
			else
			{
				i -> tree = F;
				i -> parent = NULL;
			}
		}
}

template <typename captype, typename tcaptype, typename flowtype> 
flowtype Graph<captype,tcaptype,flowtype>::maxflow()
{
	node_id i, j, current_node = NULL;
	arc_id a;
	nodeptr *np, *np_next;

//	maxflow_init();
	nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);

	while ( 1 )
	{
		if ((i=current_node))
		{
			i -> next = NULL; /* remove active flag */
			if (!i->parent) i = NULL;
		}
		if (!i)
			if (!(i = next_active())) break;
		/* 1. growth stage */
		if (i->tree == S)
		{
		/* 1a. grow source tree */
			for (a = i->first; a; a = a->next)
			{
				if (a -> r_cap)
				{
					j = a -> head;
					if (j->tree == F)
					{
						j -> tree = i -> tree;
						j -> parent = a -> sister;
						j -> TS = i -> TS;
						j -> DIST = i -> DIST + 1;
						set_active(j);
					}
					else if (j->tree == T) break;
					else if (j->TS <= i->TS && j->DIST > i->DIST)
					{
						/* heuristic - trying to make the distance from j to the source shorter */
						j -> parent = a -> sister;
						j -> TS = i -> TS;
						j -> DIST = i -> DIST + 1;
					}
				}
			}
		}
		else if (i->tree == T)
		{
			/* 1b. grow sink tree */
			for (a = i->first; a; a = a->next)
			{
				if ((a -> sister -> r_cap))
				{
					j = a -> head;
					if (j->tree == F)
					{
						j -> tree = i -> tree;
						j -> parent = a -> sister;
						j -> TS = i -> TS;
						j -> DIST = i -> DIST + 1;
						set_active(j);
					}
					else if (j->tree == S) { a = a -> sister; break; }
					else if (j->TS <= i->TS && j->DIST > i->DIST)
					{
						/* heuristic - trying to make the distance from j to the sink shorter */
						j -> parent = a -> sister;
						j -> TS = i -> TS;
						j -> DIST = i -> DIST + 1;
					}
				}
			}
		}

		TIME ++;

		if (a)
		{
			i -> next = i; /* set active flag */
			current_node = i;

			/* augmentation */
			augment(a);
			/* augmentation end */

			/* adoption */
			while (( np = orphan_first ))
			{
				np_next = np -> next;
				np -> next = NULL;

				while ((np = orphan_first))
				{
					orphan_first = np -> next;
					i = np -> ptr;
					nodeptr_block -> Delete(np);
					if (!orphan_first) orphan_last = NULL;
					if (i->tree == T) process_sink_orphan(i);
					else if (i->tree == S) process_source_orphan(i);
				}

				orphan_first = np_next;
			}
			/* adoption end */
		}
		else current_node = NULL;
	}

	delete nodeptr_block; nodeptr_block = NULL;
	maxflow_iteration ++;
	return flow;
}
