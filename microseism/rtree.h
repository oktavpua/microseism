#ifndef RTREE_H
#define RTREE_H

/**
 * \file   rtree.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  R-дерево с целочисленными координатами.
 */

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <limits>
#include <cstddef>
#include "config.h"
#include "bbox.h"

/// R-дерево.
class RTree {
    public:
	/// Максимальное число потомков в узле.
	const static int node_capacity = RTREE_BRANCHING;

	/// Добавление объекта в индекс.
	void insert(const BoundingBox &bb, size_t oid);

	/// Вывод любого объекта, пересекающегося с заданным прямоугольником.
	int pick(const BoundingBox &bb) const;

    private:
	struct Entry;
	struct LeafEntry;
	struct NodeEntry;

	typedef std::unique_ptr<Entry>     EntryPtr;
	typedef std::unique_ptr<LeafEntry> LeafPtr;
	typedef std::unique_ptr<NodeEntry> NodePtr;

	// Узел дерева.
	struct Entry {
	    BoundingBox bbox;

	    Entry() {};
	    Entry(const BoundingBox &bb) : bbox(bb) {}

	    virtual ~Entry() {};
	};

	// Листовой узел.
	struct LeafEntry : public Entry {
	    LeafEntry(const BoundingBox &bb, size_t oid)
		: Entry(bb), object_id(oid) {}

	    size_t object_id;
	};

	// Внутренний узел.
	struct NodeEntry : public Entry {
	    NodeEntry(int lvl) : level(lvl) {
		child.reserve(node_capacity + 1);
	    }

	    int pick(const BoundingBox &bb) const;

	    NodePtr insert(LeafPtr &&leaf);
	    NodePtr split();

	    int level;
	    std::vector<EntryPtr> child;
	};

	// Корень дерева.
	NodePtr root;
};

#endif
