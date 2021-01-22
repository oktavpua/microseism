#include <iomanip>
#include <algorithm>
#include "rtree.h"

using namespace std;

//---------------------------------------------------------------------------
int RTree::NodeEntry::pick(const BoundingBox &bb) const {
    if (level == 0) {
	int max_id = -1;

	for(auto c = child.begin(); c != child.end(); c++) {
	    LeafEntry *leaf = static_cast<LeafEntry*>(c->get());
	    if (leaf->bbox.overlaps(bb)) max_id = std::max<int>(max_id, leaf->object_id);
	}

	return max_id;
    } else {
	int oid = -1;

	for(auto c = child.begin(); c != child.end(); c++) {
	    NodeEntry *node = static_cast<NodeEntry*>(c->get());
	    if (node->bbox.overlaps(bb)) {
		oid = std::max(oid, node->pick(bb));
	    }
	}

	return oid;
    }
}

//---------------------------------------------------------------------------
RTree::NodePtr RTree::NodeEntry::insert(LeafPtr &&leaf) {
    bbox.stretch(leaf->bbox);

    if (level == 0) {
	child.push_back(move(leaf));
    } else {
	struct GrowArea {
	    int area;
	    int grow;

	    GrowArea()
		: area(numeric_limits<int>::max()),
		  grow(numeric_limits<int>::max())
	    {}

	    GrowArea(int a1, int a2)
		: area(a1), grow(a2 - a1) {}

	    bool operator<(GrowArea a) const {
		if (grow < a.grow) return true;
		if (grow > a.grow) return false;
		return area < a.area;
	    }
	};

	NodeEntry *best_node = 0;
	GrowArea best_params;

	for(auto c = child.begin(); c != child.end(); c++) {
	    BoundingBox bb = (*c)->bbox;
	    bb.stretch(leaf->bbox);

	    GrowArea params((*c)->bbox.area(), bb.area());

	    if (params < best_params) {
		best_params = params;
		best_node   = static_cast<NodeEntry*>(c->get());
	    }
	}

	if (auto split_node = best_node->insert(move(leaf)))
	    child.emplace_back(move(split_node));
    }

    return (child.size() > node_capacity) ? split() : NodePtr((NodeEntry*)0);
}

// Overlap, Distribution, Coverage
struct ODC {
    vector<int> list_lo;
    vector<int> list_mi;
    vector<int> list_hi;

    BoundingBox bbox_lo;
    BoundingBox bbox_hi;

    ODC() {
	list_lo.reserve(RTree::node_capacity);
	list_mi.reserve(RTree::node_capacity);
	list_hi.reserve(RTree::node_capacity);
    }

    void reset() {
	list_lo.resize(0);
	list_mi.resize(0);
	list_hi.resize(0);

	bbox_lo = bbox_hi = BoundingBox();
    }

    bool operator<(const ODC &a) const {
	if (overlap < a.overlap) return true;
	if (overlap > a.overlap) return false;

	if (distrib < a.distrib) return true;
	if (distrib > a.distrib) return false;

	return coverage < a.coverage;
    }

    void compute() {
	overlap  = bbox_lo.overlap(bbox_hi);
	distrib  = max(list_lo.size(), list_hi.size());
	coverage = bbox_lo.area() + bbox_hi.area();
    }

    private:
    int overlap;
    int distrib;
    int coverage;
};

//---------------------------------------------------------------------------
RTree::NodePtr RTree::NodeEntry::split() {
    /*
     * Enhanced linear split method from
     * [1] A. Al-Badarneh and M. Tawil, "Linear R-Tree Revisited",
     *     International Journal of Computers and Applications,
     *     Vol. 31, No. 2, 2009
     */

    static array<ODC, BoundingBox::dim> odc;

    // Distribution.
    for(auto d = odc.begin(); d != odc.end(); d++)
	d->reset();

    for(uint i = 0; i < child.size(); i++) {
	for(int d = 0; d < BoundingBox::dim; d++) {
	    int dl = child[i]->bbox.lo[d] - bbox.lo[d];
	    int dr = bbox.hi[d] - child[i]->bbox.hi[d];

	    if (dl < dr) {
		odc[d].list_lo.push_back(i);
		odc[d].bbox_lo.stretch(child[i]->bbox);
	    } else if (dl > dr) {
		odc[d].list_hi.push_back(i);
		odc[d].bbox_hi.stretch(child[i]->bbox);
	    } else {
		odc[d].list_mi.push_back(i);
	    }
	}
    }

    for(auto d = odc.begin(); d != odc.end(); d++) {
	for(auto i = d->list_mi.begin(); i != d->list_mi.end(); i++) {
	    if (d->list_lo.size() < d->list_hi.size()) {
		d->list_lo.push_back(*i);
		d->bbox_lo.stretch(child[*i]->bbox);
	    } else {
		d->list_hi.push_back(*i);
		d->bbox_hi.stretch(child[*i]->bbox);
	    }
	}
	d->compute();
    }

    auto distr = min_element(odc.begin(), odc.end());

    // Splitting.
    NodePtr new_node(new NodeEntry(level));
    new_node->bbox = distr->bbox_hi;

    for(auto i = distr->list_hi.begin(); i != distr->list_hi.end(); i++) {
	new_node->child.push_back(move(child[*i]));
    }

    bbox = distr->bbox_lo;
    sort(distr->list_lo.begin(), distr->list_lo.end());

    std::vector<EntryPtr> new_child;
    new_child.reserve(node_capacity + 1);

    for(auto i = distr->list_lo.begin(); i != distr->list_lo.end(); i++) {
	new_child.push_back(move(child[*i]));
    }

    std::swap(child, new_child);

    return new_node;
}

//---------------------------------------------------------------------------
void RTree::insert(const BoundingBox &bb, size_t oid) {
    LeafPtr new_leaf(new LeafEntry(bb, oid));

    if (root) {
	if (auto split = root->insert(move(new_leaf))) {
	    NodePtr new_root(new NodeEntry(root->level + 1));

	    new_root->bbox = root->bbox;
	    new_root->bbox.stretch(split->bbox);

	    new_root->child.push_back(move(root));
	    new_root->child.push_back(move(split));

	    root = move(new_root);
	}
    } else {
	root.reset(new NodeEntry(0));
	root->child.push_back(move(new_leaf));
	root->bbox = bb;
    }
}

//---------------------------------------------------------------------------
int RTree::pick(const BoundingBox &bb) const {
    if (root->bbox.overlaps(bb))
	return max<int>(0, root->pick(bb));

    return 0;
}

//---------------------------------------------------------------------------
std::ostream& operator<<(std::ostream &os, const BoundingBox &bb) {
    os << "[" << std::setw(3) << bb.lo[0] << ", " << std::setw(3)
       << bb.lo[1] << ", " << std::setw(3) << bb.lo[2] << "] - ["
       << std::setw(3) << bb.hi[0] << ", " << std::setw(3) << bb.hi[1]
       << ", " << std::setw(3) << bb.hi[2] << "]";
    return os;
}

