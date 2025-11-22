#include <vector>
#include "NewDiagram.h"
#include <cassert>

class UF {
    struct UFNode {
        int parent;
        Voronoi::NewDiagram::FacePtr region;
        UFNode* next;
    };

    int cnt;
    std::vector<int> sz;
    std::vector<UFNode> nodes;
    std::vector<UFNode*> head;
    std::vector<UFNode*> tail;

public:
    UF(int N)
        : cnt(N),
          sz(N),
          nodes(N),
          head(N),
          tail(N) {

        for (int i = 0; i < N; i++) {
            nodes[i].parent = i;
            nodes[i].next = nullptr;
            sz[i] = 1;
            head[i] = &nodes[i];
            tail[i] = &nodes[i];
        }
    }

    UFNode& element(int i) {
        return nodes[i];
    }

    void add_UF(int i, Voronoi::NewDiagram::FacePtr reg) {
        // Optional safety check (disable in release)
        assert(i >= 0 && i < (int)nodes.size());

        nodes[i].region = reg;
    }

    int find(int p) {
        while (p != nodes[p].parent) {
            p = nodes[p].parent;
        }
        return p;
    }

    void merge(int x, int y) {
        int i = find(x);
        int j = find(y);
        if (i == j) return;

        // union by size
        if (sz[i] < sz[j]) std::swap(i, j);

        // j joins i
        nodes[j].parent = i;
        sz[i] += sz[j];

        // concatenate lists
        tail[i]->next = head[j];
        tail[i] = tail[j];

        cnt--;
    }

    UFNode* iterate_set(int x) {
        int root = find(x);
        return head[root];  // use next to iterate
    }

    int count() const {
        return cnt;
    }

    int size() const {
        return nodes.size();
    }
};
