#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <chrono>
#include <random>

using namespace std;

struct Point {
    long long x, y, z;
    bool operator==(const Point& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

const int M = 524288;
const int MASK = M - 1;

int head[M];
struct Node {
    int idx;
    int next;
} nodes[150005];
int node_cnt = 0;

void init_hash() {
    memset(head, -1, sizeof(head));
    node_cnt = 0;
}

inline unsigned int hash_func(long long cx, long long cy, long long cz) {
    unsigned int h = (unsigned int)(cx * 73856093) ^ (unsigned int)(cy * 19349663) ^ (unsigned int)(cz * 83492791);
    return h & MASK;
}

void insert_hash(long long cx, long long cy, long long cz, int idx) {
    unsigned int h = hash_func(cx, cy, cz);
    nodes[node_cnt] = {idx, head[h]};
    head[h] = node_cnt++;
}

inline long long get_cell(long long x, long long S) {
    if (x >= 0) return x / S;
    return (x - S + 1) / S;
}

inline long long dist(const Point& a, const Point& b) {
    long long dx = a.x - b.x;
    long long dy = a.y - b.y;
    long long dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<Point> pts(N);
    for (int i = 0; i < N; ++i) {
        cin >> pts[i].x >> pts[i].y >> pts[i].z;
    }

    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    });
    pts.erase(unique(pts.begin(), pts.end()), pts.end());
    N = pts.size();

    if (N < 2) {
        return 0;
    }

    mt19937 rng(1337);
    shuffle(pts.begin(), pts.end(), rng);

    long long min_d = dist(pts[0], pts[1]);
    long long S = max(1LL, (long long)sqrt(min_d));

    init_hash();
    insert_hash(get_cell(pts[0].x, S), get_cell(pts[0].y, S), get_cell(pts[0].z, S), 0);
    insert_hash(get_cell(pts[1].x, S), get_cell(pts[1].y, S), get_cell(pts[1].z, S), 1);

    for (int i = 2; i < N; ++i) {
        long long cx = get_cell(pts[i].x, S);
        long long cy = get_cell(pts[i].y, S);
        long long cz = get_cell(pts[i].z, S);

        bool updated = false;
        for (long long dx = -1; dx <= 1; ++dx) {
            for (long long dy = -1; dy <= 1; ++dy) {
                for (long long dz = -1; dz <= 1; ++dz) {
                    unsigned int h = hash_func(cx + dx, cy + dy, cz + dz);
                    for (int p = head[h]; p != -1; p = nodes[p].next) {
                        int j = nodes[p].idx;
                        long long jcx = get_cell(pts[j].x, S);
                        long long jcy = get_cell(pts[j].y, S);
                        long long jcz = get_cell(pts[j].z, S);
                        if (jcx == cx + dx && jcy == cy + dy && jcz == cz + dz) {
                            long long d = dist(pts[i], pts[j]);
                            if (d < min_d) {
                                min_d = d;
                                updated = true;
                            }
                        }
                    }
                }
            }
        }

        if (updated) {
            S = max(1LL, (long long)sqrt(min_d));
            init_hash();
            for (int j = 0; j <= i; ++j) {
                insert_hash(get_cell(pts[j].x, S), get_cell(pts[j].y, S), get_cell(pts[j].z, S), j);
            }
        } else {
            insert_hash(cx, cy, cz, i);
        }
    }

    int min_cnt = 0;
    init_hash();
    S = max(1LL, (long long)sqrt(min_d));

    for (int i = 0; i < N; ++i) {
        long long cx = get_cell(pts[i].x, S);
        long long cy = get_cell(pts[i].y, S);
        long long cz = get_cell(pts[i].z, S);

        for (long long dx = -1; dx <= 1; ++dx) {
            for (long long dy = -1; dy <= 1; ++dy) {
                for (long long dz = -1; dz <= 1; ++dz) {
                    unsigned int h = hash_func(cx + dx, cy + dy, cz + dz);
                    for (int p = head[h]; p != -1; p = nodes[p].next) {
                        int j = nodes[p].idx;
                        long long jcx = get_cell(pts[j].x, S);
                        long long jcy = get_cell(pts[j].y, S);
                        long long jcz = get_cell(pts[j].z, S);
                        if (jcx == cx + dx && jcy == cy + dy && jcz == cz + dz) {
                            long long d = dist(pts[i], pts[j]);
                            if (d == min_d) {
                                min_cnt++;
                            }
                        }
                    }
                }
            }
        }
        insert_hash(cx, cy, cz, i);
    }

    cout << min_d << "\n" << min_cnt << "\n";
    return 0;
}