#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MAX_NODES = 1 << 21;
int cnt[MAX_NODES];
long long len[MAX_NODES];

void update(int node, int L, int R, int qL, int qR, int val, const vector<long long>& V) {
    if (qR < L || qL > R) return;
    if (qL <= L && R <= qR) {
        cnt[node] += val;
    } else {
        int mid = L + (R - L) / 2;
        update(2 * node, L, mid, qL, qR, val, V);
        update(2 * node + 1, mid + 1, R, qL, qR, val, V);
    }
    
    if (cnt[node] > 0) {
        len[node] = V[R + 1] - V[L];
    } else {
        if (L == R) len[node] = 0;
        else len[node] = len[2 * node] + len[2 * node + 1];
    }
}

struct Edge {
    long long x;
    int type;
    long long y1;
    long long y2;
    
    bool operator<(const Edge& other) const {
        if (x != other.x) return x < other.x;
        return type > other.type;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    long long W, H;
    int K;
    if (!(cin >> W >> H >> K)) return 0;
    
    vector<Edge> edges;
    vector<long long> V;
    
    for (int i = 0; i < K; ++i) {
        long long f, c, x1, y1, x2, y2;
        cin >> f >> c >> x1 >> y1 >> x2 >> y2;
        
        long long h = H / (c + 1);
        
        vector<pair<long long, long long>> x_ints;
        long long X1_L = max(0LL, f + x1);
        long long X1_R = min(W, f + x2);
        if (X1_L < X1_R) x_ints.push_back({X1_L, X1_R});
        
        long long X2_L = max(0LL, f - x2);
        long long X2_R = min(W, f - x1);
        if (X2_L < X2_R) x_ints.push_back({X2_L, X2_R});
        
        vector<pair<long long, long long>> y_ints;
        for (long long k = 0; k <= c; ++k) {
            long long Y_L, Y_R;
            if (k % 2 == 0) {
                Y_L = k * h + y1;
                Y_R = k * h + y2;
            } else {
                Y_L = (k + 1) * h - y2;
                Y_R = (k + 1) * h - y1;
            }
            Y_L = max(0LL, Y_L);
            Y_R = min(H, Y_R);
            if (Y_L < Y_R) y_ints.push_back({Y_L, Y_R});
        }
        
        for (auto& xi : x_ints) {
            for (auto& yi : y_ints) {
                edges.push_back({xi.first, 1, yi.first, yi.second});
                edges.push_back({xi.second, -1, yi.first, yi.second});
                V.push_back(yi.first);
                V.push_back(yi.second);
            }
        }
    }
    
    if (edges.empty()) {
        cout << W * H << "\n";
        return 0;
    }
    
    sort(V.begin(), V.end());
    V.erase(unique(V.begin(), V.end()), V.end());
    
    sort(edges.begin(), edges.end());
    
    long long total_area = 0;
    long long last_x = edges[0].x;
    
    for (int i = 0; i < (int)edges.size(); ) {
        int j = i;
        long long curr_x = edges[i].x;
        total_area += (curr_x - last_x) * len[1];
        
        while (j < (int)edges.size() && edges[j].x == curr_x) {
            int i1 = lower_bound(V.begin(), V.end(), edges[j].y1) - V.begin();
            int i2 = lower_bound(V.begin(), V.end(), edges[j].y2) - V.begin();
            if (i1 < i2) {
                update(1, 0, V.size() - 2, i1, i2 - 1, edges[j].type, V);
            }
            j++;
        }
        
        last_x = curr_x;
        i = j;
    }
    
    cout << W * H - total_area << "\n";
    
    return 0;
}