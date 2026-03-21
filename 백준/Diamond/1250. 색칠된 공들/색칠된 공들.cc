#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <cstdio>

using namespace std;

struct Int24 {
    unsigned char data[3];
    inline int get() const {
        return data[0] | (data[1] << 8) | (data[2] << 16);
    }
    inline void set(int v) {
        data[0] = v & 0xFF;
        data[1] = (v >> 8) & 0xFF;
        data[2] = (v >> 16) & 0xFF;
    }
};

const int MAXM = 10000005;
Int24 tree_nodes[MAXM];
Int24 L_arr[MAXM];
Int24 R_arr[MAXM];
Int24 len_arr[MAXM];
char color[MAXM];

int M = 0;

inline int get_val(int p) {
    return (p >= M) ? (p - M) : tree_nodes[p].get();
}

inline void update_node(int p) {
    int left_id = get_val(2 * p);
    int right_id = get_val(2 * p + 1);
    int l_len = len_arr[left_id].get();
    int r_len = len_arr[right_id].get();
    if (l_len > r_len) {
        tree_nodes[p].set(left_id);
    } else if (l_len < r_len) {
        tree_nodes[p].set(right_id);
    } else {
        tree_nodes[p].set(left_id < right_id ? left_id : right_id);
    }
}

inline void modify(int idx) {
    int p = idx + M;
    p >>= 1;
    for (; p > 0; p >>= 1) {
        update_node(p);
    }
}

char buf[1 << 16];
int buf_idx = 0, buf_sz = 0;
inline char next_char() {
    if (buf_idx == buf_sz) {
        buf_sz = fread(buf, 1, sizeof(buf), stdin);
        buf_idx = 0;
        if (buf_sz == 0) return 0;
    }
    return buf[buf_idx++];
}

inline int next_int() {
    int res = 0;
    char c = next_char();
    while (c <= ' ') {
        if (c == 0) return 0;
        c = next_char();
    }
    while (c > ' ') {
        res = res * 10 + (c - '0');
        c = next_char();
    }
    return res;
}

int main() {
    int N = next_int();
    int k = next_int();
    
    char prev_c = 0;
    int curr_len = 0;
    int current_ball = 1;
    int target_run_id = -1;
    
    char c = next_char();
    while (c <= ' ') {
        if (c == 0) break;
        c = next_char();
    }
    
    while (c > ' ') {
        if (c == prev_c) {
            curr_len++;
        } else {
            if (M > 0) {
                len_arr[M-1].set(curr_len);
            }
            color[M] = c;
            prev_c = c;
            curr_len = 1;
            M++;
        }
        if (current_ball == k) {
            target_run_id = M - 1;
        }
        current_ball++;
        c = next_char();
    }
    if (M > 0) {
        len_arr[M-1].set(curr_len);
    }
    
    for (int i = 0; i < M; ++i) {
        L_arr[i].set(i == 0 ? 0xFFFFFF : i - 1);
        R_arr[i].set(i == M - 1 ? 0xFFFFFF : i + 1);
    }
    
    for (int i = M - 1; i > 0; --i) {
        update_node(i);
    }
    
    int step = 1;
    while (true) {
        int max_id = get_val(1);
        if (len_arr[max_id].get() == 0) break; 
        
        len_arr[max_id].set(0);
        modify(max_id);
        
        if (max_id == target_run_id) {
            printf("%d\n", step);
            return 0;
        }
        
        int l_id = L_arr[max_id].get();
        int r_id = R_arr[max_id].get();
        
        if (l_id != 0xFFFFFF) R_arr[l_id].set(r_id);
        if (r_id != 0xFFFFFF) L_arr[r_id].set(l_id);
        
        if (l_id != 0xFFFFFF && r_id != 0xFFFFFF && color[l_id] == color[r_id]) {
            int sum_len = len_arr[l_id].get() + len_arr[r_id].get();
            len_arr[l_id].set(sum_len);
            
            len_arr[r_id].set(0);
            modify(r_id);
            
            if (target_run_id == r_id) target_run_id = l_id;
            
            int rr_id = R_arr[r_id].get();
            R_arr[l_id].set(rr_id);
            if (rr_id != 0xFFFFFF) L_arr[rr_id].set(l_id);
            
            modify(l_id);
        }
        
        step++;
    }
    
    return 0;
}