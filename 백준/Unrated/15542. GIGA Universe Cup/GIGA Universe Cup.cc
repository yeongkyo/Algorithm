#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

struct Match {
    int u, v;
    bool played;
    int gu, gv;
};

struct TeamStat {
    int id;
    int pts, gd, gf;
};

bool compare_stats(const TeamStat& a, const TeamStat& b) {
    if (a.pts != b.pts) return a.pts > b.pts;
    if (a.gd != b.gd) return a.gd > b.gd;
    return a.gf > b.gf;
}

bool equal_stats(const TeamStat& a, const TeamStat& b) {
    return a.pts == b.pts && a.gd == b.gd && a.gf == b.gf;
}

struct TeamGroup {
    int t[4];
    int sz;
    TeamGroup() { sz = 0; }
    void add(int id) { t[sz++] = id; }
};

struct GroupList {
    TeamGroup g[4];
    int sz;
    GroupList() { sz = 0; }
    void add(const TeamGroup& tg) { g[sz++] = tg; }
};

GroupList resolve_h2h(const TeamGroup& subset, const vector<Match>& matches) {
    TeamStat stats[4];
    for(int i = 0; i < subset.sz; ++i) {
        stats[i].id = subset.t[i];
        stats[i].pts = stats[i].gd = stats[i].gf = 0;
    }
    
    for(const auto& m : matches) {
        int u_idx = -1, v_idx = -1;
        for(int i = 0; i < subset.sz; ++i) {
            if(subset.t[i] == m.u) u_idx = i;
            if(subset.t[i] == m.v) v_idx = i;
        }
        if(u_idx != -1 && v_idx != -1) {
            if(m.gu > m.gv) stats[u_idx].pts += 3;
            else if(m.gu < m.gv) stats[v_idx].pts += 3;
            else { stats[u_idx].pts += 1; stats[v_idx].pts += 1; }
            stats[u_idx].gd += (m.gu - m.gv);
            stats[v_idx].gd += (m.gv - m.gu);
            stats[u_idx].gf += m.gu;
            stats[v_idx].gf += m.gv;
        }
    }
    sort(stats, stats + subset.sz, compare_stats);

    GroupList new_subsets;
    TeamGroup current_subset;
    current_subset.add(stats[0].id);
    for(int i = 1; i < subset.sz; ++i) {
        if (equal_stats(stats[i], stats[i-1])) {
            current_subset.add(stats[i].id);
        } else {
            new_subsets.add(current_subset);
            current_subset.sz = 0;
            current_subset.add(stats[i].id);
        }
    }
    new_subsets.add(current_subset);

    if (new_subsets.sz == 1) {
        GroupList res;
        res.add(subset);
        return res;
    }

    GroupList result;
    for(int i = 0; i < new_subsets.sz; ++i) {
        if(new_subsets.g[i].sz == 1) {
            result.add(new_subsets.g[i]);
        } else {
            GroupList res = resolve_h2h(new_subsets.g[i], matches);
            for(int j = 0; j < res.sz; ++j) {
                result.add(res.g[j]);
            }
        }
    }
    return result;
}

double evaluate(int target, const vector<Match>& matches) {
    TeamStat global_stats[4];
    for(int i = 0; i < 4; ++i) {
        global_stats[i].id = i;
        global_stats[i].pts = global_stats[i].gd = global_stats[i].gf = 0;
    }
    
    for(const auto& m : matches) {
        int u = m.u, v = m.v;
        if (m.gu > m.gv) global_stats[u].pts += 3;
        else if (m.gu < m.gv) global_stats[v].pts += 3;
        else { global_stats[u].pts += 1; global_stats[v].pts += 1; }
        global_stats[u].gd += (m.gu - m.gv);
        global_stats[v].gd += (m.gv - m.gu);
        global_stats[u].gf += m.gu;
        global_stats[v].gf += m.gv;
    }
    sort(global_stats, global_stats + 4, compare_stats);

    GroupList subsets;
    TeamGroup current_subset;
    current_subset.add(global_stats[0].id);
    for(int i = 1; i < 4; ++i) {
        if (equal_stats(global_stats[i], global_stats[i-1])) {
            current_subset.add(global_stats[i].id);
        } else {
            subsets.add(current_subset);
            current_subset.sz = 0;
            current_subset.add(global_stats[i].id);
        }
    }
    subsets.add(current_subset);

    GroupList final_ranking;
    for(int i = 0; i < subsets.sz; ++i) {
        if (subsets.g[i].sz == 1) {
            final_ranking.add(subsets.g[i]);
        } else {
            GroupList res = resolve_h2h(subsets.g[i], matches);
            for(int j = 0; j < res.sz; ++j) {
                final_ranking.add(res.g[j]);
            }
        }
    }

    int advanced = 0;
    double prob = 0.0;
    for(int i = 0; i < final_ranking.sz; ++i) {
        if (advanced >= 2) break;
        const TeamGroup& sub = final_ranking.g[i];
        
        if (advanced + sub.sz <= 2) {
            for(int j = 0; j < sub.sz; ++j) {
                if (sub.t[j] == target) prob = 1.0;
            }
            advanced += sub.sz;
        } else {
            for(int j = 0; j < sub.sz; ++j) {
                if (sub.t[j] == target) {
                    prob = (2.0 - advanced) / sub.sz;
                }
            }
            advanced += sub.sz;
        }
    }
    return prob;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    double P[9];
    double fact[9] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320};
    for(int p = 0; p <= 8; ++p) {
        P[p] = (fact[8] / (fact[p] * fact[8-p])) * pow(0.25, p) * pow(0.75, 8-p);
    }

    int n_records;
    if (!(cin >> n_records)) return 0;
    
    while(n_records--) {
        string s0;
        cin >> s0;
        
        int target = -1;
        for(int i = 0; i < 4; ++i) {
            if(s0[5 * (i + 1)] == '*') target = i;
        }

        vector<Match> matches;
        int u1 = -1, u2 = -1;
        for(int i = 0; i < 4; ++i) {
            string s; 
            cin >> s;
            for(int j = i + 1; j < 4; ++j) {
                string m_str = s.substr(5 * (j + 1), 4);
                Match m;
                m.u = i; m.v = j;
                if(m_str == "__-_") {
                    m.played = false;
                    m.gu = 0; m.gv = 0;
                    matches.push_back(m);
                    if(u1 == -1) u1 = matches.size() - 1;
                    else u2 = matches.size() - 1;
                } else {
                    m.played = true;
                    m.gu = m_str[1] - '0';
                    m.gv = m_str[3] - '0';
                    matches.push_back(m);
                }
            }
        }

        double ans = 0.0;
        for(int x1 = 0; x1 <= 8; ++x1) {
            for(int y1 = 0; y1 <= 8; ++y1) {
                double p1 = P[x1] * P[y1];
                if(p1 == 0) continue;
                matches[u1].gu = x1; 
                matches[u1].gv = y1;
                
                for(int x2 = 0; x2 <= 8; ++x2) {
                    for(int y2 = 0; y2 <= 8; ++y2) {
                        double p = p1 * P[x2] * P[y2];
                        if(p == 0) continue;
                        matches[u2].gu = x2; 
                        matches[u2].gv = y2;
                        
                        ans += p * evaluate(target, matches);
                    }
                }
            }
        }
        cout << fixed << setprecision(7) << ans << "\n";
    }
    return 0;
}