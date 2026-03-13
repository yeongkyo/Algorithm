#include <iostream>
#include <list>
#include <string>

using namespace std;

struct Block {
    int id;
    long long len;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    while (cin >> N && N != 0) {
        list<Block> disk;
        disk.push_back({-1, 2000000000000000LL});

        for (int i = 0; i < N; ++i) {
            char cmd;
            cin >> cmd;

            if (cmd == 'W') {
                int id;
                long long s;
                cin >> id >> s;
                
                auto it = disk.begin();
                while (s > 0 && it != disk.end()) {
                    if (it->id == -1) {
                        if (it->len > s) {
                            disk.insert(it, {id, s});
                            it->len -= s;
                            s = 0;
                        } else {
                            it->id = id;
                            s -= it->len;
                        }
                    }
                    ++it;
                }
            } else if (cmd == 'D') {
                int id;
                cin >> id;
                
                for (auto& b : disk) {
                    if (b.id == id) {
                        b.id = -1;
                    }
                }
                
                auto it = disk.begin();
                while (it != disk.end()) {
                    auto next_it = next(it);
                    if (next_it != disk.end() && it->id == -1 && next_it->id == -1) {
                        it->len += next_it->len;
                        disk.erase(next_it);
                    } else {
                        ++it;
                    }
                }
            } else if (cmd == 'R') {
                long long p;
                cin >> p;
                
                long long current_pos = 0;
                int ans = -1;
                
                for (const auto& b : disk) {
                    if (current_pos + b.len > p) {
                        ans = b.id;
                        break;
                    }
                    current_pos += b.len;
                }
                cout << ans << "\n";
            }
        }
        cout << "\n";
    }

    return 0;
}