#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Player {
    int level;
    string name;
};

struct Room {
    int baseLevel;
    vector<Player> players;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int p, m;
    cin >> p >> m;

    vector<Room> rooms;

    for (int i = 0; i < p; i++) {
        int l;
        string n;
        cin >> l >> n;

        bool entered = false;
        for (auto& room : rooms) {
            if (room.players.size() < m && l >= room.baseLevel - 10 && l <= room.baseLevel + 10) {
                room.players.push_back({l, n});
                entered = true;
                break;
            }
        }

        if (!entered) {
            Room new_room;
            new_room.baseLevel = l;
            new_room.players.push_back({l, n});
            rooms.push_back(new_room);
        }
    }

    for (auto& room : rooms) {
        if (room.players.size() == m) {
            cout << "Started!\n";
        } else {
            cout << "Waiting!\n";
        }

        sort(room.players.begin(), room.players.end(), [](const Player& a, const Player& b) {
            return a.name < b.name;
        });

        for (const auto& player : room.players) {
            cout << player.level << " " << player.name << "\n";
        }
    }

    return 0;
}