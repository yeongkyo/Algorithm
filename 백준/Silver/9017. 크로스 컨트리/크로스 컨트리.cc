#include <iostream>
#include <vector>

using namespace std;

void solve() {
    int n;
    cin >> n;
    
    vector<int> runners(n);
    vector<int> counts(201, 0);
    
    for (int i = 0; i < n; ++i) {
        cin >> runners[i];
        counts[runners[i]]++;
    }

    vector<bool> isValid(201, false);
    for (int i = 1; i <= 200; ++i) {
        if (counts[i] == 6) {
            isValid[i] = true;
        }
    }

    vector<int> scoreSum(201, 0);
    vector<int> validCount(201, 0);
    vector<int> fifthScore(201, 0);

    int currentScore = 1;
    for (int i = 0; i < n; ++i) {
        int team = runners[i];
        if (isValid[team]) {
            validCount[team]++;
            if (validCount[team] <= 4) {
                scoreSum[team] += currentScore;
            } else if (validCount[team] == 5) {
                fifthScore[team] = currentScore;
            }
            currentScore++;
        }
    }

    int minScore = 99999999;
    int winTeam = -1;
    
    for (int i = 1; i <= 200; ++i) {
        if (isValid[i]) {
            if (scoreSum[i] < minScore) {
                minScore = scoreSum[i];
                winTeam = i;
            } else if (scoreSum[i] == minScore) {
                if (fifthScore[i] < fifthScore[winTeam]) {
                    winTeam = i;
                }
            }
        }
    }
    
    cout << winTeam << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int t;
    cin >> t;
    while (t--) {
        solve();
    }
    
    return 0;
}