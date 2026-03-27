#include <bits/stdc++.h>
using namespace std;

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    int N;
    cin >> N;
    vector<int> S(N);
    for (auto &x : S)
        cin >> x;

    int freq[10001] = {0};
    int unique_cnt = 0;
    int left = 0;
    int max_len = 0;

    for (int right = 0; right < N; right++)
    {
        if (freq[S[right]] == 0)
        {
            unique_cnt++;
        }
        freq[S[right]]++;

        while (unique_cnt > 2)
        {
            freq[S[left]]--;
            if (freq[S[left]] == 0)
            {
                unique_cnt--;
            }
            left++;
        }
        max_len = max(max_len, right - left + 1);
    }
    cout << max_len << '\n';
}