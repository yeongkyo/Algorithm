#include <iostream>
#include <algorithm>

using namespace std;

int getClassValue(char c) {
    if (c == 'N') return 0;
    if (c == 'Z') return 1;
    if (c == 'Q') return 2;
    if (c == 'R') return 3;
    return -1;
}

char getValueClass(int v) {
    if (v == 0) return 'N';
    if (v == 1) return 'Z';
    if (v == 2) return 'Q';
    if (v == 3) return 'R';
    return ' ';
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    char a, op, b;
    if (cin >> a >> op >> b) {
        int valA = getClassValue(a);
        int valB = getClassValue(b);
        
        int res = max(valA, valB);
        
        if (valA == 0 && valB == 0 && op == '-') {
            res = 1;
        }
        
        cout << getValueClass(res) << "\n";
    }
    
    return 0;
}