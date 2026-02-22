#include <iostream>
#include <string>
#include <vector>

using namespace std;

// AST를 구성할 노드 구조체
struct Node {
    char type; // 'v': 변수, '+': 덧셈, '*': 곱셈
    char var;
    vector<Node*> children;
    long long total_terms; // 서브트리에서 파생될 수 있는 총 항의 개수
};

// 괄호와 기호를 읽으며 AST를 재귀적으로 구성하는 함수
Node* parse_ast(int& pos, const string& s) {
    Node* node = new Node();
    // 영소문자인 경우 변수 노드로 처리
    if (s[pos] >= 'a' && s[pos] <= 'z') {
        node->type = 'v';
        node->var = s[pos++];
        node->total_terms = 1; // 단일 변수는 무조건 항 1개
        return node;
    }
    
    pos++; // '(' 건너뛰기
    node->type = s[pos++];
    // 덧셈 노드는 0부터 합산, 곱셈 노드는 1부터 곱산
    node->total_terms = (node->type == '+') ? 0 : 1; 
    
    // 닫는 괄호가 나올 때까지 자식 노드 파싱
    while (s[pos] != ')') {
        Node* child = parse_ast(pos, s);
        node->children.push_back(child);
        if (node->type == '+') {
            node->total_terms += child->total_terms;
        } else {
            node->total_terms *= child->total_terms;
        }
    }
    pos++; // ')' 건너뛰기
    return node;
}

// AST를 순회하며 정확히 i번째 항의 문자열을 생성하는 함수
string get_term(Node* node, long long i) {
    if (node->type == 'v') {
        return string(1, node->var);
    } 
    else if (node->type == '+') {
        // 덧셈: i번째 항이 속하는 자식 노드를 찾음
        for (Node* child : node->children) {
            if (i < child->total_terms) {
                return get_term(child, i);
            }
            i -= child->total_terms;
        }
    } 
    else if (node->type == '*') {
        // 곱셈: 1차원 인덱스 i를 다차원 조합 인덱스로 변환 (Cartesian Product)
        string res = "";
        vector<long long> child_indices(node->children.size());
        long long current_i = i;
        
        // 오른쪽 자식부터 나머지 연산을 통해 각 자식별 인덱스 추출
        for (int j = (int)node->children.size() - 1; j >= 0; --j) {
            long long ways = node->children[j]->total_terms;
            child_indices[j] = current_i % ways;
            current_i /= ways;
        }
        
        // 추출된 인덱스대로 왼쪽 자식부터 순차적으로 문자열 결합
        for (int j = 0; j < (int)node->children.size(); ++j) {
            res += get_term(node->children[j], child_indices[j]);
        }
        return res;
    }
    return "";
}

// 메모리 누수 방지를 위한 트리 삭제 함수
void free_ast(Node* node) {
    for (Node* child : node->children) {
        free_ast(child);
    }
    delete node;
}

int main() {
    // 입출력 속도 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    char c;
    while (cin >> c) {
        string F = string(1, c);
        int balance = 0;
        
        // 공백을 무시하며 수식 F만 깔끔하게 추출
        if (c == '(') {
            balance = 1;
            while (balance > 0 && cin >> c) {
                F += c;
                if (c == '(') balance++;
                else if (c == ')') balance--;
            }
        }
        
        int pos = 0;
        Node* root = parse_ast(pos, F);
        long long total_terms = root->total_terms;
        
        long long current_idx = 0;
        long long k;
        
        // 명령어 0이 나올 때까지 반복
        while (cin >> k && k != 0) {
            if (k > 0) {
                for (long long i = 0; i < k; ++i) {
                    string t = get_term(root, current_idx);
                    
                    // 정규 폼(NF) 포맷팅
                    if (t.length() == 1) cout << t << "\n";
                    else cout << "(*" << t << ")\n";
                    
                    current_idx = (current_idx + 1) % total_terms;
                }
            } else {
                // k가 음수일 경우 건너뛰기만 수행
                long long skip = -k;
                current_idx = (current_idx + skip) % total_terms;
            }
        }
        
        free_ast(root);
    }
    
    return 0;
}