#include <iostream>
#include <string>
#include <cctype>
#include <vector>

using namespace std;

// 알파벳별 기본 점수를 반환하는 함수
int get_letter_value(char c) {
    switch(c) {
        case 'Q': case 'Z': return 10;
        case 'J': case 'X': return 8;
        case 'K': return 5;
        case 'F': case 'H': case 'V': case 'W': case 'Y': return 4;
        case 'B': case 'C': case 'M': case 'P': return 3;
        case 'D': case 'G': return 2;
        default: return 1; // A, E, I, L, N, O, R, S, T, U
    }
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // 15x15 보드의 좌측 상단(8x8) 사분면 보너스 정보만 하드코딩
    string board[8][8] = {
        {"3W", "", "", "2L", "", "", "", "3W"},
        {"", "2W", "", "", "", "3L", "", ""},
        {"", "", "2W", "", "", "", "2L", ""},
        {"2L", "", "", "2W", "", "", "", "2L"},
        {"", "", "", "", "2W", "", "", ""},
        {"", "3L", "", "", "", "3L", "", ""},
        {"", "", "2L", "", "", "", "2L", ""},
        {"3W", "", "", "2L", "", "", "", "2W"}
    };

    string pos_str;
    // "#"이 입력될 때까지 테스트 케이스 반복
    while (cin >> pos_str && pos_str != "#") {
        string word;
        cin >> word;

        bool horizontal = false;
        int r = 0, c = 0;

        // 문자열 파싱: 숫자로 시작하면 가로(Row 우선), 알파벳이면 세로(Col 우선)
        if (isdigit(pos_str[0])) {
            horizontal = true;
            int i = 0;
            while (i < pos_str.length() && isdigit(pos_str[i])) {
                r = r * 10 + (pos_str[i] - '0');
                i++;
            }
            c = pos_str[i] - 'A' + 1;
        } else {
            horizontal = false;
            c = pos_str[0] - 'A' + 1;
            int i = 1;
            while (i < pos_str.length() && isdigit(pos_str[i])) {
                r = r * 10 + (pos_str[i] - '0');
                i++;
            }
        }

        // 1-based 인덱스를 0-based 인덱스로 변환
        r--; c--; 

        int word_mult = 1;
        int letter_sum = 0;

        // 단어의 각 알파벳을 보드에 올려놓으며 점수 계산
        for (int i = 0; i < word.length(); ++i) {
            int curr_r = r + (horizontal ? 0 : i);
            int curr_c = c + (horizontal ? i : 0);

            // 대칭성을 이용해 좌측 상단(0~7) 인덱스로 변환
            int mapped_r = curr_r > 7 ? 14 - curr_r : curr_r;
            int mapped_c = curr_c > 7 ? 14 - curr_c : curr_c;

            string bonus = board[mapped_r][mapped_c];
            int val = get_letter_value(word[i]);

            if (bonus == "2L") {
                val *= 2;
            } else if (bonus == "3L") {
                val *= 3;
            } else if (bonus == "2W") {
                word_mult *= 2;
            } else if (bonus == "3W") {
                word_mult *= 3;
            }

            letter_sum += val;
        }

        // 최종 점수 출력 (출력 형식에 맞춤)
        cout << pos_str << " " << word << " " << letter_sum * word_mult << "\n";
    }

    return 0;
}