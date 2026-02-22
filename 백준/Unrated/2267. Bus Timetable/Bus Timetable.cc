#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>

using namespace std;

// 개별 버스 운행 일정을 담는 구조체
struct Service {
    int depart_time;
    int route_id;
    map<int, int> stop_times; // 정류장 번호 -> 도착 시간(분)
    
    // 정렬 기준: 1순위 차고지 출발 시간, 2순위 노선 번호
    bool operator<(const Service& other) const {
        if (depart_time != other.depart_time)
            return depart_time < other.depart_time;
        return route_id < other.route_id;
    }
};

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string scenario;
    bool is_first_scenario = true;

    // '#'이 입력될 때까지 시나리오 반복
    while (cin >> scenario && scenario != "#") {
        if (!is_first_scenario) {
            cout << "\n"; // 시나리오 사이 빈 줄 출력
        }
        is_first_scenario = false;

        int h, m;
        int route_id = 1;
        vector<Service> all_services;
        set<int> all_stops; // 모든 정류장 번호를 오름차순으로 관리하기 위해 set 사용

        // 0 0 이 입력될 때까지 노선(Route) 정보 읽기
        while (cin >> h >> m) {
            if (h == 0 && m == 0) break;
            
            int interval;
            cin >> interval;

            int start_time = h * 60 + m; // 출발 시간을 분 단위로 변환
            vector<pair<int, int>> stops_info;
            int travel_time, stop_num;
            int cumulative_time = 0;
            
            // 각 정류장까지의 소요 시간 및 정류장 번호 읽기
            while (cin >> travel_time >> stop_num) {
                cumulative_time += travel_time;
                if (stop_num != 0) {
                    stops_info.push_back({stop_num, cumulative_time});
                    all_stops.insert(stop_num);
                } else {
                    // 정류장 번호가 0이면 차고지 복귀이므로 해당 노선 입력 종료
                    break; 
                }
            }
            
            int total_duration = cumulative_time;
            int current_depart = start_time;

            // 차고지 복귀 시간이 자정(1440분) 이하일 때만 운행 스케줄 생성
            while (current_depart + total_duration <= 1440) {
                // 오전 6시(360분) 이후 출발인지 확인
                if (current_depart >= 360) {
                    Service s;
                    s.depart_time = current_depart;
                    s.route_id = route_id;
                    for (const auto& p : stops_info) {
                        // 정류장 도착 시간 = 차고지 출발 시간 + 누적 소요 시간
                        s.stop_times[p.first] = current_depart + p.second;
                    }
                    all_services.push_back(s);
                }
                
                // 배차 간격이 0인 경우 방어 코드 (무한 루프 방지)
                if (interval == 0) break;
                current_depart += interval;
            }
            
            route_id++;
        }

        // 전체 버스 운행 스케줄을 조건에 맞게 정렬
        sort(all_services.begin(), all_services.end());

        // 1. 시나리오 이름 출력
        cout << scenario << "\n";
        
        // 2. 시간표 헤더(노선 번호) 출력
        cout << "  |";
        for (const auto& s : all_services) {
            cout << setw(5) << s.route_id << "|";
        }
        cout << "\n";

        // 3. 각 정류장별 시간표 데이터 출력
        for (int stop_id : all_stops) {
            cout << setw(2) << stop_id << "|";
            for (const auto& s : all_services) {
                // 해당 버스가 이 정류장에 정차하는 경우
                if (s.stop_times.count(stop_id)) {
                    int t = s.stop_times.at(stop_id);
                    // hh:mm 포맷으로 출력 (한 자리 숫자일 경우 앞에 0을 채움)
                    cout << setfill('0') << setw(2) << (t / 60) << ":" 
                         << setw(2) << (t % 60) << setfill(' ') << "|";
                } else {
                    // 정차하지 않으면 공백 5칸으로 비워둠
                    cout << "     |";
                }
            }
            cout << "\n";
        }
    }

    return 0;
}