# [Diamond V] Weeping Fig - 13367 

[문제 링크](https://www.acmicpc.net/problem/13367) 

### 성능 요약

메모리: 3068 KB, 시간: 708 ms

### 분류

그래프 이론, 스토어–바그너

### 제출 일자

2026년 3월 27일 07:13:58

### 문제 설명

<p>“Roots grow rapidly, invading gardens, growing under and lifting sidewalks, patios, and driveways.” - United States Forest Service</p>

<p>Because of that, it is not recommend to be planted in residential household. However, you happen to have one in your backyard, and your fig is very nice — it does not invade any part of your household at all. It just expands into empty space around the house.</p>

<p>Your friend are constructing a new house, and want to have a weeping fig in his backyard too. But he is scared that traditional weeping fig could destroy his home, so he wants your weeping fig (which is now known among your friends to be a “nice weeping fig.”) Since the tree grows rapidly, you decided to just cut parts of your tree for your friend.</p>

<p>Even though just only one root is enough for your friend to have the same nice tree as you, it is not easyseparating one root since it has many connections to other roots. Specifically, each root has many connections to some other roots, with no pattern. Each root connection also vary in diameter, thus the effort required for separating is not constant.</p>

<p>Given the number of roots, its connection and connection separation effort, please find the minimal effort required to separate your weeping fig into two entities.</p>

### 입력 

 <p>First line has the number of test cases, T (1 ≤ T ≤ 20)</p>

<p>For each test case, the first line contains two integers N M specifying the number of roots and number of connection (1 ≤ N ≤ 500, N-1 ≤ M ≤ (N<sup>2</sup> -N)/2)</p>

<p>Following M lines contain three integer a<sub>i</sub> b<sub>i</sub> w<sub>i</sub> (1 ≤ a<sub>i</sub> , b<sub>i</sub> ≤ N, 1 ≤ w<sub>i</sub> ≤ 1000) indicating that there is a root connection from root a<sub>i</sub> to root b<sub>i</sub> that requires effort w<sub>i</sub> to separate. It is guaranteed that the input will form single weeping fig tree only.</p>

### 출력 

 <p>For each test case, output the minimal effort required on its own line.</p>

