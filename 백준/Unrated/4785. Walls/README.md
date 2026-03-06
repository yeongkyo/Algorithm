# [Unrated] Walls - 4785 

[문제 링크](https://www.acmicpc.net/problem/4785) 

### 성능 요약

메모리: 2024 KB, 시간: 40 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 6일 09:23:27

### 문제 설명

<p>There are a number of research stations on a featureless patch of desert, which can be modeled as a Cartesian plane. Each station is located at some point (x,y) where x and y are even integers. For security reasons, sufficiently long and high walls are to be constructed to separate the stations so that no station is visible from any of the other stations. A wall may only be constructed along a North-South or East-West line. A vertical wall may be built at an odd x-coordinate, and a horizontal wall may be built at an odd y-coordinate. Since the stations are located at even valued coordinates, and the walls are built along odd valued coordinates, no wall can ever touch a station. The walls are always long enough to completely separate stations on one side from stations on the other side.</p>

<p>Given a list of stations, you must determine the smallest number of walls that need to be constructed.</p>

### 입력 

 <p>There will be several test cases in the input. Each test case will begin with an integer n (2≤n≤100), which is the number of stations. The next n lines each will contain two integers x and y (0≤x,y≤36), separated by a single space, indicating the (x,y) location of a station. The x and y values are guaranteed to be even. Within a test case, all (x,y) locations will be unique. The last test case will be followed by a line with a single 0.</p>

### 출력 

 <p>For each test case, output a single integer, indicating the smallest number of walls that can prevent the given n stations from seeing each other. That is, a straight line-segment joining any two stations must be intersected by at least one wall. Output no extra spaces, and do not separate answers with blank lines.</p>

