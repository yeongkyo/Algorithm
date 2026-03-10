# [Gold IV] Safe Squares (Large) - 14309 

[문제 링크](https://www.acmicpc.net/problem/14309) 

### 성능 요약

메모리: 45996 KB, 시간: 192 ms

### 분류

다이나믹 프로그래밍, 누적 합, 매개 변수 탐색

### 제출 일자

2026년 3월 10일 15:25:40

### 문제 설명

<p>Codejamon trainers are actively looking for monsters, but if you are not a trainer, these monsters could be really dangerous for you. You might want to find safe places that do not have any monsters!</p>

<p>Consider our world as a grid, and some of the cells have been occupied by monsters. We define a safe square as a grid-aligned D × D square of grid cells (with D ≥ 1) that does not contain any monsters. Your task is to find out how many safe squares (of any size) we have in the entire world.</p>

<ul>
</ul>

### 입력 

 <p>The first line of the input gives the number of test cases, T. T test cases follow. Each test case starts with a line with three integers, R, C, and K. The grid has R rows and C columns, and contains K monsters. K more lines follow; each contains two integers R<sub>i</sub> and C<sub>i</sub>, indicating the row and column that the i-th monster is in. (Rows are numbered from top to bottom, starting from 0; columns are numbered from left to right, starting from 0.)</p>

### 출력 

 <p>For each test case, output one line containing <code>Case #x: y</code>, where <code>x</code> is the test case number (starting from 1) and <code>y</code> is the the total number of safe zones for this test case.</p>

