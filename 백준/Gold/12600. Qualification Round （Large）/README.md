# [Gold I] Qualification Round (Large) - 12600 

[문제 링크](https://www.acmicpc.net/problem/12600) 

### 성능 요약

메모리: 2020 KB, 시간: 0 ms

### 분류

애드 혹, 이분 탐색, 그리디 알고리즘, 수학, 매개 변수 탐색

### 제출 일자

2026년 3월 9일 13:16:19

### 문제 설명

<p>You've just advanced from the Qualification Round of Google Code Jam Africa 2010, and you want to know how many of your fellow contestants advanced with you. To give yourself a challenge, you've decided only to look at how many people solved each problem.</p>

<p>The Qualification Round consisted of <strong>P</strong> problems; the i<sup>th</sup> problem was fully solved by <strong>S<sub>i</sub></strong>contestants. Contestants had to solve <strong>C</strong> problems in order to advance to the next round. Your job is to figure out, using only that information, the maximum number of contestants who could have advanced.</p>

### 입력 

 <p>The first line of the input gives the number of test cases, <strong>T</strong>. T lines follow. Each will consist only of space-separated integers: first <strong>P</strong>, then <strong>C</strong>, then P integers <strong>S<sub>0</sub>...S<sub>P-1</sub></strong>.</p>

<h3>Limits</h3>

<ul>
	<li>1 ≤ T ≤ 100</li>
	<li>1 ≤ C ≤ P</li>
	<li>1 ≤ P ≤ 60</li>
	<li>0 ≤ S<sub>i</sub> ≤ 10<sup>17</sup></li>
</ul>

### 출력 

 <p>For each test case, output one line containing "Case #x: y", where x is the case number (starting from 1) and y is the maximum number of contestants who could have advanced (in other words, the maximum number of contestants who could have solved at least <strong>C</strong>problems).</p>

