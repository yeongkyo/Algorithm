# [Ruby IV] Fibonacci representations - 16029 

[문제 링크](https://www.acmicpc.net/problem/16029) 

### 성능 요약

메모리: 173900 KB, 시간: 320 ms

### 분류

수학, 다이나믹 프로그래밍, 자료 구조, 트리, 세그먼트 트리, 집합과 맵, 트리를 사용한 집합과 맵, 레드-블랙 트리

### 제출 일자

2026년 3월 27일 13:56:26

### 문제 설명

<p>Let us define the sequence of Fibonacci numbers as:</p>

<ul>
	<li>F<sub>1</sub> = 1</li>
	<li>F<sub>2</sub> = 2</li>
	<li>F<sub>n</sub> = F<sub>n−1</sub> + F<sub>n−2</sub> for n ≥ 3</li>
</ul>

<p>The first few elements of the sequence are 1, 2, 3, 5, 8, 13, 21, . . .</p>

<p>For a positive integer p, let X(p) denote the number of different ways of expressing p as a sum of different Fibonacci numbers. Two ways are considered different if there is a Fibonacci number that exists in exactly one of them.</p>

<p>You are given a sequence of n positive integers a<sub>1</sub>, a<sub>2</sub>, . . . , a<sub>n</sub>. For a non-empty prefix a<sub>1</sub>, a<sub>2</sub>, . . . , a<sub>k</sub>, we define p<sub>k</sub> = F<sub>a<sub>1</sub></sub> + F<sub>a<sub>2</sub></sub> + . . . + F<sub>a<sub>k</sub></sub>. Your task is to find the values X(p<sub>k</sub>) modulo 10<sup>9</sup> + 7, for all k = 1, . . . , n.</p>

### 입력 

 <p>The first line of the standard input contains an integer n (1 ≤ n ≤ 100 000). The second line contains n space-separated integers a<sub>1</sub>, a<sub>2</sub>, . . . , a<sub>n</sub> (1 ≤ a<sub>i</sub> ≤ 10<sup>9</sup>).</p>

### 출력 

 <p>The standard output should contain n lines. In the k-th line, print the value X(p<sub>k</sub>) modulo (10<sup>9</sup> + 7).</p>

