# [Gold II] Lucky Controller - 7491 

[문제 링크](https://www.acmicpc.net/problem/7491) 

### 성능 요약

메모리: 2084 KB, 시간: 0 ms

### 분류

다이나믹 프로그래밍, 수학

### 제출 일자

2026년 3월 6일 08:37:03

### 문제 설명

<p>Egor works as a conductor in the bus. Each day he is given a pack of tickets which he then sells. Recently he has become interested about how many tickets in the pack are lucky. He thinks that the more tickets are lucky the luckier day he will have. Now he wants to find out how lucky for him the next day is going to be. Each ticket number consist of <strong>n</strong> digits. The ticket is considered to be lucky if the sum of the first <strong>n</strong>/<strong>2</strong> digits equals to the sum of the last <strong>n</strong>/<strong>2</strong> digits. Egor knows that the numbers in the pack that he will be given can start with equal probability from any number in the interval from <strong>a</strong> to <strong>b</strong> inclusively. The pack holds <strong>k</strong> tickets. The numbers of the tickets are consecutive. Help Egor to find an expected quantity of lucky tickets in the pack.</p>

### 입력 

 <p>The input file consists of a single line with three integers <strong>a</strong>, <strong>b</strong> and <strong>k</strong> (<strong>0</strong> ≤ <strong>a</strong> ≤ <strong>b</strong> < <strong>10</strong><strong><sup>12</sup></strong>, <strong>1</strong> ≤ <strong>k</strong> ≤ <strong>100000</strong>). Integers <strong>a</strong> and <strong>b</strong> consist of same amount of digits, and this amount equals to the amount of digits in the number of each ticket. They may start with zeroes. The amount of digits in <strong>a</strong> and <strong>b</strong> is always even.</p>

### 출력 

 <p>Output the expected quantity of lucky tickets in the pack in the form of irreducible fraction. In case the result is an integer – no slash should appear in the output.</p>

