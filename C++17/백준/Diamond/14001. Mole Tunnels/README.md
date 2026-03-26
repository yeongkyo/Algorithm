# [Diamond III] Mole Tunnels - 14001 

[문제 링크](https://www.acmicpc.net/problem/14001) 

### 성능 요약

메모리: 4764 KB, 시간: 60 ms

### 분류

다이나믹 프로그래밍, 그래프 이론, 그리디 알고리즘, 트리, 최대 유량, 최소 비용 최대 유량

### 제출 일자

2026년 3월 27일 06:47:27

### 문제 설명

<p>Moles create tunnels for traveling between their holes. In this problem we investigate one tunnel system that was built by moles. It consists of n holes and n − 1 tunnels connecting them. Let us number all holes from 1 to n. Then for all i > 1, a hole number i is connected by a tunnel to the hole number ⌊ i /2 ⌋. Tunnels are bidirectional. For each hole i we know the amount of food c<sub>i</sub> in that hole. It means that there is enough food for exactly c<sub>i</sub> moles in that hole.</p>

<p>There are m moles living in the tunnel system. For each mole i you are given an integer p<sub>i</sub> — the hole number where the mole i is currently sleeping. In the morning, the first k moles wake up and want to eat, while m − k others are sleeping. Each of k woken up moles chooses some hole and crawls to it. They are quite smart, so they want to minimize the total distance traveled. The distance traveled by one mole is the total number of tunnels which it uses to get from one hole to another. The first k moles who woke up want to move in such a way, so that there is enough food for them at the holes they choose to crawl to. It means that in the hole i there are no more than c<sub>i</sub> woken up moles after all their movements are done.</p>

<p>You must find the minimum total distance for all k from 1 to m. It is guaranteed that there always exists a way for all moles to eat.</p>

### 입력 

 <p>The first line contains two integers n and m (1 ≤ n, m ≤ 10<sup>5</sup> ) — the number of holes and moles. The second line contains n integers c<sub>i</sub> (0 ≤ c<sub>i</sub> ≤ m) — the amount of food in the hole i. The third line contains m integers p<sub>i</sub> (1 ≤ p<sub>i</sub> ≤ n) — the starting positions of the moles.</p>

### 출력 

 <p>You must print m numbers. The k-th number is the minimum total distance the first k moles need to travel if they woke up first.</p>

