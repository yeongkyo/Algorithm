# [Diamond II] Colorgraph - 16544 

[문제 링크](https://www.acmicpc.net/problem/16544) 

### 성능 요약

메모리: 2828 KB, 시간: 24 ms

### 분류

자료 구조, 그래프 이론, 많은 조건 분기, 분리 집합, 스토어–바그너

### 제출 일자

2026년 3월 27일 11:06:52

### 문제 설명

<p>Lora is playing an online puzzle game. She receives an undirected graph with N vertices, numbered from 1 to N. The graph is such that between every two distinct vertices there is an edge, colored either in blue, or in red. We will call the graph redconnected, if from any vertex we can reach any other vertex, using only red edges. Similarly, a graph is blue-connected if from any vertex we can reach any other vertex, using only blue edges. We now define the state of the graph as a pair of numbers (A, B), such that:</p>

<ul>
	<li>A=1 if the graph is red-connected, and A=0 if it is not</li>
	<li>B=1 if the graph is blue-connected, and B=0 if it is not</li>
</ul>

<p>So for example, state (1, 0) describes a graph that is red-connected, but not blueconnected.</p>

<p>With a single click on a given edge, Lora can change its color (from blue to red or from red to blue). The goal of the game is, given an initial graph and a desired final state, to change the initial graph to one that is in the final state, using minimum amount of clicks (see the sample test cases for more information). Your task is to help Lora by writing a program colorgraph that finds the minimum amount of clicks needed to solve the problem. </p>

### 입력 

 <p>On the first line of the standard input is the positive integer N – the amount of vertices in the graph. After that, N lines follow, each with N space-separated numbers, describing the colors of the edges. Let us denote the j-th number on the i-th of those lines by G<sub>ij</sub>. If G<sub>ij</sub>=0, then the edge between i and j is red, and if G<sub>ij</sub>=1, then the edge between i and j is blue. It is guaranteed that G<sub>ij</sub>=G<sub>ji</sub>. For i=j, the value of G<sub>ij</sub> is irrelevant, since the graph contains no loops. On the last line are two space-separated numbers – A and B, describing the desired final state of the graph.</p>

### 출력 

 <p>If it is impossible to transform the graph into one with the desired state, you must print -1 on a single line of the standard output. In all other cases, on the first line of the standard output you must print a single integer K – the minimum amount of clicks Lora needs to make in order to transform the graph into the desired state. On each of the next K lines you must print a pair of numbers – the two endpoints of the edge that Lora should click on. If there is more than one solution, print any of them. The order, in which the edges are printed, as well as the order of the two endpoints, is irrelevant.</p>

