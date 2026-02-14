# [Unrated] Tower - 17630 

[문제 링크](https://www.acmicpc.net/problem/17630) 

### 성능 요약

메모리: 2020 KB, 시간: 8 ms

### 분류

분류 없음

### 제출 일자

2026년 2월 15일 05:56:41

### 문제 설명

<p>Farmhand Jernej gets bored in the evenings, thus he invented a simple game. He wants to build a tower from numbers on pieces of paper. He starts with a piece of paper and writes 1 on it.</p>

<p>Jernej can write another number on a piece of paper and place it on top of the tower. The new value on the top of the tower must be a valid sum of numbers on consecutive papers comprising the tower. Let's say there are currently <em>n</em> pieces of paper comprising the tower. He makes a sum of numbers in the tower within arbitrary positions [<em>l</em>, <em>u</em>], where 1 ≤ <em>l</em> ≤ <em>u</em> ≤ <em>n</em> and puts the sum on top.</p>

<p>Jernej wants to produce <em>T</em> towers with desired numbers on top. Help him find out the required steps. He also asks you to minimize the number of those steps.</p>

### 입력 

 <p>In the first line of input, you will be given a positive integer <em>T</em> (number of different towers Jernej wants to produce).</p>

<p>In each of the next <em>T</em> lines, there will be one positive integer <em>q</em> , the value Jernej wants to produce at the end. All games are independent.</p>

### 출력 

 <p>For each integer <em>q</em>:</p>

<ul>
	<li>print one line with number <em>s</em> (0 ≤ <em>s</em> ≤ 1000) - required steps to produce the desired value.</li>
	<li>In the next <em>s</em> lines, there should be 2 space-separated positive integers <em>l</em> and <em>u</em>, range bounds to produce the desired value.</li>
</ul>

