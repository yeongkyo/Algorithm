# [Diamond III] String Puzzle - 15338 

[문제 링크](https://www.acmicpc.net/problem/15338) 

### 성능 요약

메모리: 2164 KB, 시간: 104 ms

### 분류

트리, 집합과 맵

### 제출 일자

2026년 3월 27일 06:52:17

### 문제 설명

<p>Amazing Coding Magazine is popular among young programmers for its puzzle solving contests offering catchy digital gadgets as the prizes. The magazine for programmers naturally encourages the readers to solve the puzzles by writing programs. Let’s give it a try!</p>

<p>The puzzle in the latest issue is on deciding some of the letters in a string (the secret string, in what follows) based on a variety of hints. The figure below depicts an example of the given hints.</p>

<p style="text-align:center"><img alt="" src="https://onlinejudgeimages.s3-ap-northeast-1.amazonaws.com/problem/15338/1.png" style="height:138px; width:419px"></p>

<p>The first hint is the number of letters in the secret string. In the example of the figure above, it is nine, and the nine boxes correspond to nine letters. The letter positions (boxes) are numbered starting from 1, from the left to the right.</p>

<p>The hints of the next kind simply tell the letters in the secret string at some specific positions. In the example, the hints tell that the letters in the 3rd, 4th, 7th, and 9th boxes are <code>C</code>, <code>I</code>, <code>C</code>, and <code>P</code>, respectively.</p>

<p>The hints of the final kind are on duplicated substrings in the secret string. The bar immediately below the boxes in the figure is partitioned into some sections corresponding to substrings of the secret string. Each of the sections may be connected by a line running to the left with another bar also showing an extent of a substring. Each of the connected pairs indicates that substrings of the two extents are identical. One of this kind of hints in the example tells that the letters in boxes 8 and 9 are identical to those in boxes 4 and 5, respectively. From this, you can easily deduce that the substring is <code>IP</code>.</p>

<p>Note that, not necessarily all of the identical substring pairs in the secret string are given in the hints; some identical substring pairs may not be mentioned.</p>

<p>Note also that two extents of a pair may overlap each other. In the example, the two-letter substring in boxes 2 and 3 is told to be identical to one in boxes 1 and 2, and these two extents share the box 2.</p>

<p>In this example, you can decide letters at all the positions of the secret string, which are “<code>CCCIPCCIP</code>”. In general, the hints may not be enough to decide all the letters in the secret string.</p>

<p>The answer of the puzzle should be letters at the specified positions of the secret string. When the letter at the position specified cannot be decided with the given hints, the symbol ? should be answered.</p>

### 입력 

 <p>The input consists of a single test case in the following format.</p>

<pre>n a b q
x<sub>1</sub> c<sub>1</sub>
.
.
.
x<sub>a</sub> c<sub>a</sub>
y<sub>1</sub> h<sub>1</sub>
.
.
.
y<sub>b</sub> h<sub>b</sub>
z<sub>1</sub>
.
.
.
z<sub>q</sub></pre>

<p>The first line contains four integers n, a, b, and q. n (1 ≤ n ≤ 10<sup>9</sup>) is the length of the secret string, a (0 ≤ a ≤ 1000) is the number of the hints on letters in specified positions, b (0 ≤ b ≤ 1000) is the number of the hints on duplicated substrings, and q (1 ≤ q ≤ 1000) is the number of positions asked.</p>

<p>The i-th line of the following a lines contains an integer x<sub>i</sub> and an uppercase letter c<sub>i</sub> meaning that the letter at the position x<sub>i</sub> of the secret string is c<sub>i</sub>. These hints are ordered in their positions, i.e., 1 ≤ x<sub>1</sub> < · · · < x<sub>a</sub> ≤ n.</p>

<p>The i-th line of the following b lines contains two integers, y<sub>i</sub> and h<sub>i</sub>. It is guaranteed that they satisfy 2 ≤ y<sub>1</sub> < · · · < y<sub>b</sub> ≤ n and 0 ≤ h<sub>i</sub> < y<sub>i</sub>. When h<sub>i</sub> is not 0, the substring of the secret string starting from the position y<sub>i</sub> with the length y<sub>i+1</sub>−y<sub>i</sub> (or n+1−y<sub>i</sub> when i = b) is identical to the substring with the same length starting from the position h<sub>i</sub>. Lines with h<sub>i</sub> = 0 does not tell any hints except that y<sub>i</sub> in the line indicates the end of the substring specified in the line immediately above.</p>

<p>Each of the following q lines has an integer z<sub>i</sub> (1 ≤ z<sub>i</sub> ≤ n), specifying the position of the letter in the secret string to output.</p>

<p>It is ensured that there exists at least one secret string that matches all the given information. In other words, the given hints have no contradiction.</p>

### 출력 

 <p>The output should be a single line consisting only of q characters. The character at position i of the output should be the letter at position z<sub>i</sub> of the the secret string if it is uniquely determined from the hints, or ? otherwise.</p>

