# [Unrated] Puzzle - 6347 

[문제 링크](https://www.acmicpc.net/problem/6347) 

### 성능 요약

메모리: 2024 KB, 시간: 0 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 6일 09:27:51

### 문제 설명

<p>Little Barborka has just started to learn how to solve a picture puzzle. She has started with a small one containing 15 pieces. Her daddy tries to solve the puzzle too. To make it a little bit harder for himself, he has turned all puzzle pieces upside down so that he cannot see pictures on the pieces. Now he is looking for a solution of the puzzle. Normally the solution should exist but he is not sure whether Barborka has not replaced some pieces of the puzzle by pieces of another similar puzzle. Help him and write a program which reads a description of a set of puzzle pieces and decides whether it is possible to assembly the pieces into a rectangle with given side lengths or not.</p>

### 입력 

 <p>The input consists of blocks of lines. Each block except the last describes one puzzle problem. In the first line of the block there are integers n and m, 0 < n, m ≤ 6 separated by one space. The integers n, m indicate the number of rows and columns in the puzzle respectively. The description of individual puzzle pieces is in the following n * m lines of the block. Each piece is a rectangle 3 centimeters wide and 4 centimeters high with possible juts or cavities in the middle of its sides. For each side of a puzzle piece just one of the following possibilities is true (see picture):</p>

<ul>
	<li>there is no jut or cavity on the side, i.e., the side is flat - such sides can be used only on edges of the final picture when assembling the puzzle,</li>
	<li>there is one jut in the middle of the side,</li>
	<li>there is one cavity in the middle of the side.</li>
</ul>

<p style="text-align: center;"><img alt="" src="https://upload.acmicpc.net/d9459ba9-37bd-426e-9935-33bf23f82b60/-/preview/" style="width: 208px; height: 155px;"></p>

<p>As is usual, two pieces can be placed side by side only if one has a jut and the other has a cavity on corresponding sides. We will denote the flat sides by F, the sides with juts by O and the sides with cavities by I. Each piece is described by four letters characterizing its top, right, bottom, and left side. To make the task easier the pieces can be used only as they are described i.e. they cannot be turned.</p>

<p>After each block there is an empty line. The last block consists of just one line containing 0 0, i.e. two zeros separated by one space.</p>

### 출력 

 <p>The output contains the lines corresponding to the blocks in the input. A line contains <code>YES</code> if the corresponding block in the input describes a puzzle that can be correctly assembled. Otherwise it contains <code>NO</code>. There is no line in the output corresponding to the last "null" block of the input.</p>

