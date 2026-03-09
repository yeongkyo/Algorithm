# [Unrated] Polyomino Powers - 5028 

[문제 링크](https://www.acmicpc.net/problem/5028) 

### 성능 요약

메모리: 2024 KB, 시간: 244 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 9일 13:28:16

### 문제 설명

<p>A polyomino is a polyform with the square as its base form. It is a connected shape formed as the union of one or more identical squares in distinct locations on the plane, taken from the regular square tiling, such that every square can be connected to every other square through a sequence of shared edges (i.e., shapes connected only through shared corners of squares are not permitted). The most well-known polyominos are the seven tetrominos made out of four squares (see figure), famous from the TetrisⓇ game, and of course the single domino consisting of two squares from the game with the same name. Some polyomino can be obtained by gluing several copies of the same smaller polyomino translated (but not rotated or mirrored) to different locations in the plane. We call those polyomino powers.</p>

### 입력 

 <p>One line with two positive integers h, w ≤ 10.</p>

<p>Next follows an h × w matrix of characters ’.’ or ’X’, the ’X”s describing a polyomino and ’.’ space.</p>

### 출력 

 <p>A k-power with 2 ≤ k ≤ 5 copies of a smaller polyomino:</p>

<p>Output a h × w matrix on the same format as the input with the ’X”s replaced by the numbers 1 through k in any order identifying the factor pieces.</p>

<p>Furthermore, if multiple solutions exist, any will do. Otherwise, output ’No solution’ if no solution exists.</p>

