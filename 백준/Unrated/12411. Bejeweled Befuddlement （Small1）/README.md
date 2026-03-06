# [Unrated] Bejeweled Befuddlement (Small1) - 12411 

[문제 링크](https://www.acmicpc.net/problem/12411) 

### 성능 요약

메모리: 2028 KB, 시간: 0 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 6일 09:03:35

### 문제 설명

<p>John has been playing a lot Bejeweled lately, but he's becoming frustrated at his inability to get high scores. He would like a program to help him play, and you can help!</p>

<p>In Bejeweled, you have a grid of gems of various colors, such as the one below:</p>

<pre>BGR
GRR
BGY
</pre>

<p>You make a move by switching two gems that are horizontally or vertically connected to each other. So, for example, you could switch the first two gems in the second row, changing the board to this:</p>

<pre>BGR
RGR
BGY
</pre>

<p>If when you do this, any three gems in a row are the same color, they all disappear at once. So in this example, the three green gems would disappear.</p>

<pre>B.R
R.R
B.Y
</pre>

<p>It can get a little more complex, though. Additionally, multiple sets of three (or more) can be formed at once. For example:</p>

<pre>RGOY
RGOG
OOYO
RGOG
YBOB
</pre>

<p>In this case, swapping the right-most orange in the middle row with the yellow next to it...</p>

<pre>RGOY
RGOG
OOOY
RGOG
YBOB
</pre>

<p>... makes a row of 3 and a column of 5, which all disappear...</p>

<pre>RG.Y
RG.G
...Y
RG.G
YB.B
</pre>

<p>... and then, any gems that are above an open space fall down to fill the empty spaces...</p>

<pre>...Y
RG.G
RG.Y
RG.G
YB.B
</pre>

<p>... and then any newly formed rows of three or more will also disappear:</p>

<pre>...Y
...G
...Y
...G
YB.B
</pre>

<p>This continues until there are no more rows or columns of three or more similar colored gems formed. So, more generally, the process is as follows:</p>

<ul>
	<li>Switch two gems</li>
	<li>While there are any rows or columns of three or more gems of the same color
	<ul>
		<li>Remove those gems</li>
		<li>Move any gems that are above empty spaces down until there are no gems above empty spaces</li>
	</ul>
	</li>
</ul>

<p>Note: In the actual game, gems fall in from the top to replace removed gems. You do not need to take this into account.</p>

### 입력 

 <p>The first line of input will contain the number of test cases, <strong>T</strong>. Each test case will start with a line containing two integers, <strong>N</strong> and <strong>M</strong>, separated by a space. The next <strong>N</strong> lines will each contain exactly <strong>M</strong> uppercase letters, each representing one gem. Gems of the same letter are considered to be of the same color.</p>

<h3>Limits</h3>

<ul>
	<li>All gem characters will be upper-case letters.</li>
	<li>The input grid will not contain a row or a column of 3 or more gems of the same color appearing consecutively.</li>
	<li>1 ≤ <strong>T</strong> ≤ 20.</li>
	<li>3 ≤ <strong>number of rows</strong> ≤ 10.</li>
	<li>3 ≤ <strong>number of columns</strong> ≤ 10.</li>
</ul>

### 출력 

 <p>For each test case, output a line of the form "Case #<strong>C</strong>: <strong>D</strong>", where <strong>C</strong> is the number of the test case, starting from 1, and <strong>D</strong> is the maximum number of gems that can be removed as a result of a single swap.</p>

