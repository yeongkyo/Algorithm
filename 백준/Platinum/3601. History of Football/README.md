# [Platinum II] History of Football - 3601 

[문제 링크](https://www.acmicpc.net/problem/3601) 

### 성능 요약

메모리: 2020 KB, 시간: 740 ms

### 분류

다이나믹 프로그래밍, 자료 구조, 브루트포스 알고리즘, 정렬, 집합과 맵, 백트래킹, 트리를 사용한 집합과 맵

### 제출 일자

2026년 2월 23일 08:48:31

### 문제 설명

<p>Henry is a historian. He specializes in the history of sports, especially football. Whenever he sees a table of a football tournament, he saves it into his database.</p>

<p>Recently he ran across a web-site with standings of a small tournament. Unfortunately for him, the results of the games were lost, and the only available information was the amount of points gained by each team.</p>

<p>Disappointed by that, he decides to have some mathematical fun and to calculate in how many different ways the games of the championship could have ended. He doesn’t care about the scores of the games, he only cares about the winners.</p>

<p>In that tournament the following rules were applied:</p>

<ul>
	<li>Each team plays against each other team exactly once.</li>
	<li>In case of a tie each team gains 1 point.</li>
	<li>In other case the winner gains 3 points and the loser gains 0 points.</li>
</ul>

<p>For example, if Henry knows that each of 3 teams had got 3 points by the end of the tournament, the answer to his question is that there are two possible tournament tables:</p>

<table class="table table-bordered th-center td-center table-center-30">
	<caption>Possible table number 1</caption>
	<thead>
		<tr>
			<th>Team</th>
			<th>A</th>
			<th>B</th>
			<th>C</th>
			<th>Points</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<th>A</th>
			<td>-</td>
			<td>3</td>
			<td>0</td>
			<td>3</td>
		</tr>
		<tr>
			<th>B</th>
			<td>0</td>
			<td>-</td>
			<td>3</td>
			<td>3</td>
		</tr>
		<tr>
			<th>C</th>
			<td>3</td>
			<td>0</td>
			<td>-</td>
			<td>3</td>
		</tr>
	</tbody>
</table>

<table class="table table-bordered th-center td-center table-center-30">
	<caption>Possible table number 2</caption>
	<thead>
		<tr>
			<th>Team</th>
			<th>A</th>
			<th>B</th>
			<th>C</th>
			<th>Points</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<th>A</th>
			<td>-</td>
			<td>0</td>
			<td>3</td>
			<td>3</td>
		</tr>
		<tr>
			<th>B</th>
			<td>3</td>
			<td>-</td>
			<td>0</td>
			<td>3</td>
		</tr>
		<tr>
			<th>C</th>
			<td>0</td>
			<td>3</td>
			<td>-</td>
			<td>3</td>
		</tr>
	</tbody>
</table>

<p>Help Henry calculate the number of different possible tournament tables (without consideration of the scores of the games).</p>

### 입력 

 <p>Input file contains integer n, the number of teams in the championship (2 ≤ n ≤ 8). The following n lines contain one integer number each — points gained by the teams.</p>

### 출력 

 <p>Output one integer number — the number of possible tournament tables with given total points. It is guaranteed that there is at least one such tournament table.</p>

