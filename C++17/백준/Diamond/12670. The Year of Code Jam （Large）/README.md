# [Diamond III] The Year of Code Jam (Large) - 12670 

[문제 링크](https://www.acmicpc.net/problem/12670) 

### 성능 요약

메모리: 3208 KB, 시간: 72 ms

### 분류

그래프 이론, 최대 유량, 이분 그래프, 최대 유량 최소 컷 정리

### 제출 일자

2026년 3월 26일 21:05:22

### 문제 설명

<p>The year 2008 will be known as a year of change and transition, the start of a new era: we're talking, of course, about the new Google Code Jam format. The introduction of this contest has jammed so many great programming competitions together in a single year that people have started calling it <em>The Year of Code Jam</em>.<br>
Sphinny, a passionate contestant, is looking at her calendar of the year and discovering that a great number of programming contests has been scheduled. She has marked every day of the year on the calendar in one of the three ways:</p>

<ul>
	<li>White: She will not participate in a contest on this day. Either no contests are scheduled, or she has more important things to do (surely there are other good things in life!).</li>
	<li>Blue: She will definitely participate in a contest on this day.</li>
	<li>Question mark: There is a contest scheduled, but she has not decided yet whether she will participate.</li>
</ul>

<p>Note: To simplify the problem, we'll assume that there is no concept of qualification: you don't have to participate in one contest to be eligible for another.</p>

<p>Being in a world that is somewhat different from ours, Sphinny's calendar has some features we must mention: It has <strong>N</strong> months, and each month has exactly <strong>M</strong> days.</p>

<p>The picture below depicts a calendar with 5 months, 8 days in each month, 15 blue days, and 5 question marks.</p>

<p><img alt="" src="https://onlinejudgeimages.s3.amazonaws.com/problem/12669/images-12.png"></p>

<p>Looking at her beautiful calendar, Sphinny has decided that each day has up to 4 <strong>neighbors</strong> in the year: The previous day in the same month, the next day in the same month, the same day in the previous month, and the same day in the next month.</p>

<p>Sphinny wants to maximize her happiness from these contests, and she estimates the effect of the contests on her happiness as a summation of values for all the blue days. For each blue day, the value is computed as follows:</p>

<ul>
	<li>The initial value is 4.</li>
	<li>For each blue neighbour the day has, decrease the value by 1.</li>
</ul>

<p>You may think that Sphinny likes the contests, but participating on two consecutive days makes her a little tired. And for aesthetic reasons, participating on the same day in two consecutive months is also not so great.</p>

<p>Sphinny wants to plan her year now, and decide for every day with a question mark whether it should be white or blue. Her goal is simply to maximize the happiness value.</p>

<p>The following picture shows a solution for the example above. By changing two question marks to blue days, and the other three to white days, she can achieve a happiness value of 42.</p>

<p><img alt="" src="https://onlinejudgeimages.s3.amazonaws.com/problem/12669/images-13.png" style="height:204px; width:293px"></p>

### 입력 

 <p>The first line in the input file contains the number of cases <strong>T</strong>. This is followed by <strong>T</strong> cases in the following format.<br>
The first line is of the form "<strong>N M</strong>", where <strong>N</strong> and <strong>M</strong> are two numbers giving the number of months and the number of days per month.<br>
The next <strong>N</strong> lines each contain a string of length <strong>M</strong>. The <strong>j</strong>-th character in the <strong>i</strong>-th string is one of {'#', '.', '?'}, which gives the status of the j-th day in the i-th month. '#' indicates a blue day, '.' indicates a white day, and '?' indicates a day with a question mark.</p>

<p>Limits</p>

<ul>
	<li>1 ≤ <strong>T</strong> ≤ 100.</li>
	<li>1 ≤ <strong>M</strong>, <strong>N</strong> ≤ 50.</li>
</ul>

### 출력 

 <p>For each input case, you should output a line in the format:</p>

<pre>Case #X: Y
</pre>

<p>where <strong>X</strong> is the 1-based case number, and <strong>Y</strong> is the maximum happiness value.</p>

