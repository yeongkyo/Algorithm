# [Unrated] GIGA Universe Cup - 15542 

[문제 링크](https://www.acmicpc.net/problem/15542) 

### 성능 요약

메모리: 2216 KB, 시간: 44 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 13일 07:09:10

### 문제 설명

<p>Following FIFA World Cup, a larger competition called "GIGA Universe Cup" is taking place somewhere in our universe. Both FIFA World Cup and GIGA Universe Cup are two rounds competitions that consist of the first round, also known as "group league," and the second called "final tournament." In the first round, participating teams are divided into groups of four teams each. Each team in a group plays a match against each of the other teams in the same group. For example, let's say we have a group of the foIowing four teams, "Engband, Swedon, Argontina, and Nigerua." They play the following six matches: Engband - Swedon, Engband - Argontina, Engband - Nigerua, Swedon - Argontina, Swedon - Nigerua, and Argontina - Nigerua.</p>

<p>The result of a single match is shown by the number of goals scored by each team, like "Engband 1 - 0 Argontina," which says Engband scored one goal whereas Argontina zero. Based on the result of a match, points are given to the two teams as follows and used to rank teams. If a team wins a match (i.e., scores more goals than the other), three points are given to it and zero to the other. If a match draws (i.e., the two teams score the same number of goals), one point is given to each.</p>

<p>The <em>goal difference</em> of a team in given matches is the total number of goals it scored minus the total number of goals its opponents scored in these matches. For example, if we have three matches "Swedon 1 - 2 Engband," "Swedon 3 - 4 Nigerua," and "Swedon 5 - 6 Argontina," then the goal diference of Swedon in these three matches is (1 + 3 + 5) - (2 + 4+ 6) = -3.</p>

<p>Given the results of all the six matches in a group, teams are ranked by the following criteria, listed in the order of priority (that is, we first apply (a) to detennine the ranking, with ties broken by (b), with ties broken by (c), and so on).</p>

<p>(a) greater number of points in all the group matches;</p>

<p>(b) greater goal difference in all the group matches;</p>

<p>(c) greater number of goals scored in all the group matches.</p>

<p>If two or more teams are equal on the basis of the above three criteria, their place shoul be determined by the following criteria, applied in this order:</p>

<p>(d) greater number of points obtained in the group matches between the teams concerned;</p>

<p>(e) greater goal difference resulting from the group matches between the teams concerned;</p>

<p>(f) greater number of goals scored in the group matches between the teams concerned;</p>

<p>If two or more teams are stiIl equal, apply (d), (e), and (f) as necessary to each such group. Repeat this until those three rules to equal teams do not make any further resolution. FinaIy, teams that still remain equal are ordered by:</p>

<p>(g) drawing lots by the Organizing Committee for the GIGA Universe Cup.</p>

<p>The two teams coming first and second in each group qualify for the second round.</p>

<p>Your job is to write a program which, given the results of matches played so far in a group and one team specified in the group, calculates the probability that the specified team will qualify for the second round. You may assume each team has played exactly two matches and has one match to play. In total, four matches have been played and two matches are to be played.</p>

<p>Assume the probability that any team scores (exactly) p goals in any match is:</p>

<p><mjx-container class="MathJax" jax="CHTML" display="true" style="font-size: 109%; position: relative;"> <mjx-math display="true" class="MJX-TEX" aria-hidden="true" style="margin-left: 0px; margin-right: 0px;"><mjx-mfrac><mjx-frac type="d"><mjx-num><mjx-nstrut type="d"></mjx-nstrut><mjx-mrow><mjx-mn class="mjx-n"><mjx-c class="mjx-c38"></mjx-c></mjx-mn><mjx-mo class="mjx-n"><mjx-c class="mjx-c21"></mjx-c></mjx-mo></mjx-mrow></mjx-num><mjx-dbox><mjx-dtable><mjx-line type="d"></mjx-line><mjx-row><mjx-den><mjx-dstrut type="d"></mjx-dstrut><mjx-mrow><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D45D TEX-I"></mjx-c></mjx-mi><mjx-mo class="mjx-n"><mjx-c class="mjx-c21"></mjx-c></mjx-mo><mjx-mo class="mjx-n"><mjx-c class="mjx-c28"></mjx-c></mjx-mo><mjx-mn class="mjx-n"><mjx-c class="mjx-c38"></mjx-c></mjx-mn><mjx-mo class="mjx-n" space="3"><mjx-c class="mjx-c2212"></mjx-c></mjx-mo><mjx-mi class="mjx-i" space="3"><mjx-c class="mjx-c1D45D TEX-I"></mjx-c></mjx-mi><mjx-mo class="mjx-n"><mjx-c class="mjx-c29"></mjx-c></mjx-mo><mjx-mo class="mjx-n"><mjx-c class="mjx-c21"></mjx-c></mjx-mo></mjx-mrow></mjx-den></mjx-row></mjx-dtable></mjx-dbox></mjx-frac></mjx-mfrac><mjx-msup><mjx-mrow><mjx-mo class="mjx-s3"><mjx-c class="mjx-c28 TEX-S3"></mjx-c></mjx-mo><mjx-mfrac><mjx-frac type="d"><mjx-num><mjx-nstrut type="d"></mjx-nstrut><mjx-mn class="mjx-n"><mjx-c class="mjx-c31"></mjx-c></mjx-mn></mjx-num><mjx-dbox><mjx-dtable><mjx-line type="d"></mjx-line><mjx-row><mjx-den><mjx-dstrut type="d"></mjx-dstrut><mjx-mn class="mjx-n"><mjx-c class="mjx-c34"></mjx-c></mjx-mn></mjx-den></mjx-row></mjx-dtable></mjx-dbox></mjx-frac></mjx-mfrac><mjx-mo class="mjx-s3"><mjx-c class="mjx-c29 TEX-S3"></mjx-c></mjx-mo></mjx-mrow><mjx-script style="vertical-align: 1.177em;"><mjx-mi class="mjx-i" size="s"><mjx-c class="mjx-c1D45D TEX-I"></mjx-c></mjx-mi></mjx-script></mjx-msup><mjx-msup><mjx-mrow><mjx-mo class="mjx-s3"><mjx-c class="mjx-c28 TEX-S3"></mjx-c></mjx-mo><mjx-mfrac><mjx-frac type="d"><mjx-num><mjx-nstrut type="d"></mjx-nstrut><mjx-mn class="mjx-n"><mjx-c class="mjx-c33"></mjx-c></mjx-mn></mjx-num><mjx-dbox><mjx-dtable><mjx-line type="d"></mjx-line><mjx-row><mjx-den><mjx-dstrut type="d"></mjx-dstrut><mjx-mn class="mjx-n"><mjx-c class="mjx-c34"></mjx-c></mjx-mn></mjx-den></mjx-row></mjx-dtable></mjx-dbox></mjx-frac></mjx-mfrac><mjx-mo class="mjx-s3"><mjx-c class="mjx-c29 TEX-S3"></mjx-c></mjx-mo></mjx-mrow><mjx-script style="vertical-align: 1.177em;"><mjx-texatom size="s" texclass="ORD"><mjx-mn class="mjx-n"><mjx-c class="mjx-c38"></mjx-c></mjx-mn><mjx-mo class="mjx-n"><mjx-c class="mjx-c2212"></mjx-c></mjx-mo><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D45D TEX-I"></mjx-c></mjx-mi></mjx-texatom></mjx-script></mjx-msup><mjx-mtext class="mjx-n"><mjx-c class="mjx-c2C"></mjx-c></mjx-mtext></mjx-math><mjx-assistive-mml unselectable="on" display="block"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mfrac><mrow><mn>8</mn><mo>!</mo></mrow><mrow><mi>p</mi><mo>!</mo><mo stretchy="false">(</mo><mn>8</mn><mo>−</mo><mi>p</mi><mo stretchy="false">)</mo><mo>!</mo></mrow></mfrac><msup><mrow data-mjx-texclass="INNER"><mo data-mjx-texclass="OPEN">(</mo><mfrac><mn>1</mn><mn>4</mn></mfrac><mo data-mjx-texclass="CLOSE">)</mo></mrow><mi>p</mi></msup><msup><mrow data-mjx-texclass="INNER"><mo data-mjx-texclass="OPEN">(</mo><mfrac><mn>3</mn><mn>4</mn></mfrac><mo data-mjx-texclass="CLOSE">)</mo></mrow><mrow data-mjx-texclass="ORD"><mn>8</mn><mo>−</mo><mi>p</mi></mrow></msup><mtext>,</mtext></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">\[\frac{8!}{p!(8-p)!}\left(\frac{1}{4}\right)^p\left(\frac{3}{4}\right)^{8-p} \text{,}\]</span> </mjx-container></p>

<p>for <em>p</em> ≤ 8, and zero for <em>p</em> > 8 . Assume the lot in the step (g) is fair.</p>

### 입력 

 <p>The first line of the input is an integer, less than 1000, that indicates the number of subsequent records.</p>

<p>The rest of the input is the indicated number of records. A single record has the following format:</p>

<pre><empty>  <_>  <team><sub>1</sub>  <_>  <team><sub>2</sub>  <_>  <team><sub>3</sub>  <_>  <team><sub>4</sub>
<team><sub>1</sub>  <_>  <empty>  <_>  <m><sub>12</sub>    <_>  <m><sub>13</sub>    <_>  <m><sub>14</sub>
<team><sub>2</sub>  <_>  <empty>  <_>  <empty>  <_>  <m><sub>23</sub>    <_>  <m><sub>24</sub>
<team><sub>3</sub>  <_>  <empty>  <_>  <empty>  <_>  <empty>  <_>  <m><sub>34</sub>
<team><sub>4</sub>  <_>  <empty>  <_>  <empty>  <_>  <empty>  <_>  <empty></pre>

<p>In the above, <_> is a single underscore (_) and <<em>empty</em>> a sequence of exactly four underscores (____). Each of <<em>team</em>><sub>1</sub>, ... , <<em>team</em>><sub>4</sub> is either an asterisk character (*) followed by exactly three uppercase letters (e.g., *ENG), or an underscore followed by exactly three uppercase letters (e.g., _SWE). The former indicates that it is the team you are asked to calculate the probabiIty of the second round qualification for. You may assume exactly one of <<em>team</em>><sub>1</sub>, ... , <<em>team</em>><sub>4</sub> is marked with an asterisk. Each <<em>m</em>><sub><em>ij</em></sub>(1 ≤ <em>i</em> < <em>j</em> ≤ 4) is a match result between the <<em>team</em>><sub><em>i</em></sub> and <<em>team</em>><sub><em>j</em></sub>. Each match result is either __-_ (i.e., two underscores, hyphen, and another underscore) or of the form _<em>x</em>-<em>y</em> where each of <em>x</em> and <em>y</em> is a single digit (≤ 8) . The former indicates that the corresponding match has not been played, whereas the latter that the result of the match was <em>x</em> goals by <<em>team</em>><sub><em>i</em></sub> and <em>y</em> goals by <<em>team</em>><sub><em>j</em></sub>. Since each team has played exactly two matches, exactly two match results are in the former format.</p>

### 출력 

 <p>The output should consist of <em>n</em> lines where <em>n</em> is the number of records in the input. The <em>i</em>th line should show the probability that the designated team (marked with an asterisk) will qualify for the second round in the <em>i</em>th record.</p>

<p>Numbers should be printed with exactly seven digits after the decimal point. Each number should not contain an error geater than 10<sup>-7</sup>.</p>

