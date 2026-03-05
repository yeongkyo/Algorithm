# [Platinum III] 돌멩이 게임 - 19569 

[문제 링크](https://www.acmicpc.net/problem/19569) 

### 성능 요약

메모리: 30612 KB, 시간: 1324 ms

### 분류

수학, 다이나믹 프로그래밍, 게임 이론, 역추적

### 제출 일자

2026년 3월 6일 08:17:41

### 문제 설명

<p>당신은 muse와 함께 아래 규칙으로 게임을 해, 승리해야 한다.</p>

<ul>
	<li>처음에 <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D441 TEX-I"></mjx-c></mjx-mi></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mi>N</mi></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$N$</span></mjx-container>개의 돌이 있으며, 게임은 당신부터 시작한다.</li>
	<li>맨 처음에 당신은 무조건 돌을 한 개 가져가야 한다.</li>
	<li>그 다음 차례부터는 돌을 1개 이상 <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-texatom texclass="ORD"><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D465 TEX-I"></mjx-c></mjx-mi><mjx-mo class="mjx-n" space="3"><mjx-c class="mjx-c2B"></mjx-c></mjx-mo><mjx-mn class="mjx-n" space="3"><mjx-c class="mjx-c31"></mjx-c></mjx-mn></mjx-texatom></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mrow data-mjx-texclass="ORD"><mi>x</mi><mo>+</mo><mn>1</mn></mrow></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">${x+1}$</span></mjx-container>개 이하로 가져갈 수 있다. 이때, <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D465 TEX-I"></mjx-c></mjx-mi></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mi>x</mi></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$x$</span></mjx-container>는 전 차례에 상대방이 가져간 돌멩이의 개수이다.</li>
	<li>마지막 돌을 가져가는 사람이 이긴다.</li>
</ul>

<p>그런데, 사악한 muse는 돌의 개수 <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D441 TEX-I"></mjx-c></mjx-mi></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mi>N</mi></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$N$</span></mjx-container>을 자신이 이길 수밖에 없게 설정해 놓기도 한다! 따라서 당신은 돌의 개수를 보고, 이길 수 없다고 판단되면 첫 수를 두기 전에 게임을 끝내야 한다. muse는 이길 수 있는 경우에서 항상 최선의 수를 둔다. 이때, 게임에서 이길 수 있겠는가?</p>

### 입력 

 <p>돌의 개수 <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D441 TEX-I"></mjx-c></mjx-mi></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mi>N</mi></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$N$</span></mjx-container>이 주어진다. (<mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-mn class="mjx-n"><mjx-c class="mjx-c32"></mjx-c></mjx-mn><mjx-mo class="mjx-n" space="4"><mjx-c class="mjx-c2264"></mjx-c></mjx-mo><mjx-mi class="mjx-i" space="4"><mjx-c class="mjx-c1D441 TEX-I"></mjx-c></mjx-mi><mjx-mo class="mjx-n" space="4"><mjx-c class="mjx-c2264"></mjx-c></mjx-mo><mjx-msup space="4"><mjx-mn class="mjx-n"><mjx-c class="mjx-c31"></mjx-c><mjx-c class="mjx-c30"></mjx-c></mjx-mn><mjx-script style="vertical-align: 0.393em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c35"></mjx-c></mjx-mn></mjx-script></mjx-msup></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><mn>2</mn><mo>≤</mo><mi>N</mi><mo>≤</mo><msup><mn>10</mn><mn>5</mn></msup></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$2 \le N \le 10^5$</span></mjx-container>)</p>

### 출력 

 <p>먼저, 돌의 개수를 보고 당신이 이길 수 있는지 판단하여라. 이길 수 없다고 판단될 경우 NO를 출력하고 프로그램을 바로 종료해야 한다. 이길 수 있다고 판단될 경우 YES를 출력하고 게임을 진행한다.</p>

<p>수를 둘 때는 가져갈 돌의 개수를 정수로 출력해야 한다. 이때 출력하고 난 뒤, 줄을 바꾸고 버퍼를 비워야 한다.</p>

<p>당신이 수를 두고 나면 muse 역시 수를 둔다. muse가 가져간 돌의 개수를 입력받아 저장한 뒤, 다시 당신이 수를 두면 된다.</p>

<p>게임이 끝나거나 당신이 잘못된 수를 둘 경우 (예: 가져간 돌의 개수가 음수이거나, 현재 있는 돌의 개수보다 많은 경우) 다음 수에서 프로그램은 즉시 종료되며, 문제를 틀리게 된다. 당신이 이겼을 경우, 프로그램은 즉시 종료되어야 한다. 그렇지 않으면, 시간 초과 등 예상치 못한 채점 결과를 받을 수 있다.</p>

<p>이길 수 없는 게임을 이길 수 있다고 판단하고 게임을 시작할 경우, 즉시 오답 판정을 받는 것이 아닌, 게임이 모두 진행된 뒤에 오답 판정을 받는 것임에 유의하자. muse는 이길 수 있는 상황에서 항상 최선의 수를 둔다.</p>

