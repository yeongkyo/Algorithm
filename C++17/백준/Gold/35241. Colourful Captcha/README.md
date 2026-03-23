# [Gold IV] Colourful Captcha - 35241 

[문제 링크](https://www.acmicpc.net/problem/35241) 

### 성능 요약

메모리: 2024 KB, 시간: 0 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 23일 10:36:11

### 문제 설명

<p>It is 2025, and everyone knows about large language models. The world is divided: while some companies want to feed all texts in the world to their next-generation models, some other companies want to keep their texts to themselves and to their human clients.</p>

<p>A new company called <em>CloseAI</em> started to develop a captcha designed to protect human-oriented resources from large language models. To do this, they employ the fact that the latter tend to happily digest and interpret any text, while humans are mostly capable of interpreting the primary message and ignoring the non-important words. The captcha is hence a piece of ASCII art that contains a message for a human and some text in its "pixels" aimed to mislead robots.</p>

<p>The <em>CloseAI</em> representatives hired you to produce the initial prototype, which will be dealing with names of rainbow colours, that is: <code>RED</code>, <code>ORANGE</code>, <code>YELLOW</code>, <code>GREEN</code>, <code>BLUE</code> and <code>VIOLET</code>. Your program will receive two different colour names as input, <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-msub><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D436 TEX-I"></mjx-c></mjx-mi><mjx-script style="vertical-align: -0.15em; margin-left: -0.045em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c31"></mjx-c></mjx-mn></mjx-script></mjx-msub></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>C</mi><mn>1</mn></msub></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$C_1$</span></mjx-container> and <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-msub><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D436 TEX-I"></mjx-c></mjx-mi><mjx-script style="vertical-align: -0.15em; margin-left: -0.045em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c32"></mjx-c></mjx-mn></mjx-script></mjx-msub></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>C</mi><mn>2</mn></msub></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$C_2$</span></mjx-container>. It has to produce a rectangular ASCII art image that would satisfy the following conditions:</p>

<ul>
<li>the number of lines does not exceed 10;</li>
<li>all lines are of the same length;</li>
<li>this length does not exceed 100 characters;</li>
<li>all symbols are either capital English letters or dots;</li>
<li>if interpreting the image as a sequence of strings, there is at least one substring equal to <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-msub><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D436 TEX-I"></mjx-c></mjx-mi><mjx-script style="vertical-align: -0.15em; margin-left: -0.045em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c32"></mjx-c></mjx-mn></mjx-script></mjx-msub></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>C</mi><mn>2</mn></msub></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$C_2$</span></mjx-container> with no extra letters at either end: that is, if $C_2 = <code>RED</code>$, then <code>RED</code> should be present as a substring, but things such as <code>REDHAT</code> or <code>HATRED</code> do not count towards this;</li>
<li>for any other rainbow colour name other than <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-msub><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D436 TEX-I"></mjx-c></mjx-mi><mjx-script style="vertical-align: -0.15em; margin-left: -0.045em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c32"></mjx-c></mjx-mn></mjx-script></mjx-msub></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>C</mi><mn>2</mn></msub></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$C_2$</span></mjx-container>, it must not be present as a substring, even with any extra letters at either end, so that the large language model is directed towards the specific colour at all costs;</li>
<li>the image, interpreted as ASCII art, should contain the word <mjx-container class="MathJax" jax="CHTML" style="font-size: 109%; position: relative;"><mjx-math class="MJX-TEX" aria-hidden="true"><mjx-msub><mjx-mi class="mjx-i"><mjx-c class="mjx-c1D436 TEX-I"></mjx-c></mjx-mi><mjx-script style="vertical-align: -0.15em; margin-left: -0.045em;"><mjx-mn class="mjx-n" size="s"><mjx-c class="mjx-c31"></mjx-c></mjx-mn></mjx-script></mjx-msub></mjx-math><mjx-assistive-mml unselectable="on" display="inline"><math xmlns="http://www.w3.org/1998/Math/MathML"><msub><mi>C</mi><mn>1</mn></msub></math></mjx-assistive-mml><span aria-hidden="true" class="no-mathjax mjx-copytext">$C_1$</span></mjx-container> and nothing else.</li>
</ul>

<p>Obviously, using a human to check the output of your program would be too expensive. Instead of a human, a simplified model is going to be used, which:</p>

<ul>
<li>treats any four-connected\footnote{The neighbours of a cell are the four cells that are: directly above, directly below, directly to the left and directly to the right.} area of English letters as a <em>big letter</em>;</li>
<li>requires that the horizontal bounding boxes of <em>big letters</em> do not intersect or touch;</li>
<li>only counts the number of holes (also four-connected ones) in each <em>big letter</em> to match the colour name with the required one.</li>
</ul>

<p>For the purpose of the latter step, this model assumes that:</p>

<ul>
<li>the letter <code>B</code> has two holes;</li>
<li>the letters <code>A</code>, <code>D</code>, <code>O</code>, <code>P</code>, <code>Q</code>, <code>R</code> have one hole each;</li>
<li>all other letters have zero holes.</li>
</ul>

<p>Obviously, under these conditions all rainbow colour names can be distinguished easily one from another. This model is far from being human-grade though, as it would also read the following image to contain the word "RED" (one hole, zero holes, one hole):</p>

<pre>OOO..P....VV
O.O.PPP..V.V
OOO..P...VVV
</pre>

<p>Now your task is to implement a captcha generator that satisfies all these requirements.</p>

### 입력 

 <ul>
<li>One line containing the colour name which a human should read.</li>
<li>One line containing the colour name which a large language model should read.</li>
</ul>

<p>Both colours are from the provided list of colours and are not equal.</p>

### 출력 

 <p>Output the ASCII art picture satisfying the requirements.</p>

