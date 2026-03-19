# [level 5] 미로 주행 테스트 - 214295 

[문제 링크](https://school.programmers.co.kr/learn/courses/30/lessons/214295) 

### 성능 요약

메모리: 138 MB, 시간: 2609.67 ms

### 구분

코딩테스트 연습 > 2023 현대모비스 알고리즘 경진대회 본선

### 채점결과

정확성: 100.0<br/>합계: 100.0 / 100.0

### 제출 일자

2026년 03월 20일 08:30:49

### 문제 설명

<p>현대모비스에서는 운전자와 탑승자의 편의를 위한 인포테인먼트 헤드유닛이 탑재되어 있습니다. 인포테인먼트 헤드유닛은 운전자와 탑승자에게 차량내/외의 다양한 정보를 제공하면서 동시에 내비게이션, 음성인식, 텔레매틱스 서비스, 멀티미디어 기능 등 엔터테인먼트 및 편의 서비스를 통합적으로 제공하는 기기입니다. </p>

<p><code>n</code> × <code>m</code> 직사각형 격자 모양의 도로가 나 있는 미로가 있습니다. 미로에서 왼쪽 아래 구석의 좌표는 <code>(0, 0)</code>, 오른쪽 위 구석의 좌표는 <code>(n, m)</code> 입니다. 정수 좌표 <code>(a, b)</code>에는 표지판이 있습니다. 자동차의 내비게이션 기능과 차량의 주행 능력을 검증하기 위해 미로의 다양한 위치에서 출발하여 표지판까지 최단 경로를 따라 이동한 테스트 기록이 있습니다. 각 테스트는 출발점의 좌표, 남은 연료량에 따른 최대 주행 거리, 표지판 도달 여부가 기록되어 있습니다. 출발점의 위치와 표지판의 위치가 같을 수 있으며, 이때는 최대 주행 거리에 상관없이 항상 표지판에 도달한 것으로 간주합니다.</p>

<p>테스트 기록을 토대로 표지판의 위치를 알아내려 합니다. 예를 들어 <code>n</code> = 3, <code>m</code> = 5이고 테스트 기록이 다음과 같은 경우를 생각해 봅시다.</p>
<table class="table">
        <thead><tr>
<th>번호</th>
<th>출발점 좌표</th>
<th>최대 주행 거리</th>
<th>표지판 도달 여부</th>
</tr>
</thead>
        <tbody><tr>
<td>#1</td>
<td>(2, 3)</td>
<td>2</td>
<td>O</td>
</tr>
<tr>
<td>#2</td>
<td>(1, 0)</td>
<td>4</td>
<td>X</td>
</tr>
<tr>
<td>#3</td>
<td>(0, 4)</td>
<td>1</td>
<td>X</td>
</tr>
</tbody>
      </table>
<p>이때 표지판이 있을 수 있는 좌표는 (2, 4), (2, 5), (3, 3), (3, 4)의 4개입니다.<br>
<img src="https://grepp-programmers.s3.ap-northeast-2.amazonaws.com/files/production/65a44b94-6899-4a3a-ac22-31c8c4555a66/pic1.png" title="" alt="pic1.png"></p>

<p>격자의 가로 길이 <code>n</code>, 세로 길이 <code>m</code>, 테스트 기록을 나타내는 2차원 정수 배열 <code>tests</code>가 매개변수로 주어집니다. 표지판이 있을 수 있는 좌표의 개수를 return 하도록 solution 함수를 완성해 주세요.</p>

<hr>

<h5>제한사항</h5>

<ul>
<li>3 ≤ <code>n</code> ≤ 10<sup>9</sup></li>
<li>3 ≤ <code>m</code> ≤ 10<sup>9</sup></li>
<li>1 ≤ <code>tests</code>의 길이 ≤ 250,000

<ul>
<li><code>tests</code>의 원소는 <code>[x, y, d, flag]</code> 형태의 길이가 4인 정수 배열입니다.</li>
<li>출발점의 좌표가 <code>(x, y)</code>, 최대 주행 거리가 <code>d</code>이고, <code>flag</code>가 1인 경우 표지판에 도달했음을, 0인 경우 표지판에 도달하지 못했음을 의미합니다.</li>
<li>0 ≤ <code>x</code> ≤ <code>n</code></li>
<li>0 ≤ <code>y</code> ≤ <code>m</code></li>
<li>0 ≤ <code>d</code> ≤ <code>n + m</code></li>
<li>0 ≤ <code>flag</code> ≤ 1</li>
<li>표지판이 있을 수 있는 좌표가 하나 이상 존재합니다.</li>
</ul></li>
</ul>

<hr>

<h5>입출력 예</h5>
<table class="table">
        <thead><tr>
<th>n</th>
<th>m</th>
<th>tests</th>
<th>result</th>
</tr>
</thead>
        <tbody><tr>
<td>3</td>
<td>5</td>
<td>[[2, 3, 2, 1], [1, 0, 4, 0], [0, 4, 1, 0]]</td>
<td>4</td>
</tr>
<tr>
<td>99999</td>
<td>99999</td>
<td>[[0, 0, 199997, 1]]</td>
<td>9999999999</td>
</tr>
<tr>
<td>99999</td>
<td>99999</td>
<td>[[50000, 50000, 3, 0]]</td>
<td>9999999975</td>
</tr>
<tr>
<td>300</td>
<td>100</td>
<td>[[123, 28, 124, 1], [183, 22, 34, 0], [188, 81, 116, 1], [167, 53, 33, 0], [125, 55, 20, 0]]</td>
<td>6535</td>
</tr>
</tbody>
      </table>
<hr>

<h5>입출력 예 설명</h5>

<p><strong>입출력 예 #1</strong></p>

<p>문제 예시와 같습니다.</p>

<p><strong>입출력 예 #2</strong></p>

<p>총 (99,999 + 1) × (99,999 + 1) = 10,000,000,000개의 정수 좌표 중 (99999, 99999)를 제외한 모든 좌표에 표지판이 있을 수 있습니다.</p>

<p>따라서 9,999,999,999를 return 하면 됩니다.</p>

<p><strong>입출력 예 #3</strong></p>

<p>총 (99,999 + 1) × (99,999 + 1) = 10,000,000,000개의 정수 좌표 중 (50000, 50000)으로부터 거리가 3 이하인 좌표 25개를 제외한 모든 좌표에 표지판이 있을 수 있습니다.</p>

<p>따라서 9,999,999,975를 return 하면 됩니다.</p>

<p><strong>입출력 예 #4</strong></p>

<p>5개의 테스트 기록을 모두 만족하는 정수 좌표의 개수는 총 6,535개임을 알 수 있습니다.</p>

<p>따라서 6,535를 return 하면 됩니다.</p>


> 출처: 프로그래머스 코딩 테스트 연습, https://school.programmers.co.kr/learn/challenges