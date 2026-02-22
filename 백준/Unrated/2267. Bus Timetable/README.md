# [Unrated] Bus Timetable - 2267 

[문제 링크](https://www.acmicpc.net/problem/2267) 

### 성능 요약

메모리: 6656 KB, 시간: 60 ms

### 분류

분류 없음

### 제출 일자

2026년 2월 23일 08:54:02

### 문제 설명

<p>A typical bus timetable has a column of stops down the left hand side, followed by a series of columns specifying a particular service and labelled with the number of a bus route. Each entry in the table, i.e. each intersection of the row for a given stop and a column for a given service, will be either blank or contain the time that that service is scheduled to leave that stop.</p>

<p>As you can imagine, producing these timetables by hand is very difficult and error prone, so this is where you come in. Write a program to produce a timetable, given details of the various routes and services.</p>

### 입력 

 <p>Input will consist of a number of scenarios. Each scenario will start with the title of the scenario, followed by a number of lines, one for each route in the scenario. The sequence of scenarios will be terminated by a line consisting of a single ‘#’.</p>

<p>The line for a route will start with the time that the first service for that route leaves the depot (hours and minutes) followed by the interval between services (minutes). This will be followed by a series of pairs of integers representing travel times and bus stops. The list of routes will be terminated by a line containing two zeroes (0 0) and will never contain details of more than 10 routes. Note that routes are numbered implicitly, starting from 1, bus stops are also numbered from 1 and are always visited in order (apart from the return to the depot), there will never be more than 99 stops, and that the depot is effectively bus stop 0 (thus the last bus stop number in the list will always be 0). Note that no buses leave before 6:00 am and all services terminate (i.e. are back in the depot) by midnight, thus you should not generate a service that violates these constraints.</p>

### 출력 

 <p>Output will consist of one timetable per scenario. Each timetable will start with the name of the scenario on a line by itself, followed by a heading line detailing the services (sorted in order of departure time from the depot, see example) followed by as many lines as there are bus stops mentioned in the input. Each line starts with the number of the stop (two columns) followed by an entry for each service. Each entry is 6 columns wide and starts with a ‘|’ followed either by 5 spaces (if that service does not stop there) or a departure time in the form hh:mm, using a 24 hour clock (i.e., including leading zeroes). The line is terminated by a ‘|’. Note that times increase monotonically downwards, but not necessarily across. Leave a blank line between scenarios.</p>

