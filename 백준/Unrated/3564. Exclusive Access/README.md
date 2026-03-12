# [Unrated] Exclusive Access - 3564 

[문제 링크](https://www.acmicpc.net/problem/3564) 

### 성능 요약

메모리: 2172 KB, 시간: 0 ms

### 분류

분류 없음

### 제출 일자

2026년 3월 13일 07:32:31

### 문제 설명

<p>One important problem in concurrent programming is to ensure exclusive access to shared resources by multiple threads. It is also known as Mutual Exclusion protocol. A code that needs to be protected from concurrent execution is called <em>critical section</em> (CS). In order to coordinate access to CS, application threads use a set of shared variables to send information to each other. These shared variables are distinct from all the variables that are used by application code. In practice, mutual exclusion protocol is implemented as two methods — <em>enterCS</em> and <em>exitCS</em>. When application needs to execute some code in CS, it calls enterCS, then executes CS, then calls exitCS.</p>

<p>For theoretical analysis of mutual exclusion protocol one must consider running application as a whole. Each thread of application is represented as an infinite loop that repeatedly performs some work unrelated to CS, which is called <em>non-critical section</em> (NCS), then calls enterCS, then executes CS, then calls exitCS, then the loop repeats. The code inside NCS and CS is not relevant; it is considered to perform no operations related to the protocol and does not modify shared variables used by the protocol.</p>

<p>We consider a system with two concurrently running threads. Threads use a set of shared one-bit variables to implement mutual exclusion protocol. Each variable can store a value of zero or one that can be read or written by a single instruction. Shared variables are initialized to zero. Each thread has a local pointer to the instruction (IP) that it is going to execute next. Execution starts from the top of the code. During each step of execution one of the threads is arbitrarily chosen, it executes one instruction, and then changes its IP to the next instruction to execute. This infinite sequence of steps is called history. A history is called legal if either both threads execute infinitely many steps or just one thread does, while the other thread, having taken a finite number of steps, stops with IP at NCS.</p>

<p>The table below contains several algorithms in pseudo-code that attempt to implement mutual exclusion protocol. In this pseudo-code <em>id</em> is 0 for the first thread and 1 for the second. Variables <em>want</em>[0], <em>want</em>[1], and <em>turn</em> are shared between threads to implement mutual exclusion protocol. Lines marked with “<code>+</code>” implement enterCS, lines marked with “<code>-</code>” implement exitCS. Lines NCS() and CS() are placeholders for some code that works inside non-critical and critical sections respectively and is not relevant for this problem.</p>

<table class="table table-bordered" style="width: 100%;">
	<thead>
		<tr>
			<th>Algorithm 1</th>
			<th>Algorithm 2</th>
			<th>Algorithm 3</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td>
			<pre>loop forever
  NCS()
+ loop while
+   (turn == 1 - id)
CS()
- turn <- (1 - id)
end loop</pre>
			</td>
			<td>
			<pre>loop forever
  NCS()
+ want[id] <- 1
+ loop while
+   (want[1 - id] == 1)
  CS()
- want[id] <- 0
end loop</pre>
			</td>
			<td>
			<pre>loop forever
   NCS()
+  want[id] <- 1
+  turn <- (1 - id)
+  loop while
+    (want[1 - id] == 1 and
+     turn == 1 - id)
   CS()
-  want[id] <- 0
end loop</pre>
			</td>
		</tr>
	</tbody>
</table>

<p>The task is to figure out if the given algorithm satisfies three important properties:</p>

<ul>
	<li>The algorithm satisfies <em>mutual exclusion</em> if in any legal history CS is not executed concurrently by two threads (that is, there is no step where IP of both threads is at CS).</li>
	<li>The algorithm satisfies <em>deadlock freedom</em> if any legal history has infinitely many executions of CS.</li>
	<li>The algorithm satisfies <em>starvation freedom</em> if in any legal history a thread that executes infinitely many steps has infinitely many executions of CS.</li>
</ul>

<p>The property of mutual exclusion is trivial. The algorithm that simply loops forever doing nothing will satisfy it. The sample algorithms above all satisfy mutual exclusion, but the first two fail to achieve deadlock freedom. The algorithm 3 (originally created by Gary Peterson) satisfies all three properties.</p>

### 입력 

 <p>The input file starts with a line with two integer numbers — m<sub>1</sub> and m<sub>2</sub>, where m<sub>i</sub> is the number of lines of code for i-th thread (2 ≤ m<sub>i</sub> ≤ 9). It is followed by m<sub>1</sub> lines with the code for the first thread and m<sub>2</sub> lines with the code for the second thread.</p>

<p>The code for each thread contains one instruction per line. Instruction starts with an integer line number from 1 to m<sub>i</sub> (lines are numbered in ascending order and are included to aid readability), followed by instruction mnemonic, followed by a list of instruction arguments, all separated by spaces. The last arguments of instruction represent line numbers of the next instructions to execute (NIP — from 1 to m<sub>i</sub>). There are three variables shared between threads — A, B, and C. Instruction mnemonics are:</p>

<ul>
	<li>NCS — non-critical section placeholder. Its single argument is NIP.</li>
	<li>CS — critical section placeholder. Its single argument is NIP.</li>
	<li>SET — write value to the shared variable. It has three arguments v, x, and g, where v is the variable to write (A, B, or C), x is the value to write (0 or 1), and g is NIP.</li>
	<li>TEST — read and test the value of the shared variable. It has three arguments v, g<sub>0</sub>, and g<sub>1</sub> where v is the variable to read (A, B, or C), g<sub>0</sub> is NIP if the value of the variable is zero, and g<sub>1</sub> is NIP if the value of the variable is one.</li>
</ul>

<p>NCS and CS appear in the code for each thread exactly once. The code may or may not represent a simple loop, but is guaranteed to alternate executions of CS and NCS by one thread, that is, in every legal history two executions of CS by one thread always have NCS execution by the same thread in between and, vice versa, two executions of NCS by one thread have CS execution by the same thread in between.</p>

### 출력 

 <p>Write to the output file a string of three letters. Letters represent properties of mutual exclusion, deadlock freedom, and starvation freedom. Write letter Y if the corresponding property is satisfied and N otherwise.</p>

