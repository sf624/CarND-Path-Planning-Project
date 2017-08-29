Trajectory planner mainly consists of 5 follwing sections,
Which is written along line 249 - 553 in main.cpp.

1. Previous path extraction. (ll.254-262)
First 50 points of previous unused trajectory is reused for new trajectory to be sent to simlator.
That means this program only refines the trajectory after 1 seconds.

2. Reference velocity dicision. (ll.265-304)
Here reference velocity is decided by the shortest distance between the car front of ego car.
The reference velocity is used to generate target trajectory in following procedures.
The ego car tries to maintain at least 2 sec intervals between the nearest car in the same lane. (ll.296-298)
If the nearset car is within 1 sec interval, emergency braking is executed. (l.277).

3. Target trajectory candidates generation (ll.308-459)
Three different trajectory which differs in target lane are generated according to reference velocity.
The trajectories are 5 seconds long which was necessary for sufficient collision estimation in next section 4.
(Note that only first 1 seconds are actually executed in simulator, since remaining trajectories are refined,
in next execution of h.onMessage (l.206).)
If the reference velocity is over 25 mph, the trajectory intend to reach target lane in 2 seconds. (ll.408-414)

4. Lane evaluation. (ll.461-510)
The cost of each lanes are evaluated according to numbers of cars and distance from current target lane.
Also collision existence are calculated and stored in different variables.

5. Finite State Machine used for lane change dicision. (ll.513-541)
The FMS has mainly two different variable: target_lane = {0,1,2} & is_in_lane_change = {false, true}.
Which means there exists follwing 6 different states.

(A) Keep lane
(A-1) target_lane = 0 & is_in_lane_change = false
(A-2) target_lane = 1 & is_in_lane_change = false
(A-3) target_lane = 2 & is_in_lane_change = false

(B) Change lane
(B-1) target_lane = 0 & is_in_lane_change = true
(B-2) target_lane = 1 & is_in_lane_change = true
(B-3) target_lane = 2 & is_in_lane_change = true

Note that this FMS does not memorize the lane from which the car departed in lane change.
Although, the lane change works in pratice.

Lane change is executed if the minimum lane cost is smaller than current lane cost
and also there's not collision estimated. (l.527)