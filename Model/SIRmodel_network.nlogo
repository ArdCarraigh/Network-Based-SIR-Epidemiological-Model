;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Network-based SIR epidemiological model
;;
;; by Valéria Romano, Cédric Sueur and Maxime Pierron
;;
;; For more information please contact romanodepaula@gmail.com
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [nodes node]
nodes-own [node-id status days-infected] ;; status 0 = susceptible; status 1 = infected; status 2 = recovered

breed  [connections connection]
connections-own [strength]

globals [links-list clock simulations adhesion-time mimetic-coefficient resting-time rank R day]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;; set up procedures ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set simulations 1
  create-file
  reset-ticks
end

to go
  clear-turtles
  clear-all-plots
  import-network
  set clock  0
  set adhesion-time 0
  set day 0
  set rank 1
  ask nodes[
    setxy (random-xcor * 0.95) (random-ycor * 0.95) ;for visual reasons, we don't put any nodes *too* close to the edges
    set status 0
    set days-infected 0
  ]
  output-write "sim"  output-write  simulations output-write "placeholder" output-write "placeholder2" output-write "placeholder3" ; placeholders are just there to satisfy R with 5 elements per line
  output-print " "
  output-write "ind_id" output-write "latency" output-write "total_time" output-write "rank" output-write "day"
  output-print " "
  ifelse mimetism? ; switch on/off: the probability of getting infected takes (or not) into account the number of individuals already infected.
  [dependent]
  [independent]
  ifelse simulations < simulations-number ; this code control the number of total simulations
  [set simulations simulations + 1]
  [
    output-write "sim" output-write  "end" output-write "placeholder" output-write "placeholder2" output-write "placeholder3"
    output-print " "
    export-output "Output.csv" stop go
  ]
  test-simulation
  tick
end


to test-simulation
  file-open "Output.csv"
  file-type simulations
  file-print " "
  file-close
end


to create-file
; Prepare an output file for testing and analysis
  if (file-exists? "Output.csv")[
    carefully [file-delete "Output.csv"]
    [print error-message]
  ]
  ;Writing data to file:
  file-open "Output.csv"
  ; Now write file header
  file-print "Simulation"
  file-print 1
  file-close
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;; transmission procedures ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to dependent
  set resting-time 10 * (max [Caff] of nodes) ; set max R0 * max Caff so that this maximum R0 means a maximum probability of infection of 1 for the most connected individual.
                                              ; All R0 lower than the max and all individuals with lower Caff than the max will result with lower probabilities of infection.
  set mimetic-coefficient ((1 / resting-time)/(number-of-nodes-present)) * R0 ; C = lambda * R0
  ;set mimetic-coefficient 1 ; to set accordingly
  ;set resting-time 1000 ;; It considers the time to be infected without social transmission. Careful setting values here.
                       ;; If 1/resting-time is higher than the mimetic coefficient, then individuals would get infected by chance rather than being infected from others
                       ;; Since we are interested in socially-transmitted pathogens, this value has to be very low.
                       ;; If high values this would mean the mimetic-coefficient is less important than random infections.

  ifelse susceptible! = number-of-nodes-present[
    ask nodes[
      if susceptible! = number-of-nodes-present[
        if random-float 1 <= (size * (((1 / resting-time)/ number-of-nodes-present)*(susceptible!)))[ ;; Probability for the first individual to be infected. STEP 1
                                                                                                    ;; This specifies that all individuals has the same probability
                                                                                                    ;; to be first infected. Once this node is infected, the probability
                                                                                                    ;; of getting infected is no longer the same to all.
                                                                                                    ;; As soon a node is infected, the probability of subsequent nodes getting infected
                                                                                                    ;; will be dependent on the above equation.
          set adhesion-time 0
          set clock 0
          set status 1
          set color pink
          output-write node-id ;; for the first individual
          output-write adhesion-time ; 0 as it is the first to get infected
          output-write clock ; 0 for the same reason
          output-write rank ; 1 as it is the first to get infected
          output-write day
          output-print " "
          set rank rank + 1 ; adds one to the rank for the next infected
        ]
      ]
    ]
  ]
  [

  ask nodes[
    set clock clock + 1 ; clock ticks once per individual. Per day we have number-of-nodes-present ticks.
    set adhesion-time adhesion-time + 1 ; same with latency between infection
    ;if adhesion-time <= time-of-adhesion-limit[ ; if inferior to limit set in the interface, it continues to run, if  superior, then it stops
      if status = 0[
        ifelse aff?[ ;; meaning influence of sociality (social networks) or not
          ifelse any? nodes with [status = 1][
            if random-float 1 <= (size * (((1 / resting-time)/(number-of-nodes-present)) + (Cfoll * infected! * mimetic-coefficient)))[ ; Probability for other nodes getting infected. STEP 2
                                                                                                                               ; Equation page 5 Romano et al. 2016 AJP
              set status 1
              set color red
              output-write node-id
              output-write adhesion-time
              output-write clock
              output-write rank
              output-write day
              output-print " "
              set adhesion-time 0 ; reinitializes the adhseion time between two infections
              set rank rank + 1
            ]
          ]
          [
            if random-float 1 <= (size * (((1 / resting-time)/ number-of-nodes-present)*(susceptible!)))[ ; Probability to get infected when there are no more infected individuals
                                                                                                        ; later on in the outbreak
                                                                                                        ; It won't be used as we stop the simulation when there are no more infected
              set status 1
              set color pink
              output-write node-id
              output-write adhesion-time
              output-write clock
              output-write rank
              output-write day
              output-print " "
              set adhesion-time 0
              set rank rank + 1

            ]
          ]
        ]
        [
          ifelse any? nodes with [status = 1][ ;This is the second step of the simulations turnover - if sociliaty is not entirely taken into account.
            if random-float 1 <= (((1 / resting-time)+(mimetic-coefficient * infected!))/(number-of-nodes-present))[ ;take into account the number of individuals that are informed
                                                                                                                     ; but not the strength of links
              set status 1
              output-write node-id
              set adhesion-time 0
            ]
          ]
          [
            if random-float 1 <= (size * (1 / resting-time)/ number-of-nodes-present)*(susceptible!)[
              set status 1
              output-write node-id
              set adhesion-time 0
            ]
          ]
        ]
      ]
      if status = 1[ ;if the node is infected, chance to recover
        ifelse days-infected < recovery-time[ ; if the individual hasn't spent the required number of days infected for recovery
          set days-infected days-infected + 1 ; we add one to the counter
        ]
        [
          set status 2 ; the individual recover he spent the required number of days infected
          set color grey
        ]
      ]
    ;]
  ]
  ]
  if susceptible! < number-of-nodes-present[ ; change the day only if the first infection has occured
    update-plots
    set day day + 1 ; change day for every run of the [dependant] function
  ]

  ifelse susceptible! < number-of-nodes-present[
    if recovered! < number-of-nodes-present[ ; wait for recovery of all individuals to stop
      if infected! >= 1[ ;or for the outbreak to stop when there are no more infected
        ifelse adhesion-limit?[
          if adhesion-time <= time-of-adhesion-limit
          [dependent]
        ]
        [dependent]
      ]
    ]
  ]
  [dependent]

end

to independent ; no mimetic coefficient + it does not take into account the total number of already infected individuals
ask nodes[
    if status = 0[
      if random-float 1 <= ((1 / resting-time)/ number-of-nodes-present)*(susceptible!)[
        fd 1
        set status 1
        ;do-plot
        output-write adhesion-time
        set adhesion-time 0
      ]
    ]
    if status = 1[
      ifelse patch-here = patch -28 28[
        fd 0
        set status 2
      ]
      [fd 1]
    ]
]
set clock clock + 1
set adhesion-time adhesion-time + 1

ifelse susceptible! < number-of-nodes-present[
    if susceptible! >= 1[
      ifelse adhesion-limit?[
        if adhesion-time <= time-of-adhesion-limit
        [independent]
      ]
      [independent]
    ]
]
[independent]
end



to-report number-of-nodes-present
report count nodes
end

to-report susceptible!
report count nodes with [status = 0]
end

to-report infected!
report count nodes with [status = 1]
end

to-report recovered!
report count nodes with [status = 2]
end

to-report Cfoll
 report ((mean [label] of my-in-links with [[status] of other-end = 1])) ; average of links from wich individuals has a status different than 0
end

to-report Caff
 report mean [label] of my-in-links
end

to-report Csum
  report ((sum [label] of my-out-links with [[status] of other-end = 1]))
end

to import-network
  set-default-shape nodes "circle"
  import-attributes
  import-links
end

;; This procedure reads in a files that contains node-specific attributes
;; including an unique identification number
to import-attributes
  ;; use CAREFULLY to ensure the file is
  ;; closed if there is an error and to notify
  ;; user of the error
  carefully [
    ;; This opens the file, so we can use it.
    file-open attributes_
    ;; Read in all the data in the file
    ;; data on the line is in this order:
    ;; node-id attribute1 attribute2
    while [not file-at-end?]
          [;; this reads a single line into a three-item list
           let items read-from-string (word "[" file-read-line "]")
           create-nodes 1 [set node-id item 0 items
                           set size    item 1 items
                           set color   blue
                          ]
          ]
             ]
  ;; this is the error handling block of the carefully command
  [
    user-message (word "Error reading attributes.txt: " error-message)
  ]
  file-close
end

;; This procedure reads in a file that contains all the links
;; The file is simply 3 columns separated by spaces.  In this
;; example, the links are directed.  The first column contains
;; the node-id of the node originating the link.  The second
;; column the node-id of the node on the other end of the link.
;; The third column is the strength of the link.

to import-links
  carefully
  [
    ;; This opens the file, so we can use it.
    file-open links_
    ;; Read in all the data in the file
    while [not file-at-end?]
    [
      ;; this reads a single line into a three-item list
      let items read-from-string (word "[" file-read-line "]")
      ask get-node (item 0 items)
      [
        create-link-to get-node (item 1 items)
          [ set label item 2 items ;]
            set label-color 1
            set color 32
          ]
      ]
    ]
  ]
  ;; this is the error handling block of the carefully command
  [
    user-message (word "Error reading links.txt: " error-message)
  ]
  file-close
end

;; Helper procedure for looking up a node by node-id.
to-report get-node [id]
  report one-of nodes with [node-id = id]
end
@#$#@#$#@
GRAPHICS-WINDOW
286
18
734
467
-1
-1
6.2
1
10
1
1
1
0
0
0
1
-35
35
-35
35
0
0
1
ticks
30.0

BUTTON
22
23
85
56
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
18
70
214
103
simulations-number
simulations-number
0
10000
10000.0
1
1
NIL
HORIZONTAL

SWITCH
19
209
118
242
mimetism?
mimetism?
0
1
-1000

SWITCH
125
209
215
242
aff?
aff?
0
1
-1000

SWITCH
19
114
142
147
adhesion-limit?
adhesion-limit?
1
1
-1000

SLIDER
19
161
214
194
time-of-adhesion-limit
time-of-adhesion-limit
0
100
100.0
1
1
NIL
HORIZONTAL

BUTTON
105
24
182
57
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

OUTPUT
742
18
1085
478
11

MONITOR
1088
24
1167
69
NIL
susceptible!
17
1
11

MONITOR
1089
150
1167
195
NIL
recovered!
17
1
11

MONITOR
1088
88
1167
133
NIL
infected!
17
1
11

INPUTBOX
20
327
234
387
attributes_
Example_Attributes.txt
1
0
String

INPUTBOX
21
398
236
463
links_
Example_Links.txt
1
0
String

INPUTBOX
19
255
118
315
R0
10.0
1
0
Number

PLOT
1174
13
1498
210
Probability to get infected
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"total node 1" 1.0 0 -7500403 true "" "plot ([(Cfoll * infected! * mimetic-coefficient)] of node 0 + ((1 / resting-time)/(number-of-nodes-present)))"
"node 1 social" 1.0 0 -2674135 true "" "plot [(Cfoll * infected! * mimetic-coefficient)] of node 0"
"node 1 non-social" 1.0 0 -955883 true "" "plot ((1 / resting-time)/(number-of-nodes-present))"

INPUTBOX
130
255
234
315
recovery-time
3.0
1
0
Number

PLOT
1174
225
1498
436
Number of individuals per status
Days
Number of individuals
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"infected" 1.0 0 -16777216 true "" "plotxy day infected!"
"susceptible" 1.0 0 -7500403 true "" "plotxy day susceptible!"
"recovered" 1.0 0 -2674135 true "" "plotxy day recovered!"

@#$#@#$#@
## WHAT IS IT?

Network-Based SIR Epidemiological Model

## HOW IT WORKS

A first random agent gets infected, then the disease is transmitted depending on social connections linking the agents. Infected agents are then being recovered after a number of days set in the interface. A day is passed when the code has looped through each individual to decide to infect it or not during that day.

## HOW TO USE IT

Parametrize the model interface, click on setup and go.

Parameters are defined as follows:

- simulations-number: Set the number of simulations. Int. Range = [0:10000] by default, but should be adapted to the usage.
- adhesion-limit?: Enable/disable the use of the adhesion limit. Bool.
- time-of-adhesion-limit: Set the adhesion limit. When the number of individuals avoiding infection in a row is superior or equal to the adhesion limit, the code stops. Int. Range = [0:100] by default, but should be adapted to the usage.
- mimetism?: Enable/Disable the use of the number of infected individuals in computing probablity of transmission. Bool.
- aff?: Enable/Disable the use of the social network links in computing probability of transmission. Bool.
- R0: Set the basic reproductive rate of the pathogen. Float.
- recovery-time: Set the infectious period of the pathogen. Int.
- attributes_: Name of the attributes table to input to the model. String.
- links_: Name of the links table to input to the model. String.

The attributes table consists of 2 columns WITHOUT header: the individual ID and the susceptibility of individuals. The susceptibility is integrated as the size of the nodes, and should be equal to 1 by default. Individual susceptibility might vary between individuals and should be equal to 0 if the insividual is vaccinated.

The links table is an edgelist of social connections between individuals. It consists of 3 columns WITHOUT header: the ID of the individual emitting the link, the ID of the individual receiving it and the strength of the link.

Examples for both tables are provided with the model files.

On the right side of the interface, you have an output with for each simulation in chronological order:

- id of the infected individual
- time spent since last infection
- time spent since the beginning of the simulation
- the chronological rank of infection
- the day of infection

This output is written in your working directory as "Output.csv" when the code stops running (i.e. when all simulations requested have finished running).

Functions for R Programming Language are provided with the model files in order to ease import and analysis.

On the far right side of the interface you have plots to track the number of individuals in each compartments (Susceptible - Infected - Recovered) as well as the probability of infection of selected individual. The user might edit the plot of probability of infection to track the desired individuals.


## THINGS TO NOTICE

Disabling the visual update allows the model to run faster.
The speed of transmission is relative to the R0 max set in the model. By default it is set to 10 and assumes a probability of transmission of 1 when R0 = 1. The default value moght be adapted to the usage.

## THINGS TO TRY

The user is encouraged to fiddle with R0 and recovery time values in order to observe their respective effects.

## EXTENDING THE MODEL

- Create an infection process that doesn't rely on a arbitrary R0 max.
- Implement a vaccination option, with various strategies (random, trait-based, centrality-based)

## NETLOGO FEATURES

We used netlogo built-in linking functionality to represent the social connetions between individuals.

## RELATED MODELS

Other epidemiological models in the netlogo's model library:

- Biology/Disease Solo
- Biology/HIV
- Biology/Virus 

## CREDITS AND REFERENCES

Model developed by Cédric Sueur, Valéria Romano and Maxime Pierron.

URL: https://github.com/ArdCarraigh/Network-Based-SIR-Epidemiological-Model
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
