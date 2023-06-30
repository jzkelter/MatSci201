breed [particles particle]

particles-own [
  fx     ; x-component of force vector
  fy     ; y-component of force vector
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
]

globals [
  diameter
  r-min
  eps
  cutoff-dist
  dt
  kb
]

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set-default-shape turtles "circle"
  set dt .05

  set kb 1 / 100 ;; arbitrary constant to make temperature scale friendlier
  set r-min 1.12245
  set eps 1
  set cutoff-dist 2 * r-min
  set diameter .9

  setup-atoms
  init-velocity
  ask particle 0 [set color green pen-down]
  reset-timer
  reset-ticks
end


to setup-atoms
  create-particles num-atoms [
    set shape "circle"
    set size diameter
    set color blue
  ]

  if initial-config = "hexagonal" [

    let l sqrt(num-atoms) ;the # of atoms in a row
    let x-dist r-min
    let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)
    let ypos (- l * x-dist / 2) ;the y position of the first atom
    let xpos (- l * x-dist / 2) ;the x position of the first atom
    let r-num 0  ;the row number
    ask turtles [  ;set the atoms; positions
      if xpos > (l * x-dist / 2)  [  ;condition to start a new row
        set r-num r-num + 1
        set xpos (- l * x-dist / 2) + (r-num mod 2) * x-dist / 2
        set ypos ypos + y-dist
      ]
      setxy xpos ypos  ;if we are still in the same row
      set xpos xpos + x-dist
    ]
  ]

  if initial-config = "Random" [
    ask particles [
      setxy random-xcor random-ycor
    ]
    remove-overlap ;make sure atoms aren't overlapping
  ]

end

to init-velocity
  let v-avg sqrt( 3 * Kb * temp)
  ask particles [
    let a random-float 1  ; a random amount of the total velocity to go the x direction
    set vx sqrt (a * v-avg ^ 2) * positive-or-negative
    set vy sqrt( v-avg ^ 2 - vx ^ 2)  * positive-or-negative
  ]
end


to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end


to remove-overlap
  ask particles [
    while [overlapping] [
      setxy random-xcor random-ycor
    ]
  ]
end

to-report overlapping
  report any? other turtles in-radius r-min
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go
  ask particles [
    update-force-and-velocity
  ]
  ask particles [
    move
  ]
  if constant-temp? [scale-velocities]
  tick-advance dt
  update-plots
end


to-report current-temp
  report (1 / (3 * Kb)) * mean [vx ^ 2 + vy ^ 2] of particles
end


to scale-velocities
  let ctemp current-temp
  (ifelse
    ctemp = 0 and temp != 0 [
      ask particles [init-velocity]
    ]
    ctemp != 0 [
      let scale-factor sqrt( temp / ctemp )  ; if "external" temperature is higher particles will speed up and vice versa
      ask particles [
        set vx vx * scale-factor
        set vy vy * scale-factor
      ]
    ]
  )
end

to update-force-and-velocity  ; particle procedure
  let new-fx 0
  let new-fy 0
  ask other particles in-radius cutoff-dist [
    let r distance myself
    let force calc-force r
    face myself
    rt 180
    set new-fx new-fx + (force * dx)
    set new-fy new-fy + (force * dy)]

  set vx velocity-verlet-velocity vx fx new-fx
  set vy velocity-verlet-velocity vy fy new-fy
  set fx new-fx
  set fy new-fy

end

to-report calc-force [r]
  let cutoff-adjument 0.0034259508797285392 ; = force(cutoff-dist) without the adjustment
  report (- eps / (r ^ 7)) * (2 * (r ^ -6) - 1) + cutoff-adjument
end

to move  ; particle procedure
  ;; Uses velocity-verlet algorithm
  set xcor velocity-verlet-pos xcor vx fx
  set ycor velocity-verlet-pos ycor vy fy
end

to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; position and velocity
  report v + (1 / 2) * (new-a + a) * dt
end
@#$#@#$#@
GRAPHICS-WINDOW
0
165
306
472
-1
-1
33.11111111111112
1
18
1
1
1
0
1
1
1
-4
4
-4
4
1
1
1
ticks
30.0

BUTTON
0
65
110
110
NIL
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

BUTTON
115
65
215
110
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
145
20
310
53
num-atoms
num-atoms
0
30
28.0
1
1
NIL
HORIZONTAL

CHOOSER
0
10
138
55
initial-config
initial-config
"hexagonal" "Random"
1

SLIDER
140
120
310
153
temp
temp
0.00
7
5.0
.1
1
NIL
HORIZONTAL

SWITCH
0
120
134
153
constant-temp?
constant-temp?
1
1
-1000

MONITOR
225
65
311
110
temp
current-temp
3
1
11

PLOT
315
320
485
475
Green particle velocity
Vx
Vy
-1.0
1.0
-1.0
1.0
false
false
"" ""
PENS
"default" 1.0 2 -10899396 true "" "plotxy [vx] of particle 0 [vy] of particle 0"
"y-axis" 1.0 0 -7500403 true "plotxy 0 plot-y-min\nplotxy 0 plot-y-max\n\n" "\n"
"x2-axis" 1.0 0 -7500403 true "plotxy plot-x-min 0 \nplotxy plot-x-max 0 " ""

PLOT
315
10
485
160
Vx distribution
NIL
NIL
-1.0
1.0
0.0
24.0
false
false
"" ""
PENS
"default" 0.1 1 -13345367 true "" "histogram [vx] of particles"
"x=0" 1.0 0 -7500403 true "plotxy 0 plot-y-min\nplotxy 0 plot-y-max\n" ""

PLOT
315
165
485
315
Vy distribution
NIL
NIL
-1.0
1.0
0.0
12.0
false
false
"" ""
PENS
"default" 0.1 1 -13345367 true "" "histogram [vy] of particles"
"x=0" 1.0 0 -7500403 true "plotxy 0 plot-y-min\nplotxy 0 plot-y-max\n" ""

MONITOR
490
10
550
47
mean |Vx|
mean [abs vx] of particles
3
1
9

MONITOR
490
175
545
212
mean |Vy|
mean [abs vy] of particles
3
1
9

TEXTBOX
490
60
570
115
The Vx distribution for all atoms each tick
11
0.0
1

TEXTBOX
490
215
570
280
The Vy distribution for all atoms each tick
11
0.0
1

TEXTBOX
490
340
570
461
The green atom's velocity vector plotted as point in the x-y plane each tick
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model shows atoms interacting and moving according to Newton's laws. This is known as a "molecular dynamics model." In this model, one atom is colored green and traces its path to visualize how any given atom moves over time. It is a visualization of Brownian Motion.

## HOW IT WORKS

Atoms in this model interact according to the Lennard-Jones interatomic potential. The protons and electrons in each atom interact with the electrons and protons in other atoms. The opposite charge particles (protons in one atom and electrons in the other) attract and the same-charged particles repel. At close distances, the repulsion is stronger than the attraction and atoms exert a net repulsive force on each other. At long distances the attraction is stronger than the repulsion and the atoms exert a net attractive force on each other. At a specific distance the two forces balance and the net force they exert on each other is zero. 

At each tick in the model, atoms calculate the forces they feel from all the atoms around them and sum these forces into a single net-force vector. Based on this force, they are accelerated according to F=ma (a=F/m). They then move a small distance at their current velocity. 

The green particle traces the path that it follows. 

## HOW TO USE IT

When you press the SETUP button, a number of atoms are created according to the NUM-ATOMS slider. If INITIAL-CONFIG is set to "random" then the atoms will be positioned randomly. If it is set to "HCP" the atoms will start in a hexagonally close-packed structure. The initial velocities are set so that the temperature of the system equals the TEMP slider. Remember, that temperature is actually a measure of the average kinetic energy of the atoms in the system.

If CONSTANT-TEMP? is on, then the velocities of the particles are scaled each tick to keep the average temperate equal to the TEMP slider. You can think of this as the system being in thermal equilibrium with some outside body that is at that temperature. If CONSTANT-TEMP? is off, then the temperature is allowed to change away from the initial temperature You can think of this as the system being thermally insulated.

## THINGS TO NOTICE

- Notice the path that the green particle takes. Does it look random?
- Notice the distribution of velocities in the x and y directions. 


## THINGS TO TRY

- Try running the simulation at different temperatures. How does this affect the velocity distributions?
- Try starting the atoms in random vs HCP configurations. 

## EXTENDING THE MODEL


## RELATED MODELS


## NETLOGO FEATURES


## CREDITS AND REFERENCES
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

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
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
