breed [atoms atom]

atoms-own [
  fx      ; x-component of force vector
  fy      ; y-component of force vector
  vx      ; x-component of velocity vector
  vy      ; y-component of velocity vector
  pinned? ; whether the atom is pinned so it can't move
]

globals [
  r-min
  eps
  cutoff-dist
  dt
  kb
  free-atoms
  distances-list
  distances-rolling-mean
]

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set-default-shape turtles "circle"
  set dt .01
  set kb 1 / 100 ;; arbitrary constant to make temperature scale friendlier
  set r-min 1.12
  set eps 1
  set cutoff-dist 2 * r-min
  set distances-list (list r-min)
  set distances-rolling-mean r-min
  setup-atoms
  init-velocity
  reset-timer
  reset-ticks
end


to setup-atoms
  let num-atoms 23
  create-atoms num-atoms [
    set shape "circle"
    set size 1
    set color blue
    set pinned? false
  ]
  let l sqrt(num-atoms) ;the # of atoms in a row
  let x-dist r-min
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)
  let ypos min-pycor
  let min-x (- l * x-dist / 2) + 0.5
  let xpos min-x ;the x position of the first atom
  let r-num 0  ;the row number
  ask turtles [  ;set the atoms; positions
    if xpos > ((l - 1) * x-dist / 2 + 0.5)  [  ;condition to start a new row
      set r-num r-num + 1
      set xpos min-x + (r-num mod 2) * x-dist / 2
      set ypos ypos + y-dist
    ]
    setxy xpos ypos  ;if we are still in the same row
    set xpos xpos + x-dist
  ]
  ;ask max-n-of 2 (atoms with-min [ycor]) [abs pxcor] [
  ask atoms with-min [ycor] [
    set pinned? true
    set shape "circle 2"
    set color color + 1
  ]
  set free-atoms atoms with [not pinned?]
end


to init-velocity
  let v-avg sqrt( 3 * Kb * temp)
  ask atoms [
    let a random-float 1  ; a random amount of the total velocity to go the x direction
    set vx sqrt (a * v-avg ^ 2) * positive-or-negative
    set vy sqrt( v-avg ^ 2 - vx ^ 2)  * positive-or-negative
  ]
end


to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end




;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go
  increment-temp
  set distances-rolling-mean 0.99 * distances-rolling-mean + 0.01 * mean distances-list
  set distances-list (list) ; reset the list

  ask free-atoms [
    update-force-and-velocity
  ]

  ask free-atoms [
    move
  ]

  scale-velocities

  tick-advance dt
  update-plots
end

to increment-temp
  if auto-increase-temp? and ticks > 10 and (precision ticks 2 mod 10) = 0 and temp < 7 [
    set temp precision (temp + 0.1) 2
    ask atoms [init-velocity]
  ]
end

to-report current-temp
  report (1 / (3 * Kb)) * mean [vx ^ 2 + vy ^ 2] of free-atoms
end


to scale-velocities
  let ctemp current-temp
  (ifelse
    ctemp = 0 and temp != 0 [
      ask atoms [init-velocity]
    ]
    ctemp != 0 [
      let scale-factor sqrt( temp / ctemp )  ; if "external" temperature is higher atoms will speed up and vice versa
      ask atoms [
        set vx vx * scale-factor
        set vy vy * scale-factor
      ]
    ]
  )
end

to update-force-and-velocity  ; atom procedure
  let new-fx 0
  let new-fy 0
  ask other atoms [;in-radius cutoff-dist [
    let r distance myself
    let force calc-force r
    face myself
    set new-fx new-fx + (force * dx)
    set new-fy new-fy + (force * dy)

    if r < (r-min * 1.3) [
      set distances-list lput r distances-list
    ]

  ]

  set vx velocity-verlet-velocity vx fx new-fx
  set vy velocity-verlet-velocity vy fy new-fy
  set fx new-fx
  set fy new-fy

end

to move  ; atom procedure
  ;; Uses velocity-verlet algorithm
  if not pinned? [
    set xcor velocity-verlet-pos xcor vx fx
    set ycor velocity-verlet-pos ycor vy fy
  ]
end

to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; position and velocity
  report v + (1 / 2) * (new-a + a) * dt
end


to-report calc-force [r]
  ;; Remember, that the force function is the derivative of the potential energy function
  let sigma 1
  report -4 * eps * (-12 * (sigma ^ 12 / r ^ 13) + 6 * (sigma ^ 6 / r ^ 7)) ;; derivative of the LJ potential
end
@#$#@#$#@
GRAPHICS-WINDOW
178
10
534
290
-1
-1
38.7143
1
18
1
1
1
0
0
0
1
-4
4
-3
3
1
1
1
ticks
30.0

BUTTON
0
50
86
83
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
90
50
176
83
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
0
137
172
170
temp
temp
.1
7
5.9
.1
1
NIL
HORIZONTAL

PLOT
538
11
821
293
mean neighbor distance vs temp
temperature
mean neighbor distance
0.0
7.0
1.112
1.13
true
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "if ticks > 10 [plotxy temp distances-rolling-mean]"

MONITOR
412
13
533
58
mean neighbor dist
distances-rolling-mean
3
1
11

SWITCH
0
101
174
134
auto-increase-temp?
auto-increase-temp?
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

This is a molecular dynamics (MD) model using the Lennard-Jones potential function. This function models the fact that atoms attract each other when they are a small distance apart and repel each other when they are very close together. By modeling many atoms behaving according to the Lennard-Jones potential, we can see how the bulk behavior of matter at different temperatures emerges from the interactions between discrete atoms. The details of the Lennard-Jones function are discussed in the next section.


## HOW IT WORKS

MD simulations operate according to Newton's laws. The basic steps of the model are as follows. Each tick, each atom:

- Calculates the force that it feels from all other atoms using the Lennard-Jones potential
- Calculates its acceleration based on the net force and its mass using a = F / m
- Updates its velocity based on its acceleration
- Updates its position based on its velocity.

### The Lennard-Jones potential
The Lennard-Jones potential tells you the potential energy of an atom, given its distance from another atom. The derivative of the Lennard-Jones potential tells yout he force an atom feels from another atom based on their distance.

The potential is: V=4ϵ[(σ/r)^12−(σ/r)^6]. Where V is the intermolecular potential between two atoms or molecules, ϵ is depth of the potential well, σ is the distance at which the potential is zero (visualized as the diameter of the atoms), and r is the center-to-center distance of separation between both particles. This is an approximation; the potential function of a real atom depends on its electronic structure and will differ somewhat from the Lennard-Jones potential.
Atoms that are very close will be strongly pushed away from one another, and atoms that are far apart will be gently attracted. Make sure to check out the THINGS TO TRY section to explore the Lennard-Jones potential more in depth.

## HOW TO USE IT

### Simulation initialization
**initial-config**: select if you want the atoms to start out randomly positioned or in a hexagonally-close-packed structure.

**num-atoms**: select the number of atoms to start the simulation

**temp**: select the initial temperature (this will determine the average initial velocity of the atoms.

### Run the simulation

**constant-temp?**: turn this on if you want temperature to be held constant. This can be turned on and off during the simulation run. When it is on, you can move the **temp** slider to change the temperature.

## THINGS TO NOTICE

At low temperatures, the atoms will solidfy and at high temperatures they will break apart (evaporate). Also notice the the packing structure of atoms when they are in a solid and the structures that form when you cool the atoms down after evaporating compared to when they start in HCP.


## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Kelter, J. and Wilensky, U. (2005).  NetLogo Lennard-Jones Molecular Dynamics model.  http://ccl.northwestern.edu/netlogo/models/Electrostatics.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.


## COPYRIGHT AND LICENSE

Copyright 2021 Jacob Kelter and Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.
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
