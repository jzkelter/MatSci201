globals [
  electron-diameter
  charge-flow
  dt
]

breed [electrons electron]
breed [holes hole]
breed [meters meter]

holes-own [
  charge
  vx
  vy
]

electrons-own [
  charge
  vx
  vy
]

patches-own [pcharge]

;**************************************
;**********SETUP PROCEDURES************
;**************************************
to setup
  clear-all
  set dt 1
  set electron-diameter .5

  ; CREATE ELECTRONS/HOLES
  ask patches [
    set pcolor blue
    if-else doping > 0
    [create-hole]
    [create-electron]
  ]

  ; CREATE CONTACTCS
  ask patches with [pxcor = min-pxcor] [set pcolor grey]
  ask patches with [pxcor = max-pxcor] [set pcolor grey]

  set charge-flow 0  ;; used to calculate average current

  reset-ticks
end

to create-hole
  if random-float 1 < doping [
    set pcharge e-charge  ; the patch has a negative charge since a hole is leaving
    sprout-holes 1 [
      hole-init
    ]
  ]
end

to create-electron
  if random-float 1 < abs doping [
    set pcharge (- e-charge)  ; the patch has a positive charge since an electron is leaving
    sprout-electrons 1 [
      electron-init
    ]
  ]
end


to electron-init
  set color orange
  set shape "circle"
  set size electron-diameter
  set charge e-charge
  init-velocity
end

to hole-init
  set color black
  set shape "circle"
  set size electron-diameter
  set charge (- e-charge)
  init-velocity
end

to init-velocity
  let kb-over-m (1 / 10)
  let v-avg kb-over-m * sqrt temperature
  let x-y-split random-float 1
  set vx v-avg * x-y-split * positive-or-negative
  set vy v-avg * (1 - x-y-split) * positive-or-negative
end

to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end

;***********************************
;**********GO PROCEDURES************
;***********************************

to go
  generate-electron-hole-pairs

  ask charge-carriers [
    calc-force-and-velocity
  ]

  ask charge-carriers [
    move
    annihilate
  ]

;  hold-contacts-charge-constant

  tick-advance dt
  update-plots
end



;**************************************
;********** move procedures **********
;**************************************
to move

  if random-float 1 < scatter-prob * dt [
    init-velocity
  ]
  ; Move using simple Euler method
  let new-x xcor + vx * dt
  let new-y ycor + vy * dt ; There are only forces acting in x direction, so just moves with constant velocity.
  let old-pxcor pxcor
  setxy new-x new-y
  add-to-charge-flow old-pxcor
end

to add-to-charge-flow [old-pxcor] ;; turtle procedure
  (ifelse
    (crossed-right? old-pxcor) [update-charge-flow charge]
    (crossed-left? old-pxcor) [update-charge-flow (- charge)]
  )
end

; I think these can be simplified by just passing in new-x and seeing if it is > max-pxcor or < min-pxcor
to-report crossed-right? [old-pxcor]
  report (old-pxcor >= (max-pxcor - 3)) and (pxcor <= (min-pxcor + 3))
end

to-report crossed-left? [old-pxcor ]
  report (old-pxcor <= (min-pxcor + 3)) and (pxcor >= (max-pxcor - 3))
end

;to scatter ; Scattering maintains speed, just changes direction
;  let x-y-split random-float 1
;  let total-v abs vx + abs vy
;  set vx total-v * x-y-split * positive-or-negative
;  set vy total-v * (1 - x-y-split) * positive-or-negative
;end

to calc-force-and-velocity
  ; calc velocity using simple Euler method
  let fx (charge * force-from-voltage)
  set vx vx + fx * dt  ;; assuming mass is 1 so force = acceleration
end


to-report force-from-voltage
  report (voltage / world-width) ; voltage is linear and force is -du/dx. No negative because voltage is U_left - U_right, so a positive voltage is negative du/dx
end

to update-charge-flow [amount]
  set charge-flow charge-flow + amount
end


;**********************************************
;   electon-hole generation and annihilation
;**********************************************
to generate-electron-hole-pairs
  ask patches [
    if random-float 1 < exp(- band-gap / temperature) [
      sprout-electrons 1 [
        electron-init
        fd random-float 1
      ]
      sprout-holes 1 [
        hole-init
        fd random-float 1
      ]
    ]
  ]
end


to annihilate
  if random-float 1 < recombine-prob [
    ifelse breed = electrons [
      let holes-here holes in-radius electron-diameter
      if (count holes-here) > 0 [
        ask one-of holes-here [die]
        die
      ]
    ] [
      let electrons-here electrons in-radius electron-diameter
      if (count electrons-here) > 0 [
        ask one-of electrons-here [die]
        die
      ]
    ]
  ]

end

;*******************************************
;************** Constants ******************
to-report e-charge
  report -1
end

;***************************************
;**********agent set reporters
;***************************************
to-report charge-carriers
  report (turtle-set electrons holes)
end



;**********************************
;***********Plotting****************
;**********************************
to plot-x-axis
  foreach (list plot-x-min plot-x-max)
  [x-cor -> plotxy x-cor 0]
end

to-report avg-current
  report charge-flow / (ticks + 1)
end
@#$#@#$#@
GRAPHICS-WINDOW
215
10
652
292
-1
-1
13.0
1
10
1
1
1
0
1
1
1
-16
16
-10
10
1
1
1
ticks
30.0

BUTTON
18
56
93
89
NIL
setup\n
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
104
55
189
88
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
18
10
190
43
doping
doping
-1
1
0.0
.1
1
NIL
HORIZONTAL

SLIDER
18
96
190
129
temperature
temperature
0.1
3
0.1
.1
1
NIL
HORIZONTAL

SLIDER
17
135
189
168
band-gap
band-gap
1
10
10.0
1
1
NIL
HORIZONTAL

SLIDER
16
174
190
207
recombine-prob
recombine-prob
0
1
0.1
.1
1
NIL
HORIZONTAL

SLIDER
15
212
191
245
scatter-prob
scatter-prob
0
1
0.3
.1
1
(per tick)
HORIZONTAL

SLIDER
14
252
186
285
voltage
voltage
-2
2
2.0
.1
1
NIL
HORIZONTAL

PLOT
215
298
472
441
Avg. Current
NIL
NIL
0.0
10.0
-3.0
3.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot avg-current"
"pen-1" 1.0 0 -7500403 true "" "plot-x-axis"

BUTTON
12
297
181
330
watch a charge carrier
clear-drawing\nask one-of charge-carriers [pen-down watch-me]
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
13
335
131
368
stop watching
ask charge-carriers [pen-up]\nclear-drawing\nreset-perspective
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
558
343
668
388
NIL
count electrons
17
1
11

@#$#@#$#@
## WHAT IS IT?

This model shows a simplified microscopic process of electrical conduction inside of a semiconductor connected across a voltage. The model is based on a adaptation of Drude's free electron theory to semiconductors and shows how electrical current emerges from the collective movement of many electrons and holes (missing electrons) in a semiconductor.

The model also shows how the current depends on the number of available charge carriers. The number of "extrinsic" charge carriers is determined by "doping," adding impurities to the semiconductor that have either extra or fewer valence electrons than the main semiconducting element. The number "intrinsic" charge carriers is determined by temperature. The higher the temperature is, the more likely it is for an electron to be excited into the conduction band, resulting in two charge carriers: a free electron and hole left behind in the valence band. 

Furthermore, the model shows how current depends on the scattering probability of charge carriers, the probability of electrons and holes recombining when they encounter one another and, of course, the applied voltage. 

## HOW IT WORKS

### General overview
The model is basically a molecular dynamics simulation of electrons and holes moving freely under the influence of an externally applied voltage. Charge carriers only respond to the external voltage. All charge carrier interactions are ignored except that electrons and holes have a probaility of recombining (annihilating) when they meet. The Avg. Current plot shows the average current over the course of the whole simulation run. You can adjust parameters in the middle of a simulation to see what happens, but to get a good current value, you should run the simulation with a single set of parameters. 

### The specific steps the model follows
The steps for each time step in the model are:

1. **generate electron-hole pairs**: Each patch generates an electron-hole pair with probability based on the band-gap and temperature according to exp(- band-gap / temperature). 
2. **charge carriers calculate force and velocity**: all the charge carriers calculate the force they feel. In this model that is simply due to the applied voltage. Then they calculate their new velocity based.
3. **charge carriers move**: all charge carriers have a probability of scattering. This resets their velocity to random direction with speed based on temperature. Then call charge carriers move based on their current velocity. If a charge-carrier crosses the left or right boundary, it updates the *charge-flow* variable for calculating current. In this model, the charge carrier simply re-enters the system on the other side to retain charge-balance. 
4. **electrons and holes annihilate**: Any electrons and holes on the same patch annihilate each other with probability equal to the *recombine-prob* slider.

## HOW TO USE IT
### Running the Model
Press the **setup** button to initialize the model. Then press **go** to watch it run. 

### What each element of the interface does.

The **doping** slider determines how many extrinsic charge carriers there are. A positive doping means there are excess holes. This corresponds to doping the semiconductor (e.g. Si) with an element with fewer valence electrons (e.g. Ga). A negative doping means there are excess electrons. This corresponds to doping the semiconductor with an element with more valence electrons (e.g. Sb). *Note*: this must be set prior to clicking the setup button.

The **temperature** slider sets the temperature (right now in arbitrary units). Temperature affects the initial speed of charge carriers and their speed after a scattering event. It also determines how likely it is for an intrinsic electron-hole pair to form. 

The **band-gap** slider determines the band-gap between the valence and conduction bands for the semi-conductor. In the model, this determines how easy it is for an electron-hole pair to form.

The **recombine-prob** slider determines the probability that an electron will recombine and annihilate when the run into one another. 

The **scatter-prob** slider determines how likely it is for charge carriers to be scattered. 

The **voltage** slider determines the voltage applied across the semiconductor. Positive is defined left to right (positive charges will be pushed to the right if the voltage is positive). 

The **watch a charge carrier** button highlights a single charge carrier (chosen randomly) and traces its trajectory. The **stop watching** button resets the perspective and erases any drawn trajectories. 

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="test V-I" repetitions="2" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000"/>
    <metric>avg-current</metric>
    <enumeratedValueSet variable="recombine-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="band-gap">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doping">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="2"/>
    </enumeratedValueSet>
    <steppedValueSet variable="voltage" first="-1" step="0.1" last="1"/>
    <enumeratedValueSet variable="scatter-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="temp-dependence" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000"/>
    <metric>avg-current</metric>
    <enumeratedValueSet variable="recombine-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="band-gap">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doping">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="temperature" first="0.1" step="0.1" last="3"/>
    <enumeratedValueSet variable="voltage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scatter-prob">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
