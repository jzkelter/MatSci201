breed [l-atoms l-atom]            ;; lattice atoms
breed [d-atoms d-atom]            ;; diffusing atoms
breed [graph-points graph-point]  ;; these are used to draw the graph when sketching a concentration profile

graph-points-own [
  num-d-atoms        ;; the number of diffusing atoms at this graph-point's x-position
  concentration      ;; the concentration of diffusing atoms at this graph-point's x-position.
]

globals [
  graph-min-ycor       ;; the ycor which the bottom of the graph when sketching a concentration profile
  sorted-graph-points  ;; a list to hold the graph-points in sorted order for plotting purposes
]


;*********************************************
;*************SETUP PROCEDURES****************
;*********************************************
to setup
  clear-all

  set-default-shape turtles "circle"
  set graph-min-ycor (max-pycor - 25)

  ; Create lattice atoms
  ask patches with [pxcor < max-pxcor] [
    sprout-l-atoms 1 [
      set color blue
      set size 0.9
      set heading 45
      fd (sqrt 2) / 2
    ]
  ]

  ; Initialize bulk concentration
  ask patches with [pxcor <= (min-pxcor + initial-columns-populated)] [
    ;; Create an atom with probability equal to the intial bulk concentration
    if random-float 1 < initial-bulk-concentration [
      create-d-atom
    ]
  ]

  ; Initialize surface concentrations
  if constant-left-surface? [set-constant-left-surface]
  if constant-right-surface? [set-constant-right-surface]

  ; Create the graph for drawing concentration profiles
  create-axes
  create-graph-pnts
  hide-graph
  calc-concentration

  reset-ticks
end


; patch procedure to create a diffusing atom
to create-d-atom
    sprout-d-atoms 1 [
      set color yellow
      set size .3
    ]
end

;******************************************
;*************GO PROCEDURES****************
;******************************************

to go
  (ifelse
    go-mode = "simulate" [simulate]
    go-mode = "draw-profile" [draw-profile]
  )
end


; observer procedure to run the simulation. This would normally be
; the go procedure but since there are two go-modes, the code for
; running the simulation is in this procedure
to simulate
  if ticks = 0 [
    hide-graph
    ask l-atoms [show-turtle]
  ]

  if ticks >= max-ticks [stop]

  ; The core of the diffusion model is diffusing atoms executing a random walk
  ask d-atoms [random-walk]

  ; If either of the surfaces are being held at constant concentration, it needs to be set
  if constant-left-surface? [set-constant-left-surface]
  if constant-right-surface? [set-constant-right-surface]

  calc-concentration
  tick
end


; turtle procedure for atoms  to random walk
to random-walk
  ; This random-walk procedure picks a direction and then moves to that site if it is open.
  ; It is important to use this instead of `move-to one-of neighbors4`, so that boundary atoms
  ; have the same probability of moving inwards as non-boundary atoms have to move outwards.
  ; Using `move-to one-of neighbors4` would result in boundary atoms having only three neighbors
  ; and therefore a higher  probability of moving inwards than non-boundary atoms with four neighbrs
  ; have of moving outwards.

  set heading one-of [0 90 180 270]  ; turn to one of 4 neighbors
  let neighbor patch-ahead 1
  if neighbor != nobody and not any? d-atoms-on neighbor [
    move-to neighbor
  ]
end


; procedure to set the left surface concentration based on the left-surface-concentration slider
to set-constant-left-surface
  ask patches with [pxcor = min-pxcor] [
    ask d-atoms-here [die] ;; First remove all atoms on the surface
    if random-float 1 < left-surface-concentration [ ;; Create an atoms on surface with probability equal to surface-concentration
      create-d-atom
    ]
  ]
end


; procedure to set the right surface concentration based on the right-surface-concentration slider
to set-constant-right-surface
  ask patches with [pxcor = max-pxcor] [
    ask d-atoms-here [die]  ;; First remove all atoms on the surface
    if random-float 1 < right-surface-concentration [ ;; Create an atoms on surface with probability equal to surface-concentration
      create-d-atom
    ]
  ]
end


;************************************************************************
;*************CALCULATE AND PLOT CONCENTRATION PROFILE CODE**************
;************************************************************************

; procedure to calculate the concentration of d-atoms at each x-location and set the concentration
; variable of each corresponding graph-point agent
to calc-concentration
  ask graph-points [set num-d-atoms 0]
  ask d-atoms [
    let my-pxcor pxcor
    ask graph-points with [pxcor = my-pxcor] [set num-d-atoms num-d-atoms + 1]
  ]
  ask graph-points [
    set concentration num-d-atoms / world-height
    set ycor ycor-from-concentration
  ]
end

; procedure to plot the concentration based on stored concentration variables
; of the graph points. It is defined here rather than directly in the plot
; so it can be resused multiple times.
to plot-concentration
  foreach sorted-graph-points [gp -> ask gp [plotxy pxcor concentration]]
end


;********************************************************************
;*****************DRAW CONCENTRATIAON PROFILE CODE*******************
;********************************************************************


; hides the graph the user sketches with
to hide-graph
  ask turtles with [any? my-links] [set hidden? true]
  ask links [set hidden? true]
end

; shows the graph the user sketches with
to show-graph
  ask turtles with [any? my-links] [set hidden? false]
  ask links [set hidden? false]
end

; Create the graph-point agents that make up the user-sketched graph
to create-graph-pnts

  ; create a graph-point at each x-position
  ask patches with [pycor = graph-min-ycor] [
    sprout-graph-points 1 [
      set color yellow
      set size 0.5
    ]
  ]

  ; link the graph points together so that they form a single curve
  ask graph-points [
    create-links-with other graph-points with [abs (xcor - [xcor] of myself) = 1] [
      set thickness 0.3
      set color yellow
    ]
  ]

  set sorted-graph-points sort-on [pxcor] graph-points ; for plotting the graph
end


; procedure to allow the user to sketch a concentration profile
to draw-profile
  ask l-atoms [hide-turtle] ; hide the lattice atoms
  show-graph
  setup-plots ; reset the plots since we are resetting the concentration profile by doing this

  if mouse-down? and mouse-ycor >= graph-min-ycor [
    let prev-concentration 0

    ask graph-points with [pxcor = round mouse-xcor] [ ; get the graph point at the user's mouse position
      set ycor mouse-ycor  ; move up to the mouse position
      set prev-concentration concentration
      set concentration concentration-from-ycor

      if prev-concentration != concentration [ ; if the concentration has changed (mouse moved)
        let my-pxcor pxcor
        ; reset the number of d-atoms at this x-position to match the concentration
        ask d-atoms with [pxcor = my-pxcor] [die]
        ask up-to-n-of (world-height * concentration) patches with [pxcor = my-pxcor] [create-d-atom]
      ]
    ]
  ]
  reset-ticks
end


; graph-point procedure to calculate the concentration represented by its ycor
to-report concentration-from-ycor
  report precision ((graph-min-ycor - ycor) / graph-min-ycor) 2
end

; graph-point procedure to calculate what its ycor should be based on its concentration
to-report ycor-from-concentration
  report  graph-min-ycor * (1 - concentration)
end


; create the axes for the graph to sketch concentrations
to create-axes
  let origin-turtle nobody

  ; create a turtle for the origin
  crt 1 [
    setxy (min-pxcor) graph-min-ycor
    set origin-turtle self
    set size 0
  ]

  ; create a turtle for the end of the x-axis and connect it with the origin
  crt 1 [
    setxy max-pxcor graph-min-ycor
    create-link-with origin-turtle [set thickness 0.2 set color white]
    set size 0
  ]

  ; create a turtle for the end of the y-axis and connect it with the origin
  crt 1 [
    setxy min-pxcor max-pycor
    create-link-with origin-turtle [set thickness 0.2 set color white]
    set size 0
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
285
225
789
2666
-1
-1
12.1
1
10
1
1
1
0
0
1
1
0
40
-200
0
1
1
1
ticks
30.0

BUTTON
65
85
155
135
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
0
195
107
254
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

PLOT
267
9
790
222
Concentration Profile 
NIL
NIL
0.0
1.0
0.0
1.0
false
false
"clear-plot\nset-plot-x-range min-pxcor max-pxcor\n\n" ""
PENS
"current" 1.0 0 -7171555 true "" "plot-pen-reset\nplot-concentration\n"
"initial" 1.0 0 -4539718 true "" "if ticks = 0 [plot-concentration]"
"after-15" 1.0 0 -7500403 true "" "if ticks = 15 [plot-concentration]"

SWITCH
0
315
222
348
constant-left-surface?
constant-left-surface?
1
1
-1000

SWITCH
0
393
222
426
constant-right-surface?
constant-right-surface?
1
1
-1000

SLIDER
0
348
222
381
left-surface-concentration
left-surface-concentration
0
1
0.7
.01
1
NIL
HORIZONTAL

SLIDER
0
425
222
458
right-surface-concentration
right-surface-concentration
0
1
0.0
.01
1
NIL
HORIZONTAL

BUTTON
25
260
204
293
show/hide lattice atoms
ask l-atoms [set hidden? not hidden?]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

INPUTBOX
112
194
225
254
max-ticks
1000.0
1
0
Number

SLIDER
0
8
225
41
initial-bulk-concentration
initial-bulk-concentration
0
1
0.36
.01
1
NIL
HORIZONTAL

SLIDER
1
47
226
80
initial-columns-populated
initial-columns-populated
0
world-width - 1
10.0
1
1
NIL
HORIZONTAL

TEXTBOX
115
140
255
180
Select \"draw-profile\" and click and drag in the world to change the concentration profile. 
10
0.0
1

CHOOSER
0
141
105
186
go-mode
go-mode
"simulate" "draw-profile"
0

@#$#@#$#@
## WHAT IS IT?

This is a model of "interstitial" diffusion in solids. The atoms of many solids arrange themselves in an ordered structure called a lattice, and the spaces between the atoms are called "interstitial sites." Smaller atoms that are about the size of the interstitial sites can then diffuse into these sites, which will change the properties of the material. For example, steel is mostly iron with a small fraction of the interstitial sites filled by carbon atoms. 

## HOW IT WORKS
Large blue atoms are created on the corners of patches to form a lattice. The interstitial sites between the lattice atoms are then made up of patches. Small yellow atoms diffuse on the interstitial sites (on patches). Each tick, the diffusing atoms pick a random direction to jump to. If the neighboring interstitial site in that direction is empty, they jump to it, otherwise they stay put. 

The concentration profile graph shows the fraction of interstitial sites that are filled by diffusing atoms. Note, that in real materials, usually on a small fraction of interstitial sites are filled, but in this model a large number, or even all of them can be filled. The grey line shows the concentration at the beginning of the run (or right after you finished sketching a concentration profile) and the yellow line shows the current concentration profile. The x-position on the concentration profile graph lines up with x-position in the view of the atoms. 

## HOW TO USE IT

Press SETUP to initialize the model. The number on the INITIAL-COLUMNS-POPULATED slider determines how many columns of interstitial sites start with atoms on them. The initial concentration of diffusing atoms in these columsn is set by the INITIAL-BULK-CONCENTRATION slider. 

You can use the DRAW-PROFILE button to sketch a concentration profile. The interstitial sites at each x-position will be filled with a certain probability based on how high the concentration profile is that you drew. 

The GO button will run the model. 

HIDE/SHOW LATTICE ATOMS will toggle the visibility of the blue lattice atoms. 

CONSTANT-LEFT-SURFACE? and CONSTANT-RIGHT-SURFACE? determine if the left and right sides of the model are maintained at a constant concentration of diffusing atoms. It is common when materials scientists and engineers process materials to put them in an environment that keeps the concentration of diffusing atoms constant at the surface. For example, immersing a piece of iron in carbon powder would maintain a high surface concentration of diffusing carbon atoms. Conversely, putting the iron in a vacuum would maintain the surface concentration at zero. LEFT-SURFACE-CONCENTRATION and RIGHT-SURFACE-CONCENTRATION determine what each surface's concentration is maintained at, if that surface has a constant concentration. 


## THINGS TO NOTICE

- Can you predict whether the concentration profile will move up or down at each x-position? What does this depend on?
- Notice how quickly or slowly the concentration profile changes at each x-position. What does this depend on?
- The concentration profile will tend to change in one direction, but there are random fluctuations. Why? 


## THINGS TO TRY
- Try turning constant surface concentration off. What happens in the long run?
- Try setting both surfaces to different constant concentrations. What happens? Try out different combinations. 

## EXTENDING THE MODEL
- Try making the probability of jumping to a new site dependent on temperature
- Try removing some of the lattice atoms and making the probability of jumping to a new site dependent on the number of lattice atoms neighboring that site

## NETLOGO FEATURES
To get good statistics on the fraction of interstitial sites filled, this model extends really far in the y-direction. You have to scroll down to see how far the model world actually extends. This is unusual, because normally the whole model world is visible at once. 

This model also uses `mouse-down?`, `mouse-ycor`, and `mouse-xcor` to detect when the user is clicking the mouse and where it is to draw a concentration profile when the GO-MODE is set to "draw-profile". 

The Concentration Profile graph uses if statements with the number of ticks to plot the concentrations at certain points in time and not change after tha. 


## RELATED MODELS

Solid Diffusion
Random Grid Walk Example

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
NetLogo 6.4.0
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
1
@#$#@#$#@
