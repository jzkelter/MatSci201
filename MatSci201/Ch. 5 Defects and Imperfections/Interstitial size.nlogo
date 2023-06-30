breed [atoms atom]
undirected-link-breed [atom-links atom-link] ; links between atoms

atoms-own [
  fx     ; x-component of force vector from last time step
  fy     ; y-component of force vector from last time step
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
  mass   ; mass of atom
  sigma  ; distnace at which intermolecular potential between 2 atoms of this typot-E is 0 (if they are different, we average their sigmas)
  pinned? ; False if the atom isn't pinned in place, True if it is (for boundaries)
  pot-E ; Potential energy of the atom
  selected? ; whether the atom is selected or  not to change its size
  base-color  ; display color for the atom when it isn't selected
]

globals [
  eps ; used in LJ force. Well depth; measure of how strongly particles attract each other
  cutoff-dist ; each atom is influenced by its neighbors within this distance (LJ force)
  dt ; time step for the velocity verlet algorithm
  Kb ; Boltzman's constnat
  link-check-dist ; each atom links with neighbors within this distance
  prev-atom-viz-size  ; previous atom viz size
  LJ-force-cutoff-values  ; since sizes can change, we store the cutoff adjustment values for different sizes of atom pairs
  LJ-pot-cutoff-values  ; since sizes can change, we store the cutoff adjustment values for different sizes of atom pairs
  message1 ; this variable holds a turtle for displaying messages.
  message2 ; this variable holds a turtle for displaying messages.
]

;*******************************************************
;**************** Setup Procedures *********************
;*******************************************************

to setup
  clear-all
  set-default-shape turtles "circle"
  set eps 1
  set cutoff-dist 2.5
  set dt .01
  set Kb (1 / 10)
  set link-check-dist 1.5
  setup-atoms-and-links-and-force-lines
  init-velocity
  setup-interstitial

  set LJ-pot-cutoff-values n-values 41 [0]  ; initialize to zeros so zero as cutoff value will be used in the calculation on next line
  set LJ-force-cutoff-values n-values 41 [0] ; initialize to zeros so zero as cutoff value will be used in the calculation on next line


  set LJ-pot-cutoff-values map  [s -> first LJ-potential-and-force cutoff-dist s s] (range 0 2.05 .05)  ; calculate cutoff values to adjust LJ potential various sigma values
  set LJ-force-cutoff-values map  [s -> last LJ-potential-and-force cutoff-dist s s] (range 0 2.05 .05)  ; calculate cutoff values to adjust LJ for various sigma values

  crt 1 [
    setxy 2.5 2.5
    set size 0
    set message1 self
  ]
  crt 1 [
    setxy 3 2.25
    set size 0
    set message2 self
  ]
  reset-ticks
end

to setup-interstitial
  create-atoms 1 [
    setxy 0.3979209337249389  -0.1689467728642201
    set sigma 0.2
    set mass sigma ^ 2
    set base-color red
    set color red
    set pinned? false
    set selected? true
    set-size
  ]
end

to setup-atoms-and-links-and-force-lines
  let r-min  2 ^ (1 / 6)
  let x-dist r-min ; the distance between atoms in the x direction
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)

  setup-atoms x-dist y-dist

  ask atoms [
    let in-radius-atoms (other atoms in-radius cutoff-dist)
    update-links in-radius-atoms
  ]

  ask atom-links [color-links]
end

to setup-atoms [x-dist y-dist]
  let atoms-per-row 5
  let atoms-per-column 5
  create-atoms atoms-per-row * atoms-per-column [
    init-atom
  ]

  let init-xpos (- atoms-per-row * x-dist / 2)  + 0.4  ;the x position of the first atom
  let ypos (- atoms-per-column * y-dist / 2) ;the y position of the first atom
  let xpos init-xpos
  let row-number 0 ; row number, starts at 0 for easy modulo division
  ask atoms [ ; setting up the HCP structure
    if xpos >= (atoms-per-row * x-dist / 2)  [ ; condition for starting new row
      set row-number row-number + 1
      set xpos init-xpos + (row-number mod 2) * x-dist / 2
      set ypos ypos + y-dist
    ]
    setxy xpos ypos
    set xpos xpos + x-dist
  ]

  ; pin the bottom row

    ask atoms with-min [ycor] [
      set pinned? true
      set shape "circle-X"
    ]
end

to init-atom
  set shape "circle"
  set base-color blue
  set color blue
  set mass 1
  set sigma 1
  set pinned? False
  set selected? false
  set-size
end

to init-velocity
  ; initializes velocity for each atom based on the initial system-temp. All atoms start
  ; with the average velocity for the temperature split randomly between the x velocity and the y velocity
  ask atoms [
    let v-avg sqrt (2 * temp * Kb / mass)
    let a random-float 1  ; a random amount of the total velocity to go the x direction
    set vx sqrt (a * v-avg ^ 2) * positive-or-negative  ; need to square v-avg to get back the
    set vy sqrt( v-avg ^ 2 - vx ^ 2)  * positive-or-negative
  ]

end

to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end


;*******************************************************
;**************** Go Procedures *********************
;*******************************************************
to go
  simulate
  interact
end


to simulate
  update-atom-size-viz

  ask atom-links [ die ]

  ; moving happens before velocity and force update in accordance with velocity verlet
  ask atoms with [not pinned?] [move]

  ask atoms [update-force-and-velocity-and-links]

  control-temp

  ask atom-links [color-links]  ; stylizing/coloring links

  tick-advance dt
  update-plots
end


to interact
  drag-atoms-with-mouse

end


to-report current-temp
  report (1 / (2 * Kb)) * mean [mass * (vx ^ 2 + vy ^ 2)] of atoms with [not pinned?]
end


to control-temp
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

;; *****************************************************
;; *********** Molecular Dynamics Procedures ***********
;; *****************************************************

to move  ; atom procedure, uses velocity-verlet algorithm
  let new_x velocity-verlet-pos xcor vx (fx / mass)
  let new_y velocity-verlet-pos ycor vy (fy / mass)

  ifelse new_x < max-pxcor and new_x > min-pxcor [
    set xcor new_x
  ] [
    ; if the atoms would have moved off the screen, we don't move and set velocity to zero
    set vx 0
  ]

  ifelse new_y < max-pycor and new_y > min-pycor [
    set ycor new_y
  ] [
    ; if the atoms would have moved off the screen, we don't move and set velocity to zero
    set vy 0
  ]
end


to update-force-and-velocity-and-links
  let n-fx 0
  let n-fy 0
  let total-potential-energy 0
  let in-radius-atoms other atoms in-radius cutoff-dist
  ask in-radius-atoms [
    ; each atom calculates the force it feels from its
    ; neighboring atoms and sums these forces
    let r distance myself
    let indiv-pot-E-and-force (LJ-potential-and-force r sigma [sigma] of myself)
    let force item 1 indiv-pot-E-and-force
    set total-potential-energy total-potential-energy + item 0 indiv-pot-E-and-force
    face myself
    rt 180
    set n-fx n-fx + (force * dx)
    set n-fy n-fy + (force * dy)
    ]
  set pot-E total-potential-energy / 2  ; divide by 2 to not double count pot-E for each atom


  ; updating velocity and force
  if not pinned? [
    set vx velocity-verlet-velocity vx (fx / mass) (n-fx / mass)
    set vy velocity-verlet-velocity vy (fy / mass) (n-fy / mass)
    set fx n-fx
    set fy n-fy
  ]

  update-atom-color pot-E
  update-links in-radius-atoms
end


;; *****************************************************
;; ****** Lennard-Jones Potential/Force Procecres ******
;; *****************************************************

to-report LJ-potential-and-force [ r sigma1 sigma2] ; for the force, positive = attractive, negative = repulsive
  let sig (sigma1 + sigma2) / 2
  let third-power (sig / r) ^ 3
  let sixth-power third-power ^ 2
  let twelfth-power sixth-power ^ 2
  let force (-48 * eps / r ) * (twelfth-power - (1 / 2) * sixth-power) - LJ-force-cutoff sig
  let potential (4 * eps * (twelfth-power - sixth-power)) - LJ-pot-cutoff sig
  report list potential force
end


to-report calc-PE
  let U 0  ;; U stands for PE

  ask other atoms in-radius cutoff-dist [
    set U U + calc-pair-PE-with myself
  ]
  report U
end

to-report calc-pair-PE-with [other-atom]
  let PE-and-force LJ-potential-and-force (distance other-atom) sigma  [sigma] of other-atom
  report first PE-and-force
end


to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; velocity, acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end


to-report LJ-force-cutoff [sig]
  report item round (sig / .05) LJ-force-cutoff-values
end

to-report LJ-pot-cutoff [sig]
  report item round (sig / .05) LJ-pot-cutoff-values
end


;; *****************************************************
;; *********      Interaction Procedures      **********
;; *****************************************************

to drag-atoms-with-mouse

  if mouse-down? [
    let close-atoms atoms with [distance-to-mouse < 0.5]
    if any? close-atoms [
      ask min-one-of close-atoms [distance-to-mouse] [
        let oldx xcor
        let oldy ycor
        setxy mouse-xcor mouse-ycor
        ifelse calc-PE > 2 [ ; if energy would be too high, don't let the atom go there.
          setxy oldx oldy
        ] [
          set vx 0
          set vy 0
        ]
      ]
    ]
  ]
end

to-report distance-to-mouse
  report distancexy mouse-xcor mouse-ycor
end

to change-atom-size [change]
  ask atoms with [selected?] [
    set sigma max list 0.2 (sigma + change)
    set mass sigma ^ 2  ; mass is proportional to radius squared (because in 2D)
    set-size
  ]
end



;; *****************************************************
;; ********* Atom and Link Display procedures **********
;; *****************************************************

to update-atom-color [atom-PE] ; updating atom color

  ifelse color-atoms-by-potential-energy? [
    set color scale-color color  atom-PE -6 0
  ] [
    set color base-color
  ]

end

to update-links [in-radius-atoms] ; updating links
  if show-diagonal-right-links? [
    set heading 330
    link-with-atoms-in-cone in-radius-atoms
  ]
  if show-diagonal-left-links? [
    set heading 30
    link-with-atoms-in-cone in-radius-atoms
  ]
  if show-horizontal-links? [
    set heading 90
    link-with-atoms-in-cone in-radius-atoms
  ]

end

to link-with-atoms-in-cone [atom-set]
  let in-cone-atoms (atom-set in-cone link-check-dist 60)
    if any? in-cone-atoms [
      create-atom-link-with min-one-of in-cone-atoms [distance myself]
    ]
end


to color-links
  set thickness .25 ; necessary because the links die and reform every tick
  let min-eq-bond-len 1.1 * ([sigma] of end1 + [sigma] of end2) / 2
  let max-eq-bond-len 1.15 * ([sigma] of end1 + [sigma] of end2) / 2
  (ifelse
    link-length < min-eq-bond-len [
      let tmp-len sqrt(min-eq-bond-len - link-length)
      let tmp-color extract-rgb scale-color red tmp-len 1 (-.2)
      set color insert-item 3 tmp-color (125 + (1 + tmp-len) * 30) ]
    link-length > max-eq-bond-len [
      let tmp-len sqrt (link-length - max-eq-bond-len)
      let tmp-color extract-rgb scale-color yellow tmp-len 1 -.2
      set color insert-item 3 tmp-color (125 + (1 + tmp-len) * 30)]
    [ let tmp-color extract-rgb white
      set color insert-item 3 tmp-color 125 ])
end

to set-shape
  ifelse selected? [
    ;set shape "circle 2"
    set shape "circle-s"
  ] [
    set shape "circle"
  ]
end

to update-atom-size-viz
  if atom-viz-size != prev-atom-viz-size [
    ask atoms [set-size]
  ]
  set prev-atom-viz-size  atom-viz-size
end

to set-size
  set size sigma * atom-viz-size
end

; Copyright 2020 Uri Wilensky.
; See Info tab for full copyright and license.
@#$#@#$#@
GRAPHICS-WINDOW
190
10
538
359
-1
-1
48.6
1
10
1
1
1
0
0
0
1
-3
3
-3
3
1
1
1
ticks
30.0

BUTTON
0
90
80
123
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
85
90
170
123
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
130
175
163
temp
temp
0
.4
0.05
.01
1
NIL
HORIZONTAL

SWITCH
560
150
832
183
color-atoms-by-potential-energy?
color-atoms-by-potential-energy?
0
1
-1000

SWITCH
560
45
790
78
show-diagonal-right-links?
show-diagonal-right-links?
1
1
-1000

SWITCH
560
80
790
113
show-diagonal-left-links?
show-diagonal-left-links?
1
1
-1000

SWITCH
560
115
790
148
show-horizontal-links?
show-horizontal-links?
1
1
-1000

TEXTBOX
560
195
710
213
NIL
11
0.0
1

TEXTBOX
560
195
710
235
Color Key\nLinks:
12
0.0
1

TEXTBOX
565
230
740
248
high compression: dark red
11
13.0
1

TEXTBOX
565
245
835
263
low compression: light red (+ grey tone)
11
18.0
1

TEXTBOX
564
259
714
277
equilibrium: grey
11
5.0
1

TEXTBOX
564
272
834
300
low tension: light yellow (+ grey tone)
11
0.0
1

TEXTBOX
565
288
725
306
high tension: dark yellow
11
44.0
1

TEXTBOX
560
305
715
323
Atoms:
12
0.0
1

TEXTBOX
565
320
835
338
low potential energy: dark blue
11
103.0
1

TEXTBOX
565
335
850
363
high potential energy: light blue (-> white)
11
107.0
1

BUTTON
95
190
180
223
increase-size
change-atom-size .1
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
190
85
223
decrease-size
change-atom-size (- .1)
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
560
10
737
43
atom-viz-size
atom-viz-size
0
1.1
0.8
.1
1
sigma
HORIZONTAL

TEXTBOX
5
170
175
196
For changing interstitial atom
12
0.0
1

TEXTBOX
295
330
445
348
Atoms with X don't move
11
9.9
1

@#$#@#$#@
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

circle+
false
0
Circle -7500403 true true 0 0 300
Rectangle -16777216 true false 0 135 300 165
Rectangle -16777216 true false 135 -15 165 300

circle-dot
true
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 88 88 124

circle-s
false
0
Circle -7500403 true true 0 0 300
Line -1 false 210 60 120 60
Line -1 false 90 90 90 120
Line -1 false 120 150 180 150
Line -1 false 210 180 210 210
Line -1 false 90 240 180 240
Line -7500403 true 90 90 120 60
Line -1 false 120 60 90 90
Line -1 false 90 120 120 150
Line -1 false 180 150 210 180
Line -1 false 210 210 180 240

circle-x
false
0
Circle -7500403 true true 0 0 300
Polygon -16777216 true false 240 30 30 240 60 270 270 60
Polygon -16777216 true false 30 60 240 270 270 240 60 30

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
1
@#$#@#$#@
