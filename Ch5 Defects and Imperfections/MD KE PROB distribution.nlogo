;; the following breed is for the molecular-dynamics-core.nls file
breed [atoms atom]
breed [trackers tracker] ; a breed to track number of atoms with a given KE

atoms-own [
  ;; the following variables are for the molecular-dynamics-core.nls file
  fx     ; x-component of force vector from last time step
  fy     ; y-component of force vector from last time step
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
  mass   ; mass of atom
  sigma  ; distnace at which intermolecular potential between 2 atoms of this typot-E is 0 (if they are different, we average their sigmas)
  atom-PE ; Potential energy of the atom
  pinned? ; False if the atom isn't pinned in place, True if it is (for boundaries)
  base-color  ; display color for the atom when it isn't selected
]

trackers-own [
  KE-bucket
  N
]

globals [
  temp
  bucket-width  ; width of bucket for keeping track of KE distribution
  blue-cutoff
  green-cutoff

  ; plotting
  blue-trackers
  green-trackers
  red-trackers

  ;; MDC globals
  eps
  eps*4
  cutoff-dist
  dt
  kb
  LJ-force-linear-1sig
  LJ-PE-linear-1sig

]



;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set temp initial-temp
  set bucket-width 0.1
  set blue-cutoff 1.5
  set green-cutoff 3
  set-default-shape turtles "circle"
  mdc.setup-constants

  setup-trackers

  if initial-config = "Solid" [mdc.setup-atoms-natoms num-atoms]
  if initial-config = "Random" [mdc.setup-atoms-random num-atoms]

  ask one-of atoms [ mdc.setup-cutoff-linear-functions-1sig ]
  mdc.init-velocity


  set blue-trackers sort trackers with [KE-bucket <= blue-cutoff]
  set green-trackers sort trackers with [KE-bucket > blue-cutoff and KE-bucket <= green-cutoff]
  set red-trackers sort trackers with [KE-bucket > green-cutoff]


  reset-timer
  reset-ticks
end



to setup-trackers
  create-trackers (14 * (1 / bucket-width)) [
   set KE-bucket who * bucket-width
   set N 0
   hide-turtle
  ]
end
;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go

  simulate

  ask atoms [
    track-KE
    color-by-KE
  ]
end


to simulate
  ask links [hide-link]
  mdc.move-atoms
  mdc.update-force-and-velocity-and-PE-1sig


  tick-advance dt
  update-plots
end


to-report KE
  report 0.5 * mass * (vx ^ 2 + vy ^ 2)
end

to track-KE
  let bucket floor (1 + KE * (1 / bucket-width))
  ask tracker bucket [set N N + (1 * dt)] ; only count it as a dt fraction of time in that bucket
end

to color-by-KE
  (ifelse
    KE < blue-cutoff [set color blue]
    KE < green-cutoff [set color green]
    [set color red] ;; else
  )

end


;; Molecular Dynamics Core

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to mdc.setup-constants
  set dt .01
  set kb 0.1 ; just picking a random constant for Kb that makes things work reasonably
  set eps 1
  set eps*4 eps * 4
  set cutoff-dist 2.5 * r-min
end

to-report r-min
  ; this is calculated by taking the derivative of the Lennard-Jones potential with sigma = 1,
  ; setting it equal to zero and solving for r
  report 2 ^ (1 / 6)
end

to mdc.setup-cutoff-linear-functions-1sig
  ; function switches at 0.9 * cutoff-distance
  ; calculates the force and potential at the switch-point, then makes a linear function from that point to the cutoff-distance
  let switch-point cutoff-dist * 0.9
  let third-power (sigma / switch-point) * (sigma / switch-point) * (sigma / switch-point)
  let sixth-power third-power * third-power
  let switch-force (-6 * eps*4 * sixth-power / switch-point) * (2 * sixth-power - 1)
  let switch-potential (eps*4 * (sixth-power * sixth-power - sixth-power))
  let force-slope ((- switch-force) / (cutoff-dist - switch-point))
  let potential-slope ((- switch-potential) / (cutoff-dist - switch-point))

  ; this is basically y = m * x + b, with b being 10 * switch-value
  ; because the y-intercept is 10 times the distance from the cutoff-distance as the switch-point
  set LJ-force-linear-1sig ( [ r -> r * force-slope + switch-force * 10 ] )
  set LJ-PE-linear-1sig  ( [ r -> r * potential-slope + switch-potential * 10 ] )
end


to-report mdc.piecewise-linear-cutoffs-2sig [sig]
  ; function switches at 0.9 * cutoff-distance
  ; calculates the force and potential at the switch-point, then makes a linear function from that point to the cutoff-distance
  let switch-point cutoff-dist * 0.9
  let third-power (sig / switch-point) * (sig / switch-point) * (sig / switch-point)
  let sixth-power third-power * third-power
  let switch-force (-6 * eps*4 * sixth-power / switch-point) * (2 * sixth-power - 1)
  let switch-potential (eps*4 * (sixth-power * sixth-power - sixth-power))
  let force-slope ((- switch-force) / (cutoff-dist - switch-point))
  let potential-slope ((- switch-potential) / (cutoff-dist - switch-point))

  ; this is basically y = m * x + b, with b being 10 * switch-value
  ; because the y-intercept is 10 times the distance from the cutoff-distance as the switch-point
  report list [ r -> r * potential-slope + switch-potential * 10 ] [ r -> r * force-slope + switch-force * 10 ]
end


to mdc.init-atom
  set shape "circle"
  set color blue
  set mass 1
  set sigma 1
  set pinned? false
  set base-color blue
end


to mdc.setup-atoms-nrc [natoms-row natoms-column]
  let xcenter (min-pxcor + max-pxcor) / 2
  let ycenter (min-pycor + max-pycor) / 2

  let x-dist r-min ; the distance between atoms in the x direction
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)

  ; set first atom
  create-atoms 1 [
    mdc.init-atom
    setxy (xcenter - (natoms-row - 1) * x-dist / 2) (ycenter - (natoms-column - 1) * y-dist / 2)
  ]

  ; create first row
  repeat natoms-row - 1 [
    ask atoms with [xcor = max [xcor] of atoms] [
      hatch 1 [
        set heading 90
        forward r-min
      ]
    ]
  ]

  ; create columns
  let pos-neg 1 ; alternate between 1 and -1
  repeat natoms-column - 1 [
    ask atoms with [ycor = max [ycor] of atoms] [
      hatch 1 [
        set heading (pos-neg * 30)
        forward r-min
      ]
    ]
    set pos-neg pos-neg * -1
  ]
end

to mdc.setup-atoms-nrc-rotated [natoms-row natoms-column]
  mdc.setup-atoms-nrc natoms-row natoms-column

  let xcenter (min-pxcor + max-pxcor) / 2
  let ycenter (min-pycor + max-pycor) / 2
  ask atoms [ setxy (xcenter + (ycor - ycenter)) (ycenter - (xcor - xcenter)) ]
end


to mdc.setup-atoms-natoms [natoms]
  ; create a square of atoms >= actual natoms
  let l ceiling sqrt(natoms)
  mdc.setup-atoms-nrc l l

  let diff l ^ 2 - natoms
  ask max-n-of diff atoms [xcor + ycor * 10] [die]
end

to mdc.setup-atoms-natoms-rotated [natoms]
  ; create a square of atoms >= actual natoms
  let l ceiling sqrt(natoms)
  mdc.setup-atoms-nrc-rotated l l

  let diff l ^ 2 - natoms
  ask max-n-of diff atoms [xcor * 10 - ycor] [die]
end


to mdc.setup-atoms-random [natoms]
  create-atoms natoms [
    mdc.init-atom
    setxy random-xcor random-ycor
  ]
  mdc.remove-overlap
end

to mdc.remove-overlap
  ask atoms [
    while [mdc.overlapping] [
      setxy random-xcor random-ycor
    ]
  ]
end

to-report mdc.overlapping
  report any? other atoms in-radius r-min
end


to mdc.pin-bottom-row
  ask atoms with-min [ycor] [
      set pinned? true
      set shape "circle-X"
    ]
end


to mdc.init-velocity
  ask atoms with [not pinned?] [
    let v-avg sqrt (2 * temp * Kb / mass)
    let a random-float 1  ; a random amount of the total velocity to go the x direction
    set vx sqrt (a * v-avg ^ 2) * mdc.positive-or-negative
    set vy sqrt ( v-avg ^ 2 - vx ^ 2)  * mdc.positive-or-negative
  ]
end

to-report mdc.positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to mdc.move-atoms
  ifelse mdc.world-wraps-both? [
    ask atoms with [not pinned?] [mdc.move-wraps-on]
  ]
  [ ask atoms with [not pinned?] [mdc.move-wraps-off] ]
end

to-report mdc.world-wraps-both?
    report [count neighbors4 = 4] of patch max-pxcor max-pycor
end

to-report mdc.world-wraps-x?
  report [count neighbors4 with [pxcor = min-pxcor] = 1] of patch max-pxcor 0
end

to-report mdc.world-wraps-y?
  report [count neighbors4 with [pycor = min-pycor] = 1] of patch 0 max-pycor
end

to mdc.move-wraps-on  ; atom procedure
  ;; Uses velocity-verlet algorithm
  set xcor mdc.velocity-verlet-pos xcor vx (fx / mass)
  set ycor mdc.velocity-verlet-pos ycor vy (fy / mass)
end

to mdc.move-wraps-off
  ;; Uses velocity-verlet algorithm
  let new_x mdc.velocity-verlet-pos xcor vx (fx / mass)
  let new_y mdc.velocity-verlet-pos ycor vy (fy / mass)

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

to mdc.move-atoms-die-at-edge
  ask atoms with [not pinned?] [
    set xcor mdc.velocity-verlet-pos xcor vx (fx / mass)
    set ycor mdc.velocity-verlet-pos ycor vy (fy / mass)

    if xcor > max-pxcor or xcor < min-pxcor or ycor > max-pycor or ycor < min-pycor [
      die ; kills atoms when they move off the world
    ]
  ]
end


to mdc.update-force-and-velocity-and-PE-1sig  ; atom procedure
  ask atoms [
    let force-velocity-PE mdc.calculate-force-and-velocity-and-PE-1sig
    let n-fx item 0 force-velocity-PE
    let n-fy item 1 force-velocity-PE
    let sum-PE item 2 force-velocity-PE

    mdc.update-force-and-velocity n-fx n-fy
  ]
end

to-report mdc.calculate-force-and-velocity-and-PE-1sig
  let n-fx 0
  let n-fy 0
  let sum-PE 0
  ask other atoms in-radius cutoff-dist [
    ; each atom calculates the force it feels from its
    ; neighboring atoms and sums these forces
    let r distance myself
    let indiv-pot-E-and-force (mdc.LJ-potential-and-force-1sig r)
    let force last indiv-pot-E-and-force
    face myself
    rt 180
    set n-fx n-fx + (force * dx)
    set n-fy n-fy + (force * dy)
    set sum-PE sum-PE + first indiv-pot-E-and-force
  ]
  set atom-PE sum-PE / 2  ; divide by 2 to not double count pot-E for each atom

  report (list n-fx n-fy sum-PE)
end



to mdc.update-force-and-velocity [n-fx n-fy]
  ; updating velocity and force
  if not pinned? [
    set vx mdc.velocity-verlet-velocity vx (fx / mass) (n-fx / mass)
      set vy mdc.velocity-verlet-velocity vy (fy / mass) (n-fy / mass)
      set fx n-fx
      set fy n-fy
  ]
end


; this heats or cools the system based on the average temperature of the system compared to the set temp
to mdc.scale-velocities
  let ctemp mdc.current-temp
  (ifelse
    ctemp = 0 and temp != 0 [
      ask atoms [mdc.init-velocity]
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

to-report mdc.current-temp
  ; temperature = kinetic energy / boltzmann constant, with KE = (1/2) * m * v^2
  report mean [ (mass / Kb) * (1 / 2) * (vx ^ 2 + vy ^ 2) ] of atoms with [not pinned?]
end


;; *****************************************************
;; *********      Interaction Procedures      **********
;; *****************************************************

to mdc.drag-atoms-with-mouse-1sig
  if mouse-down? [
    let close-atoms atoms with [mdc.distance-to-mouse < 0.5]
    if any? close-atoms [
      ask min-one-of close-atoms [mdc.distance-to-mouse] [
        let oldx xcor
        let oldy ycor
        setxy mouse-xcor mouse-ycor
        ; this is set for mdc. for now, this is affected by the duplicate LJ
        ifelse mdc.calc-PE-1sig > 1 [ ; if energy would be too high, don't let the atom go there.
          setxy oldx oldy
        ] [
          set vx 0
          set vy 0
        ]
      ]
      display
    ]
  ]
end



to-report mdc.distance-to-mouse
  report distancexy mouse-xcor mouse-ycor
end


;; *****************************************************
;; ****** Lennard-Jones Potential/Force Procedures *****
;; *****************************************************

;; this next section is LJ for single sigma

to-report mdc.LJ-potential-and-force-1sig [r]
  ;; piecewise function in order to not have a discontinuous change at the cutoff-distance
  ;; switches to a linear function that goes to 0 at the cutoff-distance, starting at 0.9 * cutoff-distance
  ;; https://www.desmos.com/calculator/w302kgs0ed
  let switch-point cutoff-dist * 0.9
  ifelse r <= switch-point [
    ; The following is the Lennard-Jones potential and its derivative, just rearranged to optimize calculation speed.
    ; The unoptimized calculation would look like this:
    ; let force eps*4 * (-12 * (sig ^ 12 / r ^ 13) + 6 * (sig ^ 6 / r ^ 7))
    ; let potential eps*4 * ((sig / r) ^ 12 - (sig / r) ^ 6)
    let third-power (sigma / r) * (sigma / r) * (sigma / r)
    let sixth-power third-power * third-power
    let force (-6 * eps*4 * sixth-power / r) * (2 * sixth-power - 1)
    let potential (eps*4 * (sixth-power * sixth-power - sixth-power))
    report list potential force
  ]
  [
    ; Switches to a linear function near the cutoff distance to avoid a discontinuity in the function
    let force ( runresult LJ-force-linear-1sig r )
    let potential (runresult LJ-PE-linear-1sig r )

    report list potential force
  ]
end

to-report mdc.calc-PE-1sig
  let U 0  ;; U stands for PE

  ask other atoms in-radius cutoff-dist [
    set U U + mdc.calc-pair-PE-with-1sig myself
  ]
  report U
end

to-report mdc.calc-pair-PE-with-1sig [other-atom]
  let PE-and-force mdc.LJ-potential-and-force-1sig distance other-atom
  report first PE-and-force
end



;; velocity verlet used by both

to-report mdc.velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report mdc.velocity-verlet-velocity [v a new-a]  ; velocity, previous acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end
@#$#@#$#@
GRAPHICS-WINDOW
180
10
518
349
-1
-1
30
1
18
1
1
1
0
1
1
1
-5
5
-5
5
0
0
1
ticks
30

BUTTON
0
140
80
175
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
140
170
175
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
10
170
43
num-atoms
num-atoms
0
30
30
1
1
NIL
HORIZONTAL

CHOOSER
0
50
170
95
initial-config
initial-config
"Solid" "Random"
1

SLIDER
0
100
170
133
initial-temp
initial-temp
0.1
8
8
.1
1
NIL
HORIZONTAL

PLOT
520
10
790
350
KE distribution
KE
Probability
0
7
0
0.2
false
false
"" "if ticks > 0 [\nclear-plot\n]"
PENS
"average" 0.1 1 -13345367 true "" "foreach blue-trackers [t -> \n  ask t [plotxy (who * bucket-width) (N / (num-atoms * (ticks + 1)))]]"
"pen-1" 0.1 1 -10899396 true "" "foreach green-trackers [t -> \n ask t [plotxy (who * bucket-width) (N / (num-atoms * (ticks + 1)))]]"
"pen-2" 0.1 1 -2674135 true "" "foreach red-trackers [t -> \n ask t [plotxy (who * bucket-width) (N / (num-atoms * (ticks + 1)))]]"
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
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0
-0.2 0 0 1
0 1 1 0
0.2 0 0 1
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@

@#$#@#$#@
