; Future Work
; - Figure out how to scale the units of the external voltage appropriately.
; - Make scattering based on temperature instead of set externally.
; - Use meaningful units
; - Initiate velocity based on Boltzman distribution, instead of always the average for the current temperature
; - change the way hold-contacts-charge-constant works so that current calculations still work if there are an odd number of y patches.
; - change average current to update the way other running averages are so you can see it change faster
; - think about whether trapazoidal rule will get much better integrations
; - try out calc-u just using the potential of the closest meter and see how different it is.
; - improve how charge-carriers find left and right meters by using the sorted-meters list
; - retry giving charges to meters based on cloud algorithm and see how different it is


__includes ["../current-utils.nls"]

globals [
  electron-diameter
  sorted-meters
  e-charge
  extra-meters-per-side
  time-ave-built-in-potential
  time-sum-built-in-potential
  recombination-probability
  current-delay
  old-band-gap
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

meters-own [
  x ; the value used for calculating e-field and potential. It is equivalent to xcor for most  meters, but some are outside the bounds of the patches
  ionic-core-charge
  e-count
  running-ave-e-count
  h-count
  running-ave-h-count
  x-charge
  running-ave-charge
  e-field
  running-ave-field
  e-potential
  running-ave-potential
]

;**************************************
;**********SETUP PROCEDURES************
;**************************************
to setup
  clear-all
  set time-ave-built-in-potential 0
  set electron-diameter .5
  set e-charge -1
  set extra-meters-per-side 1
  set recombination-probability 1
  set current-delay 50


  ; CREATE ELECTRONS/HOLES

  make-holes
  make-electrons


  ; CREATE CONTACTCS
  ask left-contacts [set pcolor grey]
  ask right-contacts [set pcolor grey]

  ; CREATE METERS FOR TRACKING CHARGE DENSITY, E-FIELD AND POTENTIAL
  create-the-meters
  setup-battery-viz


  set charge-flow 0  ;; used to calculate average current

  reset-timer
  reset-ticks
end

to make-holes
  let p-patches patches with [pxcor < (min-pxcor + world-width / 2)]
  ask p-patches [
    set pcolor rgb (160 + 40 * p-doping-frac) (170 + 80 * p-doping-frac) (160 + 50 * p-doping-frac)
  ]
  ask n-of (p-doping-frac * count p-patches) p-patches [
    set pcharge e-charge
    sprout-holes 1 [
      hole-init
    ]
  ]
end

to make-electrons
  let n-patches patches with [pxcor >= (min-pxcor + world-width / 2)]
  ask n-patches [
    set pcolor rgb (160 + 40 * n-doping-frac) (160 + 50 * n-doping-frac) (170 + 80 * n-doping-frac)
  ]
  ask n-of (abs n-doping-frac * count n-patches) n-patches [
    set pcharge (- e-charge)
    sprout-electrons 1 [
      electron-init
    ]
  ]
end

to create-the-meters
  let meter-spacing 0.5
  ask patches with [pycor = 0 ] [
    sprout-meters 1 [
      set ionic-core-charge (sum [pcharge] of patches with [ abs (pxcor - [xcor] of myself) < 0.5]) / world-height
      set x xcor
    ]
  ]
  create-meters 1 [
    set xcor 0.5
    set ionic-core-charge 0
    set x xcor
  ]

  foreach (range 1 (extra-meters-per-side + 1))
  [x-cor ->
    create-meters 1 [
      set xcor min-pxcor - .25
      set ycor (- x-cor)
      set x min-pxcor - x-cor]

    create-meters 1 [
      set xcor max-pxcor + .25
      set ycor x-cor
      set x max-pxcor + x-cor]]


  ask charge-carriers [update-meter-charge-counts]
  ask meters [
    set hidden? true
    set running-ave-e-count e-count
    set running-ave-h-count h-count]

  set sorted-meters sort-on [x] meters
end

to electron-init
  set color orange - 2
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
  let kb-over-m 1 / 2
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
  update-battery-viz

  reset-charge-flow-if-voltage-or-scatter-prob-change
  reset-charge-flow-if-band-gap-changes

  remove-stars
  create-electron-hole-pairs

  ask charge-carriers [
    update-force-and-velocity
  ]

  ask charge-carriers [
    move
    annihilate
  ]

  hold-contacts-charge-constant

  calc-field
  calc-built-in-potential

  set flow-timer flow-timer + 1
  tick
end


to move
  ; Move using forward Euler
  let new-x xcor + vx
  let new-y ycor + vy
  facexy new-x new-y
  fd distancexy new-x new-y
end


to update-force-and-velocity
  ifelse (random-float 1 < scatter-prob )  [
    init-velocity  ; if scattered, velocity is reset based on temperature
  ] [
    let k_e 1
    let left-meter min-one-of (meters with [x < [xcor] of myself]) [[xcor] of myself - x]
    let right-meter min-one-of (meters with [x > [xcor] of myself]) [x - [xcor] of myself]
    let del-x (xcor - [x] of left-meter) / ([x] of right-meter - [x] of left-meter) ;  fraction of the way between left and right meters
    let fx k_e * charge * ((1 - del-x) *  ([e-field] of left-meter) + (del-x * ([e-field] of right-meter)))
    set vx vx + fx
  ]
end



;*********************************************************
;******** Generation and Annihilation Procedures *********
;*********************************************************
;***
to create-electron-hole-pairs
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
  if random-float 1 < recombination-probability [
    ifelse breed = electrons [
      let holes-here holes with [distance myself < electron-diameter]
      if (count holes-here) > 0 [
        ask one-of holes-here [die]
        ;die
        set breed turtles
        set shape "star"
        set color yellow
        set size 1
      ]
    ] [
      let electrons-here electrons with [distance myself < electron-diameter]
      if (count electrons-here) > 0 [
        ask one-of electrons-here [
          set shape "star"
          set color yellow
          wait 0.01
          die]
        die
      ]
    ]
  ]

end

to remove-stars
  if erase-stars? [
    ask turtles with [shape = "star"] [die]
  ]
end


;*********************************************************
;* Charge-neutrality and current calculation Proceedures *
;*********************************************************
to hold-contacts-charge-constant
  let left-charge-balance xcor-charge-sum min-pxcor
  let net-charge 0
  update-charge-flow (net-charge - left-charge-balance)
  (ifelse
    left-charge-balance < net-charge [
      create-holes abs (net-charge - left-charge-balance) [
        hole-init
        set xcor min-pxcor + random-float 1 - 0.5
        set ycor random-ycor]
    ]
    left-charge-balance > net-charge [
      ask n-of (left-charge-balance - net-charge) holes with [pxcor = min-pxcor] [die]
    ]
  )

  let right-charge-balance xcor-charge-sum max-pxcor
  update-charge-flow (right-charge-balance - net-charge)
  (ifelse
    right-charge-balance < net-charge [
      ask n-of (net-charge - right-charge-balance) electrons with [pxcor = max-pxcor] [die]
    ]
    right-charge-balance > net-charge [
      create-electrons (right-charge-balance - net-charge) [
        electron-init
        set xcor max-pxcor + random-float 1 - 0.5
        set ycor random-ycor]
    ]
  )
end


to update-charge-flow [amount]
  set charge-flow charge-flow + amount
end





;***************************************
;**********agent set reporters
;***************************************
to-report charge-carriers
  report (turtle-set electrons holes)
end

to-report left-contacts
  report patches with [pxcor = min-pxcor]
end

to-report right-contacts
  report patches with [pxcor = max-pxcor]
end

to-report bulk-meters
  report meters with [x = xcor]
end

to-report boundary-meters
  report meters with [x != xcor]
end


;**************************************************
;**********E-FIELD AND E-POTENTIAL PROCEDURES******
;**************************************************


to-report sign [num]
  report ifelse-value num > 0 [1] [-1]
end

to-report xcor-charge-mean [x-cor]
  report xcor-charge-sum x-cor / world-height
end

to-report xcor-charge-sum [x-cor]
  report ((sum [charge] of charge-carriers with [abs(xcor - x-cor) <= .5]) + (sum [pcharge] of patches with [pxcor = x-cor]))
end


to calc-field
  calc-charge-carrier-distributions
  calc-charges
  calc-e-fields
  calc-potentials
end


to calc-charge-carrier-distributions
  ask meters [
    set e-count 0
    set h-count 0]

  ask charge-carriers [update-meter-charge-counts]

  ask meters [
    set running-ave-h-count weighted-running-ave running-ave-h-count h-count
    set running-ave-e-count weighted-running-ave running-ave-e-count e-count]
end

to update-meter-charge-counts
  let meter-num pxcor - min-pxcor
  let my-xcor xcor
  ask min-one-of meters [abs (x - my-xcor)][
    ;  ask item meter-num sorted-meters [
    ifelse [breed = electrons] of myself
    [set e-count e-count + 1]
    [set h-count h-count + 1]]
end

to calc-charges
    ask bulk-meters [
      set x-charge ionic-core-charge + ((e-charge * e-count ) + (- e-charge * h-count )) / world-height
      set running-ave-charge (weighted-running-ave running-ave-charge x-charge)
  ]

  ask boundary-meters [
    set x-charge (voltage / 2) * (- sign x)
    set running-ave-charge (weighted-running-ave running-ave-charge x-charge)
  ]

end


to calc-e-fields
  ask meters [
    let my-xcor x
    set e-field sum [e-field-force x-charge my-xcor x] of other meters
    set running-ave-field (weighted-running-ave running-ave-field e-field)
  ]
end

to-report e-field-force [q xcor1 xcor2]
  report q / ((xcor1 - xcor2) ^ 2) * sign (xcor1 - xcor2)
end


to calc-potentials
  let last-val 0
  let last-x (min [x] of meters) - 1
  foreach sorted-meters
  [m ->
    ask m [
      set e-potential last-val - (e-field * (x - last-x))
      set running-ave-potential (weighted-running-ave running-ave-potential e-potential)
      set last-val e-potential
      set last-x x]
  ]
end


to-report weighted-running-ave [old-ave new-val]
  let new-weight .05
  report (1 - new-weight) * old-ave + new-weight * new-val
end


to-report built-in-potential
  report [running-ave-potential] of one-of meters with [x = 2] - [running-ave-potential] of one-of meters with [x = -1]
end

to calc-built-in-potential
  let bip built-in-potential
  set time-ave-built-in-potential weighted-running-ave time-ave-built-in-potential bip
  set time-sum-built-in-potential time-sum-built-in-potential + bip
end

to reset-charge-flow-if-band-gap-changes
  if old-band-gap != band-gap [
    set old-band-gap band-gap
    reset-current-meter
  ]
end

to reset-current-meter
  set charge-flow 0
  set flow-timer 0
end

;*******************************************
;***********Plotting and Viz****************
;*******************************************
to plot-x-axis
  foreach (list plot-x-min plot-x-max)
  [x-cor -> plotxy x-cor 0]
end



to plot-y-axis
  foreach (list plot-y-min plot-y-max)
  [y-cor -> plotxy 0.5 y-cor]
end

to set-proper-x-range
;  set-plot-x-range (min-pxcor - extra-meters-per-side) (max-pxcor + extra-meters-per-side)
  set-plot-x-range min-pxcor max-pxcor
end


to update-battery-viz
  if old-voltage != voltage [
    ifelse voltage > 0 [
      ask patches with [pxcor = max-pxcor] [
        set pcolor red
      ]

      ;; now set up the Battery-negative
      ask patches with [pxcor = min-pxcor] [
        set pcolor black
      ]
      ask cathodes [set xcor min-pxcor]
      ask anodes [set xcor max-pxcor]

    ] [
      ask patches with [pxcor = max-pxcor] [
        set pcolor black
      ]

      ;; now set up the Battery-negative
      ask patches with [pxcor = min-pxcor] [
        set pcolor red
      ]

      ask cathodes [set xcor max-pxcor]
      ask anodes [set xcor min-pxcor]

    ]
  ]

  ask (turtle-set cathodes anodes) [
    ifelse voltage = 0 [
      set size 0
    ] [
      set size 0.25 + 0.4 * abs voltage
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
278
30
966
379
-1
-1
17.0
1
10
1
1
1
0
0
1
1
-19
20
-9
10
1
1
1
ticks
30.0

BUTTON
6
59
61
92
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
63
59
118
92
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
1

PLOT
260
750
975
870
Potential
NIL
NIL
-10.0
10.0
-5.0
5.0
false
false
"set-proper-x-range" "clear-plot\nset-proper-x-range\n"
PENS
"Instanteous Potential" 1.0 0 -4539718 true "" "foreach sorted-meters \n[ m -> ask m [plotxy x e-potential]]\n"
"ave-potential" 1.0 0 -16777216 true "" "foreach sorted-meters \n[ m -> ask m [plotxy x running-ave-potential]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"

PLOT
261
380
972
500
Hole (black) and Electron (orange) Distributions
NIL
NIL
0.0
5.0
0.0
30.0
false
false
"set-proper-x-range\n\n" "clear-plot\nset-proper-x-range\n\n"
PENS
"holes" 1.0 0 -16777216 true "" "foreach sorted-meters\n[ m -> ask m [plotxy (x - 0) running-ave-h-count]]\n\n"
"electrons" 1.0 0 -3844592 true "" "foreach sorted-meters \n[ m -> ask m [plotxy (x - 0) running-ave-e-count]]\n"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"
"holes-instantaneous" 1.0 0 -7500403 true "" "; foreach sorted-meters \n; [ m -> ask m [plotxy x h-count]]"
"electrons-instantaneous" 1.0 0 -408670 true "" ";foreach sorted-meters \n;[ m -> ask m [plotxy x e-count]]"
"pen-5" 1.0 0 -4539718 true "" ";plot-x-bounds"

PLOT
261
503
973
623
Charge Density
NIL
NIL
0.0
0.0
-1.0
1.0
false
false
"" "clear-plot\nset-proper-x-range"
PENS
"ave-charge-density" 1.0 1 -16777216 true "" "foreach sorted-meters\n[ m -> ask m [if xcor = pxcor [plotxy (x - 0.5) running-ave-charge]]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"
"instantaneous-charge" 1.0 0 -3026479 true "" ";foreach sorted-meters\n;[ m -> ask m [plotxy x x-charge]]"

PLOT
260
627
974
747
E-Field
NIL
NIL
0.0
10.0
-7.0
2.0
false
false
"set-plot-x-range min-pxcor max-pxcor" "clear-plot\nset-proper-x-range\n;set-plot-y-range round((min [running-ave-field] of meters) - 1) 3"
PENS
"running-ave" 1.0 0 -16777216 true "" "foreach sorted-meters\n[ m -> ask m [plotxy x running-ave-field]]"
"x-axis" 1.0 0 -4539718 true "" "plot-x-axis"
"y-axis" 1.0 0 -4539718 true "" "plot-y-axis"
"instantaneous" 1.0 0 -4539718 true "" "foreach sorted-meters \n[ m -> ask m [plotxy x e-field]]"

SLIDER
7
10
106
43
p-doping-frac
p-doping-frac
0
1
0.5
.1
1
NIL
HORIZONTAL

SLIDER
112
10
212
43
n-doping-frac
n-doping-frac
0
1
0.5
.1
1
NIL
HORIZONTAL

SLIDER
124
110
236
143
band-gap
band-gap
4
10
10.0
1
1
NIL
HORIZONTAL

SLIDER
5
240
198
273
voltage
voltage
-2
2
0.7
.1
1
NIL
HORIZONTAL

SLIDER
5
205
200
238
scatter-prob
scatter-prob
0
1
0.02
.01
1
(per-tick)
HORIZONTAL

MONITOR
279
752
391
797
built-in-potential
time-ave-built-in-potential
2
1
11

BUTTON
6
276
140
309
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
142
275
234
308
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

TEXTBOX
279
11
331
45
P-side
16
0.0
1

TEXTBOX
910
10
976
30
N-side
16
0.0
1

SWITCH
121
59
237
92
erase-stars?
erase-stars?
0
1
-1000

PLOT
0
404
232
554
Vx electrons(orange) & holes(black)
Vx
count
-4.0
4.0
0.0
80.0
false
false
"" "\n "
PENS
"holes distribution" 0.2 1 -16777216 true "" "histogram [vx] of holes"
"y-axis" 1.0 0 -7500403 true "" "foreach (list plot-y-min plot-y-max)\n[y-cor -> plotxy 0 y-cor]"
"hole mean" 1.0 0 -16777216 true "" "plot-pen-reset\nlet mean-vx mean [vx] of holes\nforeach (list plot-y-min plot-y-max)\n[y-cor -> plotxy mean-vx y-cor]"
"electrons distribution" 0.2 1 -3844592 true "" "histogram [vx] of electrons"
"electron mean" 1.0 0 -3844592 true "" "plot-pen-reset\nlet mean-vx mean [vx] of electrons\nforeach (list plot-y-min plot-y-max)\n[y-cor -> plotxy mean-vx y-cor]"

SLIDER
7
110
120
143
temperature
temperature
.1
3
0.1
.1
1
NIL
HORIZONTAL

MONITOR
6
146
164
191
ave intrinsic pairs / tick
exp(- band-gap / temperature) * count patches
3
1
11

MONITOR
0
350
59
395
NIL
current
2
1
11

BUTTON
65
355
209
388
NIL
reset-current-meter
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

This model simulates a *p-n* junction using simplified microscopic process. Charge carriers are accellerated classically and the main phenomenological characteristics of the  *p-n* junction emerge. Users can asses the development of charge carrier distribution profiles, fields, and current-voltage characteristics.


## HOW IT WORKS

### General overview
The model consists of a *p*-type semiconductor (left side) and an *n*-type semiconductor (right side) that are physically connected. The number of initial free electrons and holes are determined by the *p*-doping and *n*-doping sliders. After the model starts, electron-hole pairs can also be generated through thermal excitation. There are also two contact regions on either side of the junction.

Electrons and holes move based on Newton's laws subject to the surrounding electric field and are subject to probabilistic scattering. When electrons and holes encounter one another they recombine. The electric field is calculated based on the charges of all the electrons and holes and the ionic cores (patches) as well as the applied external voltage.

### Setup
Each patch represents an atom in the semiconductor lattice. Patches on the left are inititalized with a free hole with probability proportional to the *p*-doping percentage set by the slider. Ppatches on the right start with a free electron on them with probability proportional to the *n*-doping percentage set by the slider. Other settings are available below.

### Go (what happens each tick)
#### Electron-hole pair generation
1. With probability equal to exp(-*E*<sub>g</sub>/*T*), each patch generates an electron-hole pair.

#### Charge-carrier behavior
Each tick, the charge carriers:

1. Update the force they experience in the *x*-direction based on the electric field. Since forces in the *y*-direction should average to zero over time, *y*-direction forces are neglected.
2. Update velocity: Scatter with some probability (scattering resets velocity to a random direction with speed based on temperature). Otherwise update velocity based on calcaulted field values.
3. Move based on current velocity.
4. Recombine with charge carrier of opposite charge if located on the same patch.

#### Calculating current
Charge-carriers move, but charge-balance is never violated. In this model, charge balance is maintained at the boundaries and used to calcuate current:

1. Net charge at each contact is calculated by summing the charges of the patches (ionic cores) and the charge carriers.
2. If net charge isn't zero at the contact, then charge carriers are either removed (equivalent to continuing around the circuit to a battery terminal) or added (equivalent to entering this part of the circuit).
3. Charge-flow is updated based on the number of charges added or removed. This is then used to calculate current.


#### Calculating the electric field.
It is possible for each charge carrier to calculate the forces it feels from all other charge carriers in the simulation. However, this is very computationally expensive. Instead, we use a mesh model that approximates the charge density at each point on a mesh (a grid) and then charge carriers feel the force from the electric field calculated from the mesh. In this model the mesh is 1 dimensional in the *x*-direction, because all forces in the *y*-direction should average to zero over time. There are a set of agents called "meters" which keep track of variables for this calculation. There is one meter for each patch in the x direction. Each tick they:

1. Count the number of electrons and holes at their *x*-coordinate. This is displayed on the "Hole and Electron Distributions" graph.
2. Calculate the average charge for their *x*-coordinate based on the number of charge carriers. There are also two "boundary meters" that set their charge based on the external applied voltage. This is displayed on the "Charge Density" graph.
3. Calculate the eletric field at their *x*-cooridnate based on the charge density. This is displayed in the "*E*-field" graph.
4. Calculate the electric potential at their *x*-coordinate by integrating the electric field from left to right up to their *x*-coordinate.

## HOW TO USE IT

### Running the Model
Press the **setup** button to initialize the model. Then press **go** to watch it run.

### What each element of the interface does.
The ***p*-doping%** slider determines how many holes there are on the p-side of the junction. It is the probability that a patch starts with a free hole.

The ***n*-doping%** slider determines how many electrons there are on the n-side of the junction. It is the probability that a patch starts with a free electron.

The **scatter-prob** slider determines how likely it is for charge carriers to be scattered from their current trajectory.

The **band-gap** slider determines the band-gap between the valence and conduction bands for the semi-conductor. In the model, this determines how easy it is for an electron-hole pair to form.

The **temperature** slider sets the temperature (currently in arbitrary units). Temperature affects the initial speed of charge carriers and their speed after a scattering event. It also determines how likely it is for an intrinsic electron-hole pair to form.

The **recombination-probability** slider determines the probability that an electron will recombine and annihilate when the run into one another.

The **external-voltage** slider determines the voltage applied across the semiconductor. Positive is defined left to right (positive charges will be pushed to the right if the voltage is positive). Voltage is created in this model by creating positive or negative charge density at the interfaces which is then converted to an electric field and electric potential as explained the **How it Works** section.


The **current-delay** input determines how many ticks to wait until starting to calculate current. This allows the model to equilibrate before calculating current to get a better calculatin faster.

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

minus
false
14
Rectangle -1 true false 0 90 300 210

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

plus
false
0
Rectangle -1 true false 105 0 195 300
Rectangle -1 true false 0 105 300 195

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
  <experiment name="experiment" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <exitCondition>ticks &gt;= 1000</exitCondition>
    <metric>ticks - current-delay</metric>
    <metric>time-sum-built-in-potential</metric>
    <metric>charge-flow</metric>
    <metric>external-voltage</metric>
    <metric>time-sum-built-in-potential / (ticks / time-step)</metric>
    <metric>external-voltage</metric>
    <metric>charge-flow / (ticks - current-delay)</metric>
    <enumeratedValueSet variable="n-doping%">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-doping%">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="external-voltage" first="-1" step="0.25" last="4"/>
    <enumeratedValueSet variable="current-delay">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="band-gap">
      <value value="25"/>
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
