
;; the following breed is for the molecular-dynamics-core.nls file
breed [atoms atom]
undirected-link-breed [atom-links atom-link] ; links between atoms

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

  ;; the following variable is for the atom-editing-procedures.nls file
  selected? ; whether the atom is selected or  not to change its size
]

globals [
  the-interstitial

  ;; globals for mdc. procedures
  eps
  eps*4
  cutoff-dist
  dt
  kb

  LJ-force-linear-1sig
  LJ-PE-linear-1sig

  LJ-force-linear-2sig ; since sizes can change, we store the linear constants for different sizes of atom pairs
  LJ-PE-linear-2sig ; since sizes can change, we store the linear constants for different sizes of atom pairs

  ;; globals for atom-editing-procedures
  prev-atom-viz-size  ; previous atom viz size

  ;; visualizing atoms and bonds globals
  link-check-dist ; each atom links with neighbors within this distance
  min-sigma-for-links
]



;*******************************************************
;**************** Setup Procedures *********************
;*******************************************************

to setup
  clear-all
  mdc.setup-constants
  mdc.setup-cutoff-linear-functions-2sig
  mdc.setup-atoms-nrc 5 5
  mdc.pin-bottom-row
  ask atoms [aep.init-atom]
  setup-interstitial

  mdc.update-force-and-velocity-and-PE-2sig
  mdc.init-velocity

  vab.setup-links

  reset-ticks
end


to setup-interstitial
  create-atoms 1 [
    ; setxy 0.5612310241546858  0.3240268828732776
    setxy 0.021186034506860775 0.6295806514946801
    set shape "circle"
    set color red
    set sigma 0.2
    set mass sigma
    set pinned? false
    set selected? true
    set base-color red
    aep.set-size
    set the-interstitial self
  ]
end


;*******************************************************
;**************** Go Procedures ************************
;*******************************************************

to go
  simulate
  interact
end


to simulate
  aep.update-atom-size-viz

  ask atom-links [ die ]

  ; moving happens before velocity and force update in accordance with velocity verlet
  mdc.move-atoms

  mdc.update-force-and-velocity-and-PE-2sig
  vab.update-atom-color-and-links
  mdc.scale-velocities
  vab.color-links  ; stylizing/coloring links

  tick-advance dt
  update-plots
end

to interact
  mdc.drag-atoms-with-mouse-2sig
end



;; *************************************************************
;; ********** mdc (molecular dynamics core) procedures**********
;; *************************************************************

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

to mdc.setup-cutoff-linear-functions-2sig
  set LJ-PE-linear-2sig n-values 41 [0]  ; initialize to zeros so zero as cutoff value will be used in the calculation on next line
  set LJ-force-linear-2sig n-values 41 [0] ; initialize to zeros so zero as cutoff value will be used in the calculation on next line

  set LJ-PE-linear-2sig map  [s -> first mdc.piecewise-linear-cutoffs-2sig s] (range 0 2.05 .05)  ; calculate cutoff values to adjust LJ potential various sigma values
  set LJ-force-linear-2sig map  [s -> last mdc.piecewise-linear-cutoffs-2sig s] (range 0 2.05 .05)  ; calculate cutoff values to adjust LJ for various sigma values
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

to mdc.update-force-and-velocity-and-PE-2sig
  ask atoms [
    let force-velocity-PE mdc.calculate-force-and-velocity-and-PE-2sig
    let n-fx item 0 force-velocity-PE
    let n-fy item 1 force-velocity-PE
    let sum-PE item 2 force-velocity-PE

    mdc.update-force-and-velocity n-fx n-fy
  ]
end

to-report mdc.calculate-force-and-velocity-and-PE-2sig
  let n-fx 0
  let n-fy 0
  let sum-PE 0
  ask other atoms in-radius cutoff-dist [
    ; each atom calculates the force it feels from its
    ; neighboring atoms and sums these forces
    let r distance myself
    let indiv-pot-E-and-force (mdc.LJ-potential-and-force-2sig r sigma [sigma] of myself)
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

to mdc.drag-atoms-with-mouse-2sig
  if mouse-down? [
    let close-atoms atoms with [mdc.distance-to-mouse < 0.5]
    if any? close-atoms [
      ask min-one-of close-atoms [mdc.distance-to-mouse] [
        let oldx xcor
        let oldy ycor
        setxy mouse-xcor mouse-ycor
        ; this is set for mdc. for now, this is affected by the duplicate LJ
        ifelse mdc.calc-PE-2sig > 1 [ ; if energy would be too high, don't let the atom go there.
          setxy oldx oldy
        ] [
          set vx 0
          set vy 0
        ]
      ]
      ; display  ;; On some models, need this to have it work properly
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



;; this next section is LJ for multiple sigma

to-report mdc.LJ-potential-and-force-2sig [ r sigma1 sigma2] ; for the force, positive = attractive, negative = repulsive
  let sig (sigma1 + sigma2) / 2

  ;; piecewise function in order to not have a discontinuous change at the cutoff-distance
  ;; switches to a linear function that goes to 0 at the cutoff-distance, starting at 0.9 * cutoff-distance
  ;; https://www.desmos.com/calculator/w302kgs0ed
  let switch-point cutoff-dist * 0.9
  ;; https://www.desmos.com/calculator/w302kgs0ed
  ifelse r <= switch-point [
    ; The following is the Lennard-Jones potential and its derivative, just rearranged to optimize calculation speed.
    ; The unoptimized calculation would look like this:
    ; let force eps*4 * (-12 * (sig ^ 12 / r ^ 13) + 6 * (sig ^ 6 / r ^ 7))
    ; let potential eps*4 * ((sig / r) ^ 12 - (sig / r) ^ 6)
    let third-power (sig / r) * (sig / r) * (sig / r)
    let sixth-power third-power * third-power
    let force (-6 * eps*4 * sixth-power / r) * (2 * sixth-power - 1)
    let potential (eps*4 * (sixth-power * sixth-power - sixth-power))
    report list potential force
  ]
  [
    ; Switches to a linear function near the cutoff distance to avoid a discontinuity in the function
    let force ( runresult (mdc.find-LJ-force-linear-2sig sig) r )
    let potential ( runresult (mdc.find-LJ-PE-linear-2sig sig) r )
    report list potential force
  ]
end

to-report mdc.find-LJ-force-linear-2sig [sig]
  ;; finds the appropriate anonymous procedure based on sigma
  report item round (sig / .05) LJ-force-linear-2sig
end

to-report mdc.find-LJ-PE-linear-2sig [sig]
  ;; finds the appropriate anonymous procedure based on sigma
  report item round (sig / .05) LJ-PE-linear-2sig
end

to-report mdc.calc-PE-2sig
  let U 0  ;; U stands for PE

  ask other atoms in-radius cutoff-dist [
    set U U + mdc.calc-pair-PE-with-2sig myself
  ]
  report U
end

to-report mdc.calc-pair-PE-with-2sig [other-atom]
  let PE-and-force mdc.LJ-potential-and-force-2sig (distance other-atom) sigma  [sigma] of other-atom
  report first PE-and-force
end


;; velocity verlet used by both

to-report mdc.velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report mdc.velocity-verlet-velocity [v a new-a]  ; velocity, previous acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end

;; *************************************************************
;; *************** aep (atom editing procedure) ****************
;; *************************************************************

to aep.init-atom
  set selected? false
  aep.set-size
end



to aep.update-atom-size-viz
  if atom-viz-size != prev-atom-viz-size [
    ask atoms [aep.set-size]
  ]
  set prev-atom-viz-size  atom-viz-size
end


to aep.select-atoms
  if mouse-down? [
    ask atoms with [distancexy mouse-xcor mouse-ycor < (1 / 2)] [
      set selected? not selected?
      aep.set-shape
    ]
    wait 0.1
  ]
end

to aep.set-shape
  ifelse selected? [
    ;set shape "circle 2"
    set shape "circle-s"
  ] [
    set shape "circle"
  ]
end


to aep.change-atom-size [change]
  ask atoms with [selected?] [
    set sigma precision (sigma + change) 2
    (ifelse ; limit how small/big it can get
      sigma < 0.1 [set sigma 0.1]
      sigma > 3 [set sigma 3]
    )
    ; mass is proportional to radius (it maybe should be sigma ^ 2 to be proportional to area,
    ;     but atoms don't actually behave that way. Heavier atoms in the same period actually have
    ;     a smaller radius)
    set mass sigma
    aep.set-size
  ]
end


to aep.set-size
  set size sigma * atom-viz-size
end


;; *************************************************************
;; ********** vab (visualize atoms and bonds) procedures********
;; *************************************************************
to vab.setup-links
  set link-check-dist 1.5
  set min-sigma-for-links 0.5
  ask atoms with [size >= min-sigma-for-links] [
    vab.update-links in-radius-linkable-atoms
  ]

  vab.color-links
end

to-report in-radius-linkable-atoms
  report other atoms in-radius cutoff-dist with [sigma >= min-sigma-for-links]
end

;*******************************************************
;**************** Go Procedures ************************
;*******************************************************

to vab.update-atom-color-and-links
  ask atoms [
    vab.update-atom-color atom-PE

    if size >= 0.4 [
      vab.update-links in-radius-linkable-atoms
    ]
  ]
end


;; *****************************************************
;; ************* Atom Display procedures ***************
;; *****************************************************

to vab.update-atom-color [pot-E] ; updating atom color

  ifelse color-atoms-by-PE? [
    set color scale-color color  pot-E -6 0
  ] [
    set color base-color
  ]

end

to vab.update-links [in-radius-atoms] ; updating links
  if show-diagonal-right-links? [
    set heading 330
    vab.link-with-atoms-in-cone in-radius-atoms
  ]
  if show-diagonal-left-links? [
    set heading 30
    vab.link-with-atoms-in-cone in-radius-atoms
  ]
  if show-horizontal-links? [
    set heading 90
    vab.link-with-atoms-in-cone in-radius-atoms
  ]

end

to vab.link-with-atoms-in-cone [atom-set]
  let in-cone-atoms (atom-set in-cone link-check-dist 60)
    if any? in-cone-atoms [
      create-atom-link-with min-one-of in-cone-atoms [distance myself]
    ]
end


to vab.color-links
  ask atom-links [
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
  ]
end

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
10
80
43
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
10
170
43
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
50
175
83
temp
temp
0
.2
0.02
.01
1
NIL
HORIZONTAL

SWITCH
560
150
832
183
color-atoms-by-PE?
color-atoms-by-PE?
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
0
1
-1000

SWITCH
560
80
790
113
show-diagonal-left-links?
show-diagonal-left-links?
0
1
-1000

SWITCH
560
115
790
148
show-horizontal-links?
show-horizontal-links?
0
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
110
180
143
increase-size
aep.change-atom-size .1
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
110
85
143
decrease-size
aep.change-atom-size (- .1)
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
90
190
116
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

PLOT
0
205
185
355
Total PE of system
NIL
NIL
0.0
10.0
-65.0
-60.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum [atom-pe] of atoms"

MONITOR
100
155
185
200
PE of System
sum [atom-pe] of atoms
2
1
11

MONITOR
0
155
90
200
interstitial size
[sigma] of the-interstitial
2
1
11

@#$#@#$#@
## WHAT IS IT?

This is a model showing the effects of interstitial atoms in a crystal. An interstitial atom is one that occupies a lattice site in between the sites of the majority atoms making up the crystal. 

## HOW IT WORKS

The base model is a molecular dynamics model using the Lennard-Jones potential for interatomic forces. But, the size parameter, sigma, is now the average of the sigma for each atom instead of a single paramter for all atoms. This means atoms can be of different sizes and still interact. 

## HOW TO USE IT

Press `setup` and `go` and then use the `decrease-size` and `increase-size` buttons to change the size of the red interstitial atom. 

## THINGS TO NOTICE

Notice how changing the size of the interstitial atom changes the potential energy of the system. 

## THINGS TO TRY

Try making the interstitial atom very big. What happens to the crystal structure and potential energy? Why? 


## RELATED MODELS



## CREDITS AND REFERENCES

This model uses the integration method known as Velocity Verlet. See the [wikipedia entry](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet) for more information

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Kelter, J. and Emery, J. (2023).  NetLogo Substitution Size model.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2022 Uri Wilensky.

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
