breed [atoms atom]
breed [fl-ends fl-end] ; turtles at the ends of the force lines, point in direction the force is acting
undirected-link-breed [fl-links fl-link] ; force line links
undirected-link-breed [wall-links wall-link] ; force line links
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

  ex-force-applied? ; is an external force directly applied to this atom? False if no, True if yes

  ;; the following variable is for the atom-editing-procedures.nls file
  selected? ; whether the atom is selected or  not to change its size
]

globals [
  temp
  force-mode
  auto-increment-force?
  new-atom-sigma
  create-floor-and-ceiling?

  ;; AEP globals
  prev-atom-viz-size  ; previous atom viz size
  message1 ; this variable holds a turtle for displaying messages.
  message2 ; this variable holds a turtle for displaying messages.
  click-mode

  ;; MP globals
  prev-lattice-view ; the lattice view in the previous time step
  upper-left-fl ; upper left force line - shear
  left-fl ; left force line - compression
  right-fl; right force line - tension
  left-edge; where the right side of the sample is (xcor) - tension
  right-edge ; where the right side of the sample is (xcor) - compression
  orig-length ; original length of sample
  prev-length ; length of sample in previous time step
  unpinned-min-ycor ; min unpinnned ycor for shear
  top-neck-atoms ; agentset of atoms on the top of the neck (thin region) (tension).
                 ; Used in calculating stress
  bottom-neck-atoms ; agentset of atoms on the bottom of the neck (thin region) (tension).
                    ; Used in calculating stress
  num-forced-atoms ; number of atoms receiving external force directly
  auto-increment-force ; force to counteract LJ forces in the x-direction (tension)
  cross-section
  ceiling-ycor
  floor-ycor

  ;; the following global is for the atom-editing-procedures.nls file
  atom-viz-size


  ;; MDC globals
  eps
  eps*4
  cutoff-dist
  dt
  kb

  LJ-force-linear-1sig
  LJ-PE-linear-1sig

  LJ-force-linear-2sig ; since sizes can change, we store the linear constants for different sizes of atom pairs
  LJ-PE-linear-2sig ; since sizes can change, we store the linear constants for different sizes of atom pairs

  ;; viz globals
  link-check-dist ; each atom links with neighbors within this distance
  min-sigma-for-links
]


;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  mp.setup-constants
  set temp 0.001
  set dt .06
  set force-mode "Shear"
  set auto-increment-force? false
  set create-floor-and-ceiling? false
  set click-mode "select-atoms"
  mdc.setup-cutoff-linear-functions-2sig

  mdc.setup-atoms-nrc atoms-per-row atoms-per-column
  ask atoms [
    mp.init-atom
    aep.init-atom
  ]
  mp.setup-force-mode-shape-and-pinned


  mp.update-lattice-view
  mdc.init-velocity


  vab.setup-links
  aep.setup-messages

  mp.setup-dislocation
  mp.setup-force-lines
  mp.identify-force-atoms
  mp.setup-floor-and-ceiling

  mp.setup-cross-section
  mp.setup-auto-increment-force

  reset-ticks
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;



to go
  if lattice-view != prev-lattice-view [ mp.update-lattice-view ]

  ask atom-links [ die ]

  ; moving happens before velocity and force update in accordance with velocity verlet
  mdc.move-atoms-die-at-edge

  mp.identify-force-atoms

  mp.update-force-and-velocity-and-PE-2sig

  mdc.scale-velocities

  vab.update-atom-color-and-links

  mp.calculate-fl-positions

  vab.color-links  ; stylizing/coloring links

  tick-advance dt
  update-plots
end


to interact
  if mouse-down? [
    (ifelse
      click-mode = "drag-atoms" [mdc.drag-atoms-with-mouse-2sig]
      click-mode = "delete-atoms" [aep.delete-atoms]
      click-mode = "add-atoms" [aep.add-atoms new-color new-atom-sigma]
      click-mode = "select-atoms" [aep.select-atoms]
    )
    display
  ]
end



;; AEP procedures

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to aep.init-atom
  set selected? false
  aep.set-size
end

to aep.setup-messages
  crt 1 [
    setxy 6 2.5
    set size 0
    set message1 self
  ]
  crt 1 [
    setxy 6.5 2
    set size 0
    set message2 self
  ]
end

;;*******************************************************
;;**************** Go Procedures ************************
;;*******************************************************

to aep.update-atom-size-viz
  if atom-viz-size != prev-atom-viz-size [
    ask atoms [aep.set-size]
  ]
  set prev-atom-viz-size  atom-viz-size
end


;; *****************************************************
;; *********      Interaction Procedures      **********
;; *****************************************************

; ncolor & nsigma taken as input arguments, they are interface elements in Point Defects
to aep.add-atoms [ncolor nsigma]
  if mouse-down? and not any? atoms with [distancexy mouse-xcor mouse-ycor < .2] [
    let closest-atom min-one-of atoms [distancexy mouse-xcor mouse-ycor]
    let new-atom-force last [mdc.LJ-potential-and-force-2sig (distancexy mouse-xcor mouse-ycor) sigma nsigma] of closest-atom

    ifelse abs new-atom-force < 30 [

      create-atoms 1 [
        mdc.init-atom
        aep.init-atom
        set sigma nsigma
        set mass sigma  ; mass is proportional to radius (only kind of true in vertical direction of periodic table but oh well)
        aep.set-size
        set base-color read-from-string ncolor
        set color base-color
        setxy mouse-xcor mouse-ycor
      ]
      repeat 10 [go]
    ] [
      print "messages"
      ask message1 [set label "Adding that atom there will make things explode!"]
      ask message2 [set label "(if you are very precise it is possible to add an interstitial)"]
      display
      wait 2
      ask message1 [set label ""]
      ask message2 [set label ""]
    ]
  ]
end


to aep.delete-atoms
  if mouse-down? [
    ask atoms with [xcor <= mouse-xcor + .5 and xcor > mouse-xcor - .5
      and ycor <= mouse-ycor + .433 and ycor > mouse-ycor - .433 ] [die]
    display
  ]

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
      sigma < 0.2 [set sigma 0.2]
      sigma > 1.5 [set sigma 1.5]
    )
    ; mass is proportional to radius (it maybe should be sigma ^ 2 to be proportional to area,
    ;     but atoms don't actually behave that way. Heavier atoms in the same period actually have
    ;     a smaller radius)
    set mass 1
    aep.set-size
    set base-color read-from-string new-color
    set color base-color
  ]
end


to aep.set-size
  set size sigma * atom-viz-size
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; MP procedures

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to mp.setup-constants
  set dt .06
  ;set kb 0.0002 ; this is a very small number right now
  ;set eps 0.07
  set kb 0.1 ; just picking a random constant for Kb that makes things work reasonably
  set eps 1
  set eps*4 eps * 4
  set cutoff-dist 2.5 * r-min
  set atom-viz-size .9
end

to mp.setup-tension-col
  ; making a symmetrical sample for tension mode
  if force-mode = "Tension" and atoms-per-column mod 2 = 0 [
    set atoms-per-column atoms-per-column + 1
  ]
end

to mp.setup-force-mode-shape-and-pinned
  let x-dist r-min ; the distance between atoms in the x direction
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)

  ; values used in assigning atom positions
  let ymax max [ycor] of atoms
  let xmax max [xcor] of atoms
  let y-min min [ycor] of atoms
  let x-min min [xcor] of atoms


  (ifelse force-mode = "Shear"[
    set unpinned-min-ycor min [ycor] of atoms with [ycor > y-min + 3.5 * y-dist]
    ask atoms with [
      (
        (xcor = x-min or
          xcor = x-min + (x-dist / 2) or
          xcor = xmax or
          xcor = xmax - (x-dist / 2)
        )
        and
        ( ycor < unpinned-min-ycor)
      )
    ] [
      set pinned? true
      ]
    ask atoms with-min [ycor] [
      set pinned? true
    ]
    ]
    force-mode = "Tension"[
      ask atoms with [xcor = x-min] [die] ; creating the symmetrical shape
      set x-min min [xcor] of atoms
      ; deleting a rectangular area of atoms to create the tension shape
      ask atoms with [
        (ycor >= ymax - 1.25 or ycor <= y-min + 1.25) and
         xcor <= xmax - 2.5 and
         xcor >= x-min + 2.5
      ] [ die ]

      ask atoms with [xcor = x-min or xcor <= x-min + (x-dist / 2) ] [set pinned? true]
      ask atoms with [ xcor <= xmax and xcor >= xmax - (x-dist * 1.75) ][
        set ex-force-applied? true
        set shape "circle-dot"
      ]
      set num-forced-atoms count atoms with [ex-force-applied?]
    ]
    force-mode = "Compression" [
      ask atoms with [xcor = xmax or xcor = xmax - (x-dist / 2) ] [set pinned? true]
    ]
  )

  ask atoms with [pinned?] [ set shape "circle-x"]
end

to mp.init-atom
  set ex-force-applied? false
end

to mp.setup-dislocation
  if create-dislocation? = true [ ; creating the dislocation
    let start-xcor (median [xcor] of atoms)
    let start-ycor (median [ycor] of atoms)
    if force-mode = "Shear" [ set start-ycor unpinned-min-ycor ]
    ask min-one-of atoms [(xcor - start-xcor) ^ 2 + (ycor - start-ycor) ^ 2] [
      set heading -30
      ask atoms in-cone (r-min * atoms-per-column) 1 [die]
    ]
  ]

  ;; old code
;  let x-dist r-min ; the distance between atoms in the x direction
;  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)
;
;  if create-dislocation? = true [ ; creating the dislocation
;    let dis-xcor (median [xcor] of atoms)
;    let dis-ycor (median [ycor] of atoms)
;    if force-mode = "Shear" [ set dis-ycor unpinned-min-ycor ]
;    let ii 0
;    while [ dis-ycor <= ceiling (max [ycor] of atoms) ] [
;      ask atoms with [ ( ycor <= dis-ycor + y-dist * .25 ) and ( ycor >= dis-ycor - y-dist * 0.25 ) ] [
;        if ( xcor <= dis-xcor + x-dist * .25 ) and ( xcor >= dis-xcor - x-dist * 0.25 ) [ die ]
;      ]
;
;      set dis-ycor dis-ycor + y-dist
;      set dis-xcor dis-xcor - x-dist / 2
;      set ii ii + 1
;    ]
;   ]
end

to mp.setup-force-lines
  let x-dist r-min ; the distance between atoms in the x direction
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2)

  ; values used in assigning atom positions
  let ymax max [ycor] of atoms
  let xmax max [xcor] of atoms
  let y-min min [ycor] of atoms
  let x-min min [xcor] of atoms

  (ifelse force-mode = "Tension"  [ ; set up force lines
    create-fl-ends 2
    set right-fl xmax
    set left-edge x-min
    set orig-length right-fl - left-edge
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor right-fl
      set ycor ymax + (y-dist * 2) ]
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor right-fl
      set ycor y-min - (y-dist * 2)
      create-fl-link-with one-of other fl-ends with [xcor = right-fl]]
    ask fl-ends [
      set color white
      set heading 90
    ]
    if force-mode = "Tension" [
      set prev-length orig-length
    ]
  ]
    force-mode = "Compression" [
      create-fl-ends 2
      set left-fl x-min
      set right-edge xmax
      set orig-length right-edge - left-fl
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor left-fl
        set ycor max-pycor - (y-dist * 2) ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor left-fl
        set ycor min-pycor + (y-dist * 2)
        create-fl-link-with one-of other fl-ends with [xcor = left-fl]]
      ask fl-ends [
        set color white
        set heading 90
      ]
    ]
    force-mode = "Shear" [
      create-fl-ends 2
      set upper-left-fl min [xcor] of atoms with [ ycor >= unpinned-min-ycor ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor upper-left-fl
        set ycor ymax + (y-dist * 2) ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor upper-left-fl
        set ycor unpinned-min-ycor
        hide-turtle
        create-fl-link-with one-of other fl-ends]
      ask fl-ends [
        set color white
        set heading 90 ]
  ])
  ask fl-links [
    set color white
  ]
end



to mp.setup-floor-and-ceiling
  if create-floor-and-ceiling? = true [
    set ceiling-ycor max [ycor] of atoms + r-min
    set floor-ycor min [ycor] of atoms - r-min

    crt 1 [
      set xcor min-pxcor
      set ycor ceiling-ycor
      ht
      hatch 1 [
        set xcor max-pxcor

        create-wall-link-with myself [
          set thickness .3
        ]
      ]
    ]

    crt 1 [
      set xcor min-pxcor
      set ycor floor-ycor
      ht
      hatch 1 [
        set xcor max-pxcor
        create-wall-link-with myself [
          set thickness .3
        ]
      ]
    ]
  ]
end

to mp.setup-cross-section
  if force-mode = "Tension" [
    let xcenter (max-pxcor + min-pxcor) / 2
    let top-neck-ycor [ycor] of one-of atoms with [xcor >= xcenter and xcor <= xcenter + r-min ] with-max [ycor]
    let bottom-neck-ycor [ycor] of one-of atoms with [xcor >= xcenter and xcor <= xcenter + r-min ] with-min [ycor]
    set cross-section top-neck-ycor - bottom-neck-ycor
  ]
end

to mp.setup-auto-increment-force
  if force-mode = "Tension" and auto-increment-force? [
    set force-applied 0
    ask atoms with [ ex-force-applied? ]  [
      set vx 0
      set vy 0
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to mp.update-lattice-view
  (ifelse lattice-view = "large-atoms" [
    set atom-viz-size .9
    ask atoms [
      show-turtle
      set size .9
    ]
  ]
  lattice-view = "small-atoms" [
      set atom-viz-size .6
    ask atoms [
       show-turtle
       set size .6
    ]
  ]
  [; lattice-view = hide-atoms
      ask atoms [ hide-turtle ]
  ])
  set prev-lattice-view lattice-view
end

to mp.calculate-fl-positions ; (calculate new force line positions)
  (ifelse force-mode = "Shear" [
      set upper-left-fl min [xcor] of atoms with [ ycor >= unpinned-min-ycor ]
      ask fl-ends [ set xcor upper-left-fl]
    ]
    force-mode = "Tension" [
      set right-fl max [xcor] of atoms
      ask fl-ends with [xcor > 0] [ set xcor right-fl ]
    ]
    force-mode = "Compression" [
      set left-fl min [xcor] of atoms
      ask fl-ends with [xcor < 0] [ set xcor left-fl ]
    ]
    )
  ifelse (force-applied + auto-increment-force) = 0 [
    ask fl-ends [ hide-turtle ]
    ask fl-links [ hide-link ]
  ]
  [
    ask fl-ends [ show-turtle ]
    ask fl-links [ show-link ]
  ]
end

; find the atoms closest to the force line that will be the ones receiving the external force
to mp.identify-force-atoms
  (ifelse force-mode = "Shear" [
    ask atoms [ set ex-force-applied?  false ]
    let forced-atoms atoms with [ ycor >= unpinned-min-ycor and (distancexy upper-left-fl ycor) <= r-min]
    set num-forced-atoms count forced-atoms
    ask forced-atoms [
      set ex-force-applied?  true
      set shape "circle-dot"
    ]
    ]
    force-mode = "Compression" [
      ask atoms [ set ex-force-applied?  false ]
      let forced-atoms atoms with [ (distancexy left-fl ycor) <= r-min]
      set num-forced-atoms count forced-atoms
      ask forced-atoms [
        set ex-force-applied?  true
    ]
  ]) ; for tension, the same atoms in the left shoulder of the sample always receive the force
end

to mp.update-force-and-velocity-and-PE
  ask atoms [
    let [n-fx n-fy sum-PE] mdc.calculate-force-and-velocity-and-PE-1sig

    set n-fy (n-fy + mp.ceiling-or-floor-force)

    set [n-fx n-fy] mp.update-external-force n-fx n-fy

    mdc.update-force-and-velocity n-fx n-fy
  ]
end


to mp.update-force-and-velocity-and-PE-2sig
  ask atoms [
    let [n-fx n-fy sum-PE] mdc.calculate-force-and-velocity-and-PE-2sig

    set n-fy (n-fy + mp.ceiling-or-floor-force)

    set [n-fx n-fy] mp.update-external-force n-fx n-fy

    mdc.update-force-and-velocity n-fx n-fy
  ]
end


to-report mp.ceiling-or-floor-force
  let force 0
  if create-floor-and-ceiling? = true [
    ;; apply force in y-direction due to floor and ceiling
    (ifelse
      ycor > (ceiling-ycor - (2 * r-min)) [  ; ceiling force
        set force item 1 (mdc.LJ-potential-and-force-1sig (ceiling-ycor - ycor))
      ]
      ycor < (floor-ycor + (2 * r-min)) [  ; floor force
        set force item 1 (mdc.LJ-potential-and-force-1sig (ycor - ceiling-ycor))
      ]
    )
  ]

  report force

end

to-report mp.update-external-force [ n-fx n-fy ]
  ;; atom procedure. n-fx and n-fy are the new forces the atoms feels this tick
  if not pinned? [
      ; adjust the forces to account for any external applied forces

    let ex-force 0  ; external force is 0 by default

      if ex-force-applied? [
        ifelse force-mode = "Tension" and auto-increment-force? [
          set ex-force ( - n-fx + 0.001 )  ; In tension, we want ex-force to cancel out any force plus a little more
          set auto-increment-force auto-increment-force + ex-force ; auto-increment-force reports total ex-force in the x-direction
          set n-fy 0 ; in tension, we want atoms with external force to not have any force in y-direction to prevent movement in that direction
        ] [
        ; if force-mode is not tension and this is a forced atom, simply apply external-force
          set shape "circle-dot"
          set ex-force ( force-applied / num-forced-atoms )
        ]
      ]
      if shape = "circle-dot" and not ex-force-applied? [ set shape "circle" ]
      set n-fx ex-force + n-fx
    ]
  report list n-fx n-fy
end

to-report mp.strain ; tension only
  report ((right-fl - left-edge) - orig-length) / orig-length
end

to-report mp.stress ; tension only
  report (((-1 * force-applied) + auto-increment-force) / cross-section)
end

to-report mp.report-indiv-ex-force
  report (force-applied + auto-increment-force) / num-forced-atoms
end

to-report mp.report-total-ex-force
  report force-applied + auto-increment-force
end




;; MDC procedures
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
    let [n-fx n-fy sum-PE] mdc.calculate-force-and-velocity-and-PE-1sig

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
    let [n-fx n-fy sum-PE] mdc.calculate-force-and-velocity-and-PE-2sig
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

to-report mdc.calc-pair-PE-and-force-with-2sig [other-atom]
  report mdc.LJ-potential-and-force-2sig (distance other-atom) sigma  [sigma] of other-atom
end


;; velocity verlet used by both

to-report mdc.velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report mdc.velocity-verlet-velocity [v a new-a]  ; velocity, previous acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end


;; VIZ procedures
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


;; *****************************************************
;; ************* Link Display procedures ***************
;; *****************************************************

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
@#$#@#$#@
GRAPHICS-WINDOW
220
10
693
484
-1
-1
27.4
1
10
1
1
1
0
0
0
1
-8
8
-8
8
1
1
1
ticks
30.0

BUTTON
0
130
100
163
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
110
130
210
163
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
400
210
433
force-applied
force-applied
0
30
0.0
.1
1
N
HORIZONTAL

SWITCH
0
10
210
43
create-dislocation?
create-dislocation?
0
1
-1000

SWITCH
695
60
920
93
color-atoms-by-PE?
color-atoms-by-PE?
1
1
-1000

SWITCH
695
100
925
133
show-diagonal-right-links?
show-diagonal-right-links?
0
1
-1000

SWITCH
695
140
925
173
show-diagonal-left-links?
show-diagonal-left-links?
1
1
-1000

SWITCH
695
180
925
213
show-horizontal-links?
show-horizontal-links?
1
1
-1000

SLIDER
0
45
210
78
atoms-per-row
atoms-per-row
5
20
11.0
1
1
NIL
HORIZONTAL

SLIDER
0
80
210
113
atoms-per-column
atoms-per-column
5
20
10.0
1
1
NIL
HORIZONTAL

MONITOR
0
435
210
480
external force per forced atom (N)
mp.report-indiv-ex-force
3
1
11

TEXTBOX
695
220
845
238
NIL
11
0.0
1

TEXTBOX
695
220
845
248
Color Key\nLinks:
11
0.0
1

TEXTBOX
705
245
880
263
high compression: dark red
11
13.0
1

TEXTBOX
705
260
975
278
low compression: light red (+ grey tone)
11
18.0
1

TEXTBOX
704
274
854
292
equilibrium: grey
11
5.0
1

TEXTBOX
704
287
974
315
low tension: light yellow (+ grey tone)
11
0.0
1

TEXTBOX
705
303
865
321
high tension: dark yellow
11
44.0
1

TEXTBOX
695
320
850
338
Atoms:
11
0.0
1

TEXTBOX
705
335
975
353
low potential energy: dark blue
11
103.0
1

TEXTBOX
705
350
990
378
high potential energy: light blue (-> white)
11
107.0
1

TEXTBOX
705
365
965
446
pinned atoms (do not move): black cross\natoms affected by an external force: \nblack dot, near a white line with \narrows on the end
11
0.0
1

CHOOSER
45
295
155
340
new-color
new-color
"red" "violet" "green" "orange" "blue"
1

BUTTON
0
255
95
288
decrease-size
aep.change-atom-size (- .1)\nrepeat 5 [go]
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
100
255
212
288
increase-size
aep.change-atom-size .1\nrepeat 5 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
695
10
920
55
lattice-view
lattice-view
"small-atoms" "large-atoms" "hide-atoms"
1

BUTTON
0
205
85
245
edit atoms
interact
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

TEXTBOX
90
205
230
255
Click atoms to select them, then change their sizes with the buttons.
11
0.0
1

TEXTBOX
40
385
190
403
Apply an external force
11
0.0
1

TEXTBOX
0
175
215
201
----------------------------------
11
0.0
1

TEXTBOX
0
355
205
381
----------------------------------
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model allows the user to observe the effects of external forces on a close-packed 2D crystal lattice. It also gives a qualitative image of stress/strain fields around edge dislocations. An edge dislocation is a type of crystal defect that is comprised of an extra half-plane inserted into a lattice. They are an important feature of almost all materials because they facilitate deformation. In this model, an edge dislocation can be initialized within the material, and shear, tension, or compression forces can be applied to the material.

This is a Molecular Dynamics simulation, meaning an interatomic potential governs the energy and motion of each atom. Here, we are using the Lennard-Jones potential. Atoms attempt to minimize their potential energy by moving to an equilibrium distance away from the atoms surrounding them, which means that the atoms are moved in response to the force they feel from the surrounding atoms. The position and velocity of the atoms are updated every time step with the velocity Verlet algorithm.

## HOW IT WORKS

### Interatomic Forces

The Lennard-Jones potential shows that there is an equilibrium distance between two atoms where the potential energy of each atom is minimized. If the atoms are farther than this equilibrium distance, they will attract each other. If they are closer, they will repel each other. The potential is given by:

V = 4 * eps ((sigma / r)<sup>12</sup> - (sigma / r)<sup>6</sup>)

In this formula:
- "V" is the potential between two particles.
- "eps", or epsilon, is the depth of the potential well.
- "sigma" is the distance where the inter-particle potential is zero.
- "r" is the distance between the centers of the two particles.

Taking the derivative of this potential equation yields an equation for force. A negative force means the atoms repel each other and a positive force means they attract.

In this model, each atom calculates the force felt from its neighbors in a 5 unit radius. This is to reduce the number of calculations and thus allow the simulation to run faster; after this distance, the force is extremely weak, and therefore it is acceptable to ignore interactions beyond this point. These forces are summed and then divided by the particle's mass to yield an acceleration, which is passed into the velocity Verlet algorithm. This algorithm is a numerical method that integrates Newton's equations of motion. Both the velocity of an atom and the position of an atom after one time step are calculated.

### External Forces

External forces are applied by adding a term to the total LJ force felt from neighboring atoms. This force varies based on the FORCE-MODE. Atoms acted on directly by an external force are marked with a black "X". Pinned atoms, or atoms that do not move, are marked with a black dot. These atoms are excluded from updating position and velocity via the velocity Verlet algorithm.

In SHEAR mode, F-APP will be exerted left-to-right on atoms in the upper left of the material, while the bordering atoms on the bottom half of the crystal are pinned and do not move. The pinning of the bottom atoms represents the bottom of the material being clamped or otherwise held in place, so F-APP cannot simply rotate the lattice.

In TENSION mode, F-APP is exerted leftwards on the left side of the sample. The sample shape is different in tension mode to replicate the samples used in tensile testing. The atoms in the left shoulder are all acted on directly by the external force. If AUTO-INCREMENT-TENSION? is turned on, a force to "counteract" or "equalize" the forces from Lennard-Jones interactions is first applied to shoulder atoms. Then, an additional leftwards force is applied. This creates a smooth stress-strain curve and allows for deformation to occur in the neck region instead of the shoulder region. Shoulder region deformation would be negligible in a real life tensile test. In this simulation, the atoms farthest to the right (border of the right shoulder) are pinned. Note, this is the only mode in which the STRESS - STRAIN CURVE is plotted.

In COMPRESSION mode, F-APP is exerted rightwards on the left side of the sample. The bordering atoms on the right side are pinned.

To help visualize external forces exerted, the white line bookended by two arrows (or one, if in shear mode) shows where the external force is applied. The arrows point in the direction the force acts. In SHEAR and COMPRESSION, atoms within 1 unit of the force line are affected.

In this simulation, the temperature is controlled in order to reduce erratic motion from the atoms. This is done by calculating the thermal speed of the material, which entails finding the mean speed of the atoms. The thermal speed is related to temperature, so the thermal speed using the user-set system temperature is calculated. This is defined as the target speed and the current thermal speed of the material is the current speed. The target speed is divided by the current speed and the velocity of each atom is multiplied by this scaling factor.

## HOW TO USE IT

Before clicking the SETUP button, you first need to select the FORCE-MODE (SHEAR / TENSION / COMPRESSION), whether you’d like a dislocation initialized (CREATE-DISLOCATION?), and the size of the crystal (ATOMS-PER-ROW / ATOMS-PER-COLUMN). If you are in the TENSION FORCE-MODE, you will also need to select whether you'd like the force to be automatically increased or you'd like to control the force manually (AUTO-INCREMENT-TENSION?).

FORCE-MODE allows the user to apply shear, tension, or compression forces to the crystal. Refer back to the External Forces subsection of the HOW IT WORKS section for a description of how the forces are applied.

AUTO-INCREMENT-TENSION? (TENSION mode only) will automatically increase F-APP by .0005 N if the current sample length is equal to the sample length during the previous time step or if the current sample length is shorter than the sample length during the previous time step. This is useful for producing a stress-strain curve. You can start at F-APP = 0 N and allow the simulation to run until fracture occurs.

CREATE-DISLOCATION? allows the user to initialize an edge dislocation in the top half of the lattice. It is not how the dislocation would exist within the material in an equilibrium state. Although the dislocation may exist in a metastable state within this simulation, that is due to the pinning of atoms. In an unpinned material, the dislocation would usually propagate out in lattices of these nano sizes due to surface effects and the Lennard-Jones potential governing atomic motion (however, this propagation is temperature dependent so at low enough temperatures it may be possible to form a metastable dislocation).

ATOMS-PER-ROW sets the number of atoms horizontally in the lattice, while ATOMS-PER-COLUMN sets the number of atoms vertically in the lattice.

These are all of the settings that need to be selected prior to SETUP. The other settings can also be selected before SETUP, but they are able to be changed mid-simulation, while the aforementioned settings are not. The functions of the other settings are listed below.

LATTICE-VIEW provides three options for observing the crystal lattice. LARGE-ATOMS shows atoms that nearly touch each other in equilibrium which helps to visualize close packing. SMALL-ATOMS shows atoms with a reduced diameter which can be used with links to see both atomic movement and regions of tension and compression in the crystal. HIDE-ATOMS allows the user to hide the atoms completely and only use the links to visualize deformation within the material.

SYSTEM-TEMP sets the temperature of the system.

F-APP is the external applied force on the material. It is the total force applied, not the individual force on each atom. The actual numbers are arbitrary, since the Lennard-Jones force has not been calibrated to model a real material. The units are in Newtons.

DELETE-ATOMS allows the user to delete atoms by clicking on them in the display. If the button is pressed, clicking on the View will delete atoms. If it is not pressed, clicking will do nothing.

UPDATE-ATOM-COLOR? controls whether the color of the atoms is updated. The color is an indicator of the potential energy of each atom from Lennard-Jones interactions. A lighter color means the atom has a higher potential whereas a darker color indicates a lower potential.

DIAGONAL-RIGHT-LINKS?, DIAGONAL-LEFT-LINKS?, and HORIZONTAL-LINKS? all provide additional ways to visualize what is happening in the lattice. They are not meant to represent bonds between atoms. The options controlling DIAGONAL-RIGHT-LINKS? and DIAGONAL-LEFT-LINKS?  are particularly useful for identifying where the extra half plane is located in the plane. The links are colored based on their deviation from an equilibrium range of lengths. If the link (the distance between two atoms) is shorter than the equilibrium range, the link will be colored a shade of red and the atoms are compressed in this direction. If the link is longer than the equilibrium range, the link will be colored a shade of yellow. If a link is within the equilibrium range, it is colored grey. See the Color Key on the interface for a more thorough explanation.

The monitor EXTERNAL FORCE PER FORCED ATOM displays the individual force that each atom directly being influenced by the external force is experiencing (in the case of TENSION mode when AUTO-INCREMENT-TENSION? is on, this is the average individual force). CURRENT F-APP is the total external force on the sample (in SHEAR, COMPRESSION, and TENSION without AUTO-INCREMENT-TENSION?, this is the same as the F-APP slider value). SAMPLE LENGTH is the length of the sample from end to end. It is in terms of the equilibrium interatomic distance between two atoms (rm). The stress-strain curve only works for TENSION mode and can be used with AUTO-INCREMENT-TENSION?.

## THINGS TO NOTICE

SHEAR Mode

As the material deforms, how does the edge dislocation travel?

Where are the areas of tension and compression around the edge dislocation?

TENSION Mode

When the material deforms, does it do so randomly or are there observable preferences for deformation in the material?

How does the stress-strain curve correspond to elastic deformation within the material?

How does changing the temperature affect the deformation patterns within the material?

COMPRESSION Mode

When the material deforms, are there sections of atoms that maintain their original shape? Or is the deformation completely random?

## THINGS TO TRY

In SHEAR mode, initialize an edge dislocation and apply the smallest force possible to deform the material. Does the material continue to deform after the edge dislocation propagates out? What force is require to deform the material after the initial edge dislocation has propagated out?

In TENSION mode, samples with larger numbers of atoms per row and smaller numbers of atoms per column are generally best for observing deformation/fracture (For example, 18 atoms per row and 13 atoms per column). Increasing the number of atoms per column also works well; you just want to maintain a sample with a long "neck" area. To produce a stress-strain curve, set F-APP to 0 N and turn AUTO-INCREMENT-TENSION? on. While running the simulation at a lower temperature creates a smoother stress-strain curve, the slip behavior differs for low and high temperatures, so it is worthwhile to run the simulation with both.

While running the simulation, pay attention to the different directions of links. Are tension and compression concentrated in certain areas? Do they differ in different directions?

Create vacancies with the sample in tension and compression by deleting atoms. Does this change how the material deforms?

Create a second edge dislocation in the shear mode by deleting atoms. How does this change material deformation?

## EXTENDING THE MODEL

Add a slider to vary eps (epsilon). How does changing eps affect deformation? Are smaller or larger forces needed to deform the material as eps increases? Does the material observably deform differently?

Color the atoms according to a different property than their potential energy. Suggestions include according to the magnitude of force felt or the direction of the net force on each atom.

Apply forces in different directions than the ones provided. Does the material deform in the same way? Why or why not?

(Advanced) Add in another type of atom. Give this atom different properties, such as a different mass, or a different radius. You will need to store a separate eps and sigma for each atom and their different types of interactions. For example, if you have A and B atoms, you will need to have an A-A eps and sigma, an A-B eps and sigma, and a B-B eps and sigma.

## RELATED MODELS

The Lennard-Jones model

## NETLOGO FEATURES

When a particle moves off of the edge of the world, it doesn't re-appear by wrapping onto the other side (as in most other NetLogo models). We would use this world wrapping feature to create periodic boundary conditions if we wanted to model a bulk material.

Link coloring uses the transparency feature in order to create greyish shades. In order to make links in tension (yellow) and compression (red) that are very close to equilibrium (grey) visually close to equilibrium, they are more transparent in order to produce a grey hue. The farther away a link is from equilibrium, the more opaque it is.

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Cetin, S.,  Kelter, J. and Wilensky, U. (2020).  NetLogo Dislocation Motion and Deformation model.  http://ccl.northwestern.edu/netlogo/models/DislocationMotionandDeformation.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2020 Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

<!-- 2020 Cite: Cetin, S.,  Kelter, J. -->
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

circle-+
false
0
Circle -7500403 true true 0 0 300
Rectangle -16777216 true false 0 120 315 165
Rectangle -16777216 true false 135 -15 180 300

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
