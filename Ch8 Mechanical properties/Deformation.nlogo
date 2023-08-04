__includes [ "../nls-files/ch8.nls" "../nls-files/molecular-dynamics-core.nls" "../nls-files/visualize-atoms-and-bonds.nls" ]

breed [atoms atom]

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
]

globals [
  ; these are here because the interface elements don't exist
  temp
  force-mode
  create-dislocation?
  create-floor-and-ceiling?
]


;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set force-mode "Tension"
  set temp .07
  setup-constants
  setup-tension-col
  mdc.setup-atoms-nrc atoms-per-row atoms-per-column
  ask atoms [ch8.init-atom]
  setup-force-mode-shape-and-pinned
  update-lattice-view
  mdc.init-velocity
  ask atom 1 [ mdc.setup-offsets-1sig ]

  vab.setup-links

  setup-dislocation
  setup-force-lines
  setup-floor-and-ceiling

  setup-cross-section
  setup-auto-increment-force
  reset-ticks
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go
  if lattice-view != prev-lattice-view [ update-lattice-view ]
  set auto-increment-force 0
  ask atom-links [ die ]
  ; moving happens before velocity and force update in accordance with velocity verlet
  mdc.move-atoms-die-at-edge
  identify-force-atoms
  update-force-and-velocity-and-PE
  mdc.scale-velocities
  vab.update-atom-color-and-links
  calculate-fl-positions
  vab.color-links  ; stylizing/coloring links
  tick-advance dt
  update-plots
end


; Copyright 2020 Uri Wilensky.
; See Info tab for full copyright and license.
@#$#@#$#@
GRAPHICS-WINDOW
240
10
668
439
-1
-1
20.0
1
10
1
1
1
0
1
0
1
-10
10
-10
10
1
1
1
ticks
30.0

BUTTON
0
85
86
118
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
94
85
179
118
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
205
172
238
force-applied
force-applied
0
2
0.0
.1
1
N
HORIZONTAL

SWITCH
680
65
910
98
color-atoms-by-potential-energy?
color-atoms-by-potential-energy?
1
1
-1000

CHOOSER
680
10
845
55
lattice-view
lattice-view
"large-atoms" "small-atoms" "hide-atoms"
1

SWITCH
680
100
910
133
show-diagonal-right-links?
show-diagonal-right-links?
0
1
-1000

SWITCH
680
135
910
168
show-diagonal-left-links?
show-diagonal-left-links?
0
1
-1000

SWITCH
680
170
910
203
show-horizontal-links?
show-horizontal-links?
0
1
-1000

SLIDER
0
10
172
43
atoms-per-row
atoms-per-row
5
15
9.0
1
1
NIL
HORIZONTAL

SLIDER
0
45
172
78
atoms-per-column
atoms-per-column
5
15
9.0
1
1
NIL
HORIZONTAL

MONITOR
690
220
900
265
external force per forced atom (N)
report-indiv-ex-force
3
1
11

PLOT
680
270
915
420
Stress-Strain Curve - Tension
strain
stress
0.0
0.1
0.0
0.4
false
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "if force-mode = \"Tension\" [ plotxy strain stress ]"

MONITOR
680
425
806
470
total external force (N)
report-total-ex-force
5
1
11

MONITOR
817
425
933
470
sample length (rm)
right-edge - left-fl
5
1
11

TEXTBOX
938
425
1083
495
<- in terms of the \nequilibrium interatomic \ndistance between \ntwo atoms (rm)
11
0.0
1

TEXTBOX
673
480
823
498
NIL
11
0.0
1

TEXTBOX
0
250
150
278
Color Key\nLinks:
11
0.0
1

TEXTBOX
10
275
185
293
high compression: dark red
11
13.0
1

TEXTBOX
10
290
280
308
low compression: light red (+ grey tone)
11
18.0
1

TEXTBOX
9
304
159
322
equilibrium: grey
11
5.0
1

TEXTBOX
9
317
279
345
low tension: light yellow (+ grey tone)
11
0.0
1

TEXTBOX
10
333
170
351
high tension: dark yellow
11
44.0
1

TEXTBOX
0
355
155
373
Atoms:
11
0.0
1

TEXTBOX
10
370
280
388
low potential energy: dark blue
11
103.0
1

TEXTBOX
10
385
295
413
high potential energy: light blue (-> white)
11
107.0
1

TEXTBOX
10
400
270
481
pinned atoms (do not move): black cross\natoms affected by an external force: \nblack dot, near a white line with \narrows on the end
11
0.0
1

SWITCH
5
165
215
198
auto-increment-force?
auto-increment-force?
0
1
-1000

BUTTON
38
127
142
160
delete-atoms
delete-atoms
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

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
