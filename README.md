# Practical 2: Constraints

## Handout date: 8/Mar/2022.

## Deadline: 17/Mar/2022 08:59.

The second practical generalizes and extends the first practical, by working with constraint-based velocity and position resolution. We still only work with rigid bodies. The objectives of the practical are:

1. Implement impulse-based velocity resolution by constraints (Lecture 6).
 <br />
2. Implement position correction by constraints (Lecture 7).
 <br />
3. Generalize collision resolution to work with the implemented constraint-resolution framework.
 <br />
4. Extend the framework with some chosen effects.  

This is the repository for the skeleton on which you will build your second practical. Using CMake allows you to work and submit your code in all platforms. The entire environment is in C++, but most of the "nastier" coding parts have been simplified; for the most part, you only code the mathemtical-physical parts. The environment is otherwise identical to the $1^{st}$ practical, including all compilation instructions.


## Background

The practical is mostly about implementing procedures for velocity and position resolution to satisfy constraints. This resolution integrates into the game-physics engine as follows:

In each scene iteration:

1. Integrate (Practical 1)
2. Detect collisions
3. Resolve Collisions (Practical 2 generalization - using constraints rather then direct particular solution)
4. Correct velocities to be tangent to user constraints
5. Fix positions to satisfy user constraints

The constraints we deal with in this practical are purely *holonomic* and bivariate. That is, each constraint $C$ is of the form $C(x_1,x_2)$, where $x_1$ and $x_2$ are two positions on two meshes (can be the same mesh). We distinguish between user constraints, read from file, and collision constraints, created on-the-fly in each iteration. Moreover, we distinguish between equality constraints $C=0$ and inequality constraints $C \geq 0$; the latter only matters when they are *violated*, and otherwise should not affect velocities or positions. In practice, we use some tolerance $\tau$ that measures validity (opting for $|C|<\tau$ for equality constraints), rather then adhere to perfect $0$ which is unattainable numerically.


### Working with Rigid Bodies

The constraints are expressed using any two points on a body (which happen to be vertices in user constraints); nevertheless, the bodies are rigid, and therefore the only movement degrees of freedom for a mesh are its COM position $p$, orientation quaternion $q$, linear COM velocity $v$ and angular velocity $\omega$. Resolving constraints should only work and change these variables, and not touch any individual vertex. Specifically, never alter ```currV``` directly; rather, recompute it from ```origV``` in the end of a time-step iteration after having corrected $p$ and $q$.

### Velocity Resolution

For equality constraints, the total velocities $\overline{v}_1$ and $\overline{v}_2$ should always satisfy $Jv=0$, where $J$ is the gradient of the constraint, and $v$ is a vector comprising $v_1, \omega_1, v_2,\omega_2$ in order (check: vector length is $12$ variables). If $Jv \neq 0$, You will be computing $\Delta v$ to satisfy $J(v+\Delta v)=0$, using the Lagrange multiplier method learnt in class (Lecture 6). This requires setting up an (inverse) mass matrix of $12 \times 12$, with the body masses and (inverse) inertia tensors in order. Use $0$ for inverse mass and inverse inertia tensor for fixed bodies, which will simulate the correct effect (as if they have infinite resistance). Note that the inertia tensor should transform like in the first practical; essentially your constraint-based collision resolution should be almost equivalent to what you implemented explicitly before, just in this more general framework (*do not* implement a particual solution to collision that is different from the constraints mechanism!).

Note: the part in the mass matrix corresponding to the linear velocity has the scalar masses $m_1$ and $m_2$ repeated $3$ times each in the diagonal of the matrix, for the $x,y,z$ components of the respective velocities.

The coefficient of restitution is given for collisions constraints in order to induce elastic velocity bias; you should use it as instructed in class (user constraints set it to $0$ by default).

The user constraints that are read from file attach two vertices from two meshes in a distance that has to be maintained. That is, the constraint is $C(x_1,x_2) = \left|x_1-x_2\right| - d_{12}$, where $d_{12}$ is computed for the position at time $t=0$. You should derive $J$ for that constraint on your own (you have a hint for it in Lecture 7; for intuition, you are supposed to get that the velocities of both vertices should only move in a way that doesn't change this distance, like it's a fixed rod).

### Position Correction

Position correction is similar to velocity correction, except that we take the easy route (in the basic practical requirements), and only correct *linearly*. That is, we do not change $q$, only $p$ of every body. That means the mass matrix is only $6 \times 6$ of body masses, without any inertia tensor components, and the Jacobian only contains derivatives relating to linear movement. That generalizes the linear-interpenetration resolution for collisions. Note that this means totally different $J, M, \lambda$ for this step, which do not relate to those computed in the velocity correction stage! The theoretical details are in lecture 7. We do not employ stiffness in this practical.

See below for details on where to do all that in the code.

### Extensions

The above will earn you $80\%$ of the grade. To get a full $100\%$, you must choose *a single* extension out of these 3 extension options, and augment the practical with it. Some choices will require minor adaptations to the GUI or the function structure which are easy to do. The extension will earn you $20\%$, and the exact grading will commensurate with the difficulty. Note that this means that all extensions are equal in grade; if you take on a hard extension, it's your own challenge to complete it well.

1. Make the fixed-distance user constraint more flexible to some extent, and therefore a two-sided inequality constraint (for instance, the fixed rod could then compress or stretch up to $20\%$ from the original $d_{12}$). **Level: easy**.

2. Fix the linear position-correction hack by adding $q$ orientation correction to constraints. For this you will need the derivatives of the position of a point w.r.t. $q$, which are not trivial; look [here](http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf) for inspiration. **Level: intermediate-hard**.

3. Add another type of original constraint, which has to be concretely exemplified. For instance, bending or some limitation on rotation. **level: intermediate**.

You may invent your own extension as substitute to **one** in the list above, but it needs approval on the Lecturer's behalf **beforehand**.


## Installation

*This installation is exactly like that of Practical 1, repeated here for completeness*

The skeleton uses the following dependencies: [libigl](http://libigl.github.io/libigl/), and consequently [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), for the representation and viewing of geometry, and [libccd](https://github.com/danfis/libccd) for collision detection. libigl viewer is using [dear imGui](https://github.com/ocornut/imgui) for the menu. Everything is bundled as either submodules, automatically installs with cmake, and you do not have to take care of any installation details. To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/INFOMGP-Practical2.git
```

to compile the environment, go into the `practical2` folder and enter in a terminal (macOS/Linux):

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

In windows, you need to use [cmake-gui](https://cmake.org/runningcmake/). Pressing twice ``configure`` and then ``generate`` will generate a Visual Studio solution in which you can work. The active soution should be ``practical2_bin``. *Note*: it only seems to work in 64-bit mode. 32-bit mode might give alignment errors.

### Using the dependencies

You do not need to acquaint yourself much with any dependency, nor install anything auxiliary not mentioned above. For the most part, the dependencies are parts of code that are background, or collision detection code, which is not a direct part of the practical. The most significant exception is [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for the representation and manipulation of vectors and matrices. However, it is a quite a shallow learning curve. It is generally possible to learn most necessary aspects (multiplication of matrices and vectors, intialization, etc.) just by looking at the existing code. However, it is advised to go through the "getting started" section on the Eigen website (reading up to and including "Dense matrix and array manipulation" should be enough).

### Working with the repository

All the code you need to update is in the ``practical2`` folder. Please do not attempt to commit any changes to the repository. <span style="color:red">You may ONLY fork the repository for your convenience and work on it if you can somehow make the forked repository PRIVATE afterwards</span>. Solutions to the practical which are published online in any public manner will **disqualify** the students! submission will be done in the "classical" department style of submission servers, published separately.

## The coding environment for the tasks

You will find the environment almost identical to the first practical, with these main differences:

1. A constraint file is being read, and has to be given as the third argument to the executable.

2. The ```handleCollision()``` needs to be written by you to work with constraints; see comments within.

3. The ```updateScene()``` function is already written down to work with the game-engine loop.

4. Most of the work is in ```Constraints.h```, where you have to fill in the velocity and position correction functions. This implements a class that is invokes from the scene class.

The code you have to complete is always marked as:

```cpp
/***************
TODO
***************/
```

The description of the function will tell you what exactly you need to put in. In some functions, you will have to complete parts you already did in the first practical (to avoid "spoilers")---it's a simple copy and paste (if you did it correctly the last time).

### Input

The TXT file that describes the scene, where you have several examples in the`data` subfolder, is the same. For completeness, the format of the file is:

```
#num_objects
object1.mesh  density1  youngModulus1 PoissonRatio1 is_fixed1    COM1     q1
object2.mesh  density2  youngModulus2 PoissonRatio2 is_fixed2    COM2     q2
.....
```

Where:

1. ``objectX.mesh`` - an MESH file (automatically assumed in the `data` subfolder) describing the geometry of a tetrahedral mesh. The original coordinates are translated automatically to have $(0,0,0)$ as their COM.
<br />
2. ``density`` - the uniform density of the object. The program will automatically compute the total mass by the volume.
<br />
3. ``is_fixed`` - if the object should be immobile (fixed in space) or not.
<br />
4. ``COM`` - the initial position in the world where the object would be translated to. That means, where the COM is at time $t=0$.
<br />
5. ``q`` - the initial orientation of the object, expressed as a quaternion that rotates the geometry to $q*object*q^{-1}$ at time $t=0$.
<br />
6. ``youngModulus1`` and  ``PoissonRatio1`` should be ignored for now; we will use them in the $3^{rd}$ practical.

The user attachment constraints file, given as the third argument, has to have the following format:

```
#num_constraints
mesh_i1 vertex_i1 mesh_j1 vertex_j1 
mesh_i2 vertex_i2 mesh_j2 vertex_j2 
.....
```

Each row is a constraint attaching the vertex ```vertex_i1``` of mesh ```mesh_i1``` to ```vertex_j2``` of mesh ```mesh_j2```. Every row is an independent such constraint. You can find TXT files in the data folder with similar name to the scenes they accompany. You can of course write new ones. Note that the meshes start indexing from $1$---if you put a constraint to mesh $0$, it will get attached to the platform (which should still work).


### User interface

![screenshot of viewer](practical2_interface.png "screenshot of viewer")

The viewer presents the loaded scene, and you may interact with the viewing with the mouse: rotate with the left button pressed and moving around (the "[" and "]" buttons change the behaviour of the trackball), zoom with the mousewheel, and translate with the right button pressed and dragging. Some other options are printed to the output when the program starts.

The menu also controls the visual features, and the setting of the coefficient of restitution and the time step. They can be updated at any point in the simulation. You might add more parameters with some extensions. Everything is set up in `main()`.

The simulation can be run in two modes: continuously, toggled with the `space` key (to stop/run), and step by step, with the `S` key. This behavior is already encoded. The visual update of the scene from the objects is also already encoded.

The main difference is that user attachement constraints are highlighted as yellow cylinders. They are just markers to a constraint and not real physical objects in the scene (so they can collide etc.).

Note that the ```demo``` folder contains compiled demos for windows and OSX; they are to be used as inspiration, because every solution can be a bit different (butterfly effect).


## Submission

The entire code of the practical has to be submitted in a zip file to the lecturer by E-mail. The deadline is **17/Mar/2022 08:59**. Late submissions will not be acceptable unless approved explicitly.

The practical must be done **in pairs**. Doing it alone requires a-priori permission. Any other combination (more than 2 people, or any number not in $\mathbb{N}$) is not allowed.

The practical will be checked in class during slots exactly as in Practical 1. Come prepared with a release version, and make sure it runs on all "-scene" and "-scene-constraints" pairs in the data folder (at least attempt to). The public sheet is [here](https://docs.google.com/spreadsheets/d/1Zeo2mPu_wv9xjMrAFdzhcTyMGLOpsJqtaf_4Y4oSp2I/edit#gid=2004257992).

##Frequently Asked Questions

Here are detailed answers to common questions. Please read through whenever ou have a problem, since in most cases someone else would have had it as well.

<span style="color:blue">Q:</span> I am getting "alignment" errors when compiling in Windows.
<span style="color:blue">A:</span> Delete everything, and re-install using 64-bit configuration in `cmake-gui` from a fresh copy. If you find it doesn't work out-of-the-box, contact the Lecturer. Do not install other non-related things, or try to alter the cmake.

<span style="color:blue">Q:</span> Why is the demo not working out of the box?
<span style="color:blue">A:</span>: with the same parameters as your input program: practical2_bin "folder_name_without_slash" "name of txt scene files".


# Success!
