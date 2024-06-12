# Méthode des éléments finis en Python

## Motivation
L'ensemble des résolutions numériques d'équations aux dérivées partielles, effectuées avec les différents programmes proposés ici, est loin d'être exhaustif. On se limite par ailleurs à des problèmes à deux dimensions, des conditions aux limites et des géométries particulières.
Cependant, la fonction **solve_edp** dans le fichier **elements_finis.py** (voir dans le dossier **Code**), utilisée dans tous les exemples traités, ne nécessite que les packages mentionnés ci-dessous.
```
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve
```
Il est ainsi possible de tester facilement la méthode des éléments finis, sans recourir à des environnements très complets et très performants comme **FEniCS** (https://fenicsproject.org/), mais dont l'installation et la prise en main sont tout de même assez complexes.
Pour le maillage, quelques méthodes "ad hoc" sont proposées. La mise en oeuvre de ces méthodes est facilitée par les fonctions de triangulation disponibles dans **matplotlib.tri** (triangulation de Delaunay, par exemple). Certains exemples font aussi appel au générateur de maillage **GMSH**, dont la version en **Python** est aisément installable via la commande ci-dessous (https://pypi.org/project/gmsh/).
```
pip install gmsh
```
D'autres maillages avec d'autres géométries peuvent être testés, à vous la main :-)

## Description des programmes

Une résolution numérique de l'équation $-\Delta u+au=f$ est proposée dans ce notebook, avec des exemples où le domaine de définition $\Omega$ est la section d'un coaxial (rayon interne r1 et rayon externe r2) et des conditions aux limites de Dirichlet aux frontières (u1 sur le cercle de rayon r1 et u2 sur le cercle de rayon r2).

Pour une comparaison avec les résultats obtenus dans Matlab (toolbox PDE), les données issues de Matlab sont sauvegardées dans les fichiers ci-dessous


*   **nodes.csv** : liste des noeud du maillage avec leurs coordonnées, obtenue dans Matlab
*   **elements.csv** : liste des triangles dont les sommets sont désignés par des index dans la liste des noeuds
*   **solution.csv** : solution de l'approximation obtenue dans Matlab
*   **labels.csv** : labels associés aux noeuds (0 : inconnue, 1 : condition aux limites avec r1, 2 : condition aux limites avec r2)
*   **params.csv** : paramètres r1, r2, u1, u2 (a et f sont des fonctions de x et de y)

La fonction **coaxial_matl** réalise cette comparaison. Les fonctions **coaxial_circ** et **coaxial_hexa** utilisent uniquement Python pour générer le maillage et obtenir l'approximation de la solution. Le fichier archive **domaine_coaxial.zip** (dans le dossier **Code**) contient les programmes permettant de tester ces fonctions.
 
La fonction **solve_edp** dans **elements_finis.py** peut s'appliquer à d'autres géométries que celles données dans les exemples. Cependant, il faudra aussi définir préalablement le maillage, avec la construction des listes **triangles**, **points** et **labels**. En complément de ce choix d'autres géométries, il est par ailleurs possible d'imposer aussi des conditions de Neumann homogènes (voir le fichier **domaine_carre.zip** dans le dossier **Code**, ainsi que les illustrations).

Pour générer automatiquement le maillage, le générateur **GMSH** (https://gmsh.info/) est disponible. La version en Python (https://pypi.org/project/gmsh/) a été testée. Les fichiers archive **domaine_carre_gmsh.zip** et **domaine_coaxial_gmsh.zip** (dans le dossier **Code**) contiennent les programmes permettant de tester la résolution avec le maillage obtenu (voir aussi les illustrations).

## Version Notebook

Le fichier notebook fourni genère les fichiers **CSV** permettant la comparaison avec **Matlab** pour une configuration particulière des paramètres. Dans l'état, il peut être directement utilisable dans **Google Colaboratory**, par exemple.

## Version Code

Les programmes donnés en exemple peuvent être lancés directement dans une console **Python**. Si on dispose de **Matlab**, les fichiers **carre_matl.m** et **coaxial_matl.m** permettent éventuellement de tester d'autres configurations, en générant les fichiers **CSV** associés.

## Site intéressant
Dans le programme **elements_finis.py** (fonction **solve_edp**), les matrices dites de masse et de rigidité sont calculées simultanément. Le détail des calculs est décrit dans le site ci-dessous.

https://bthierry.pages.math.cnrs.fr/course-fem/lecture/elements-finis-triangulaires/contributions-elementaires/

Dans la version PDF ci-dessous, on retrouve le détail des calculs dans le chapitre 3.4 "Calcul des Matrices Élémentaires".

https://bthierry.pages.math.cnrs.fr/course-fem/download/FEM.pdf

## Illustrations
### Comparaison avec Matlab

Le maillage est récupéré à partir des données générées dans Matlab.

![](Images/coaxial_matl_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=2$ et $f(x,y)=10x$.

![](Images/coaxial_matl_sol.png)

On compare la solution calculée avec la fonction **solve_edp** et celle de Matlab, pour valider le programme **Python** fourni.

![](Images/coaxial_matl_err.png)

### Maillage circulaire

Le maillage est généré avec **Python**.

![](Images/coaxial_circ_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=2$ et $f(x,y)=10x$.

![](Images/coaxial_circ_sol.png)

### Maillage hexagonal

Le maillage est généré avec **Python**.

![](Images/coaxial_hexa_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=2$ et $f(x,y)=10x$.

![](Images/coaxial_hexa_sol.png)

## Autres géométries

### Comparaison avec Matlab

Le maillage est récupéré à partir des données générées dans Matlab.

![](Images/carre_matl_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=1$ et $f(x,y)=\dfrac{1}{20}\left(\left(x-\dfrac{D}{2}\right)^2+\left(y-\dfrac{D}{2}\right)^2\right)$.

![](Images/carre_matl_sol.png)

On compare la solution calculée avec la fonction **solve_edp** et celle de Matlab, pour valider le programme **Python** fourni.

![](Images/carre_matl_err.png)

### Maillage à base de carrés

Le maillage est généré avec **Python**.

![](Images/carre_pyth_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=1$ et $f(x,y)=\dfrac{1}{20}\left(\left(x-\dfrac{D}{2}\right)^2+\left(y-\dfrac{D}{2}\right)^2\right)$.

![](Images/carre_pyth_sol.png)

### Maillage du domaine carré avec GMSH

Le maillage est généré avec la version Python de GMSH. Le paramètre **h** de la fonction **carre_gmsh** est la dimension moyenne d'un élément du maillage.

![](Images/carre_gmsh_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=1$ et $f(x,y)=\dfrac{1}{20}\left(\left(x-\dfrac{D}{2}\right)^2+\left(y-\dfrac{D}{2}\right)^2\right)$.

![](Images/carre_gmsh_sol.png)

### Maillage du domaine coaxial avec GMSH

Le maillage est généré avec la version Python de GMSH. Le paramètre **h** de la fonction **coaxial_gmsh** est la dimension moyenne d'un élément du maillage.

![](Images/coaxial_gmsh_maillage.png)

La solution est calculée avec la fonction **solve_edp**.

On impose $a(x,y)=2$ et $f(x,y)=10x$.

![](Images/coaxial_gmsh_sol.png)
