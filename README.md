# Méthode des éléments finis en Python

## Description des programmes

L'équation résolue numériquement est $-\Delta u+au=f$, avec un domaine de définition $\Omega$ qui est la section d'un coaxial (rayon interne r1 et rayon externe r2) et des conditions aux limites de Dirichlet aux frontières (V1 sur le cercle de rayon r1 et V2 sur le cercle de rayon r2).

Pour une comparaison avec les résultats obtenus dans **Matlab** (toolbox **PDE**), les données issues de **Matlab** sont sauvegardés dans les fichiers ci-dessous


*   **nodes.csv** : liste des noeud du maillage avec leurs coordonnées, obtenue dans Matlab
*   **elements.csv** : liste des triangles dont les sommets sont désignés par des index dans la liste des noeuds
*   **solution.csv** : solution de l'approximation obtenue dans Matlab
*   **labels.csv** : labels associés aux noeuds (0 : inconnue, 1 : condition aux limites avec r1, 2 : condition aux limites avec r2)
*   **params.csv** : paramètres r1, r2, V1, V2, a, f

La fonction **coaxial_matlab** réalise cette comparaison.

Les fonctions **coaxial_circ** et **coaxial_hexa** utilisent uniquement **Python** pour générer le maillage et obtenir l'approximation de la solution.

## Version Notebook

Le fichier notebook fourni genère les fichiers **CSV** permettant la comparaison avec **Matlab** pour une configuration particulière des paramètres. Dans l'état, il peut être directement utilisable dans **Google Colaboratory**, par exemple.

## Version Code

Les programmes de test peuvent être lancés directement dans une console **Python**. Si on dispose de **Matlab**, le fichier **elements_finis.m** permet éventuellement de tester d'autres configurations, en générant les fichiers **CSV** associés.

## Site intéressant
Dans le programme **elements_finis.py** (fonction **sys_lin**), les matrices dites de masse et de rigidité sont calculées simultanément. Le détail des calculs est décrit dans le site ci-dessous.

https://bthierry.pages.math.cnrs.fr/course-fem/lecture/elements-finis-triangulaires/contributions-elementaires/