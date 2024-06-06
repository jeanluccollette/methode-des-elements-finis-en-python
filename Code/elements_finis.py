import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve


def jacobiens(triangles, points):
    liste_jacobiens = np.zeros((triangles.shape[0]))
    for tri in range(triangles.shape[0]):
        vec1 = [points[triangles[tri, 1], 0]-points[triangles[tri, 0], 0],
                points[triangles[tri, 1], 1]-points[triangles[tri, 0], 1]]
        vec2 = [points[triangles[tri, 2], 0]-points[triangles[tri, 0], 0],
                points[triangles[tri, 2], 1]-points[triangles[tri, 0], 1]]
        liste_jacobiens[tri] = np.abs(np.linalg.det([vec1, vec2]))
    return liste_jacobiens


def gradients(triangles, points):
    liste_gradients = np.zeros((triangles.shape[0], 3, 3))
    for k in range(triangles.shape[0]):
        a = np.zeros((3, 3))
        for m in range(3):
            a[m, :] = [points[triangles[k, m], 0],
                       points[triangles[k, m], 1], 1]
        for m in range(3):
            b = np.zeros((3, 1))
            b[m] = 1
            sol = np.linalg.solve(a, b)
            liste_gradients[k, m, :] = sol.T
    return liste_gradients


def solve_edp(triangles, points, labels, a, f, cond_lim):

    print('Génération des listes')
    liste_jacobiens = jacobiens(triangles, points)
    liste_gradients = gradients(triangles, points)

    print('Génération des matrices de masse et de rigidité')
    index_k = np.zeros(labels.shape[0], dtype=np.int64)
    inconnus = np.where(labels == 0)[0]
    index_k[inconnus] = np.arange(inconnus.shape[0])
    B = np.zeros(len(inconnus))
    ligne = []
    colonne = []
    valeur = []
    for n_tri in range(triangles.shape[0]):
        triangle = triangles[n_tri]
        jac = liste_jacobiens[n_tri]
        grad = liste_gradients[n_tri, :, 0:2]
        ordre_tri = 0
        for noeud in triangle:
            grad_i = grad[ordre_tri, 0:2]
            k = index_k[noeud]
            if labels[noeud] == 0:
                B[k] = B[k] + f[noeud]*jac/12
                ligne.append(k)
                colonne.append(k)
                valeur.append(jac * (np.dot(grad_i, grad_i)/2+a[noeud]/12))
                for n_vois in [-1, 1]:
                    noeud_voisin = triangle[(ordre_tri+n_vois) % 3]
                    k_vois = index_k[noeud_voisin]
                    grad_j = grad[(ordre_tri+n_vois) % 3, 0:2]
                    val = jac * (np.dot(grad_i, grad_j)/2+a[noeud_voisin]/24)
                    B[k] = B[k] + f[noeud_voisin]*jac/24
                    if labels[noeud_voisin] == 0:
                        ligne.append(k)
                        colonne.append(k_vois)
                        valeur.append(val)
                    for clim in range(len(cond_lim)):
                        if labels[noeud_voisin] == cond_lim[clim][0]:
                            B[k] = B[k] - cond_lim[clim][1]*val
            ordre_tri = ordre_tri+1
    A = csr_array((valeur, (ligne, colonne)),
                  shape=(len(inconnus), len(inconnus)))

    print('Résolution du système linéaire')
    sol = spsolve(A, B)

    return sol
