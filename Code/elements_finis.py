import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve


def coeffDM(trs, pts):
    liste_coeffD = np.zeros((trs.shape[0], 3, 3))
    liste_coeffM = np.zeros((trs.shape[0], 3, 3))
    Dphi = np.array([[-1, -1], [1, 0], [0, 1]])
    for tr in range(trs.shape[0]):
        xp = pts[trs[tr, :], 0]
        yp = pts[trs[tr, :], 1]
        detJp = (xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0])
        Kp = np.abs(detJp)/2
        Bp = np.array([[yp[2]-yp[0], yp[0]-yp[1]],
                      [xp[0]-xp[2], xp[1]-xp[0]]])/detJp
        liste_coeffM[tr, :, :] = (
            Kp/12)*np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
        liste_coeffD[tr, :, :] = Kp*Dphi@Bp.T@Bp@Dphi.T
    return liste_coeffD, liste_coeffM


def solve_edp(triangles, points, labels, a, f, cond_lim):

    print('Génération des listes')
    liste_coeffD, liste_coeffM = coeffDM(triangles, points)

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
        ordre_tri = 0
        for noeud in triangle:
            k = index_k[noeud]
            if labels[noeud] == 0:
                B[k] = B[k] + f*liste_coeffM[n_tri, ordre_tri, ordre_tri]
                ligne.append(k)
                colonne.append(k)
                valeur.append(liste_coeffD[n_tri, ordre_tri, ordre_tri] +
                              a*liste_coeffM[n_tri, ordre_tri, ordre_tri])
                for n_vois in [-1, 1]:
                    noeud_voisin = triangle[(ordre_tri+n_vois) % 3]
                    k_vois = index_k[noeud_voisin]
                    val = liste_coeffD[n_tri, ordre_tri, (
                        ordre_tri+n_vois) % 3]+a*liste_coeffM[n_tri, ordre_tri, (ordre_tri+n_vois) % 3]
                    B[k] = B[k] + f * \
                        liste_coeffM[n_tri, ordre_tri, (ordre_tri+n_vois) % 3]
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
