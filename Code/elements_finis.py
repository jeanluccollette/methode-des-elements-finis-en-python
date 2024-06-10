import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve


def coeffDK(trs, pts):
    liste_Dp = np.zeros((trs.shape[0], 3, 3))
    liste_Kp = np.zeros((trs.shape[0]))
    Dphi = np.array([[-1, -1], [1, 0], [0, 1]])
    for tr in range(trs.shape[0]):
        xp = pts[trs[tr, :], 0]
        yp = pts[trs[tr, :], 1]
        detJp = (xp[1]-xp[0])*(yp[2]-yp[0]) - (xp[2]-xp[0])*(yp[1]-yp[0])
        Kp = np.abs(detJp)/2
        Bp = np.array([[yp[2]-yp[0], yp[0]-yp[1]],
                      [xp[0]-xp[2], xp[1]-xp[0]]])/detJp
        liste_Dp[tr, :, :] = Kp*Dphi@Bp.T@Bp@Dphi.T
        liste_Kp[tr] = Kp
    return liste_Dp, liste_Kp


def solve_edp(triangles, points, labels, a, f, cond_lim):

    print('Génération des listes')
    liste_Dp, liste_Kp = coeffDK(triangles, points)

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
        x = np.sum(points[triangle, 0])/3
        y = np.sum(points[triangle, 1])/3
        a_xy = a(x, y)
        f_xy = f(x, y)
        Dp = liste_Dp[n_tri, :, :]
        Kp = liste_Kp[n_tri]
        s_tri = 0
        for noeud in triangle:
            if labels[noeud] == 0:
                k = index_k[noeud]
                B[k] = B[k] + f_xy*Kp/3
                ligne.append(k)
                colonne.append(k)
                valeur.append(Dp[s_tri, s_tri] + a_xy*Kp/6)
                for s_vois in [-1, 1]:
                    noeud_voisin = triangle[(s_tri+s_vois) % 3]
                    k_vois = index_k[noeud_voisin]
                    val = Dp[s_tri, (s_tri+s_vois) % 3] + a_xy*Kp/12
                    if labels[noeud_voisin] == 0:
                        ligne.append(k)
                        colonne.append(k_vois)
                        valeur.append(val)
                    for clim in range(len(cond_lim)):
                        if labels[noeud_voisin] == cond_lim[clim][0]:
                            B[k] = B[k] - cond_lim[clim][1]*val
            s_tri = s_tri+1
    A = csr_array((valeur, (ligne, colonne)),
                  shape=(len(inconnus), len(inconnus)))

    print('Résolution du système linéaire')
    sol = spsolve(A, B)

    return sol
