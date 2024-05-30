import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve


def adj_triangles(triangles, points):
    liste_adj_triangles = []
    for noeud in range(points.shape[0]):
        index_triangles = np.where(np.any(triangles == noeud, axis=1))[0]
        liste_adj_triangles.append(index_triangles)
    return liste_adj_triangles


def jacobiens(triangles, points):
    liste_jacobiens = []
    for tri in range(triangles.shape[0]):
        vec1 = [points[triangles[tri, 1], 0]-points[triangles[tri, 0], 0],
                points[triangles[tri, 1], 1]-points[triangles[tri, 0], 1]]
        vec2 = [points[triangles[tri, 2], 0]-points[triangles[tri, 0], 0],
                points[triangles[tri, 2], 1]-points[triangles[tri, 0], 1]]
        liste_jacobiens.append(np.abs(np.linalg.det([vec1, vec2])))
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


def sys_lin(triangles, labels, a, f, V1, V2, liste_jacobiens, liste_gradients):
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
                B[k] = B[k] + f*jac/6
                ligne.append(k)
                colonne.append(k)
                valeur.append(jac * (np.dot(grad_i, grad_i)/2+a/12))
                for n_vois in [-1, 1]:
                    noeud_voisin = triangle[(ordre_tri+n_vois) % 3]
                    k_vois = index_k[noeud_voisin]
                    grad_j = grad[(ordre_tri+n_vois) % 3, 0:2]
                    val = jac * (np.dot(grad_i, grad_j)/2+a/24)
                    if labels[noeud_voisin] == 0:
                        ligne.append(k)
                        colonne.append(k_vois)
                        valeur.append(val)
                    elif labels[noeud_voisin] == 1:
                        B[k] = B[k] - V1*val
                    elif labels[noeud_voisin] == 2:
                        B[k] = B[k] - V2*val
                    else:
                        pass
            ordre_tri = ordre_tri+1
    A = csr_array((valeur, (ligne, colonne)),
                  shape=(len(inconnus), len(inconnus)))
    return A, B


def init_maillage_coaxial_matlab():
    print('Récupération du maillage')
    points = pd.read_csv('nodes.csv').values
    triangles = pd.read_csv('elements.csv').values-1
    print('Récupération des labels')
    labels = pd.read_csv('labels.csv').values.flatten()
    return triangles, points, labels


def coaxial_matlab():
    params = pd.read_csv('params.csv').values.flatten()
    r1 = params[0]
    r2 = params[1]
    V1 = params[2]
    V2 = params[3]
    a = params[4]
    f = params[5]
    print('Paramètres')
    print(
        f'r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')

    print('Récupération de la solution')
    sol_matl = pd.read_csv('solution.csv').values.flatten()

    print('Construction du maillage')
    triangles, points, labels = init_maillage_coaxial_matlab()
    print('Construction des listes')
    liste_jacobiens = jacobiens(triangles, points)
    liste_gradients = gradients(triangles, points)

    print('Génération des matrices de masse et de rigidité')
    A, B = sys_lin(triangles, labels, a, f, V1, V2,
                   liste_jacobiens, liste_gradients)

    print('Résolution du système linéaire')
    sol_0 = spsolve(A, B)

    print('Tracé des solutions pour comparaison')
    plt.figure(figsize=(8, 6))
    plt.triplot(points[:, 0], points[:, 1],
                triangles, lw=0.5, ms=100)
    plt.plot(points[labels == 0, 0], points[labels == 0, 1], 'bo', ms=1)
    plt.plot(points[labels == 1, 0], points[labels == 1, 1], 'ro', ms=2)
    plt.plot(points[labels == 2, 0], points[labels == 2, 1], 'go', ms=2)
    plt.title('Maillage Matlab')
    plt.axis('equal')

    inconnus = np.where(labels == 0)[0]
    cond_dir_1 = np.where(labels == 1)[0]
    cond_dir_2 = np.where(labels == 2)[0]
    sol_1 = V1*np.ones(len(cond_dir_1))
    sol_2 = V2*np.ones(len(cond_dir_2))
    sol = np.concatenate((sol_1, sol_2, sol_0))
    trace_points = np.concatenate((cond_dir_1, cond_dir_2, inconnus))

    ax = plt.figure(figsize=(8, 6)).add_subplot(projection='3d')
    ax.scatter(points[trace_points, 0], points[trace_points, 1], sol, s=1)
    ax.set_title(
        f'Solution Python  r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')

    cond = len(sol_1)+len(sol_2)
    ax = plt.figure(figsize=(8, 6)).add_subplot(projection='3d')
    ax.scatter(points[trace_points[cond:-1], 0],
               points[trace_points[cond:-1], 1], sol[cond:-1]-sol_matl[cond:-1], s=1)
    ax.set_title(
        f'Erreur max = {np.max(np.abs((sol[cond:-1]-sol_matl[cond:-1]))):.2e} (Python vs. Matlab)')
    plt.show()


def init_maillage_coaxial_circ(r1=0.2, r2=1.0):
    n_angles = 50
    n_rayons = 30
    r = np.linspace(r1, r2, n_rayons)

    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    angles = np.repeat(angles[..., np.newaxis], n_rayons, axis=1)
    angles[:, 1::2] += np.pi / n_angles

    x = (r*np.cos(angles)).flatten()
    y = (r*np.sin(angles)).flatten()

    points = np.array([[0, 0]])
    labels = np.array([-1], dtype=np.int64)
    for k in range(x.shape[0]):
        points = np.append(points, [[x[k], y[k]]], axis=0)
        if np.abs(x[k]**2+y[k]**2 - (r1**2)) < 1e-10:
            labels = np.append(labels, 1)
        elif np.abs(x[k]**2+y[k]**2 - (r2**2)) < 1e-10:
            labels = np.append(labels, 2)
        else:
            labels = np.append(labels, 0)

    # Perform Delaunay triangulation
    tri = Delaunay(points)
    triangles = tri.simplices
    noeud = 0
    liste_adj_triangles = adj_triangles(triangles, points)
    index_tri = liste_adj_triangles[noeud]
    triangles = np.delete(triangles, index_tri, axis=0)

    return triangles, points, labels


def coaxial_circ():
    r1 = 0.2
    r2 = 1.0
    V1 = 12
    V2 = 10
    a = 10.0
    f = -10.0
    print('Paramètres')
    print(
        f'r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')

    print('Construction du maillage')
    triangles, points, labels = init_maillage_coaxial_circ()
    print('Construction des listes')
    liste_jacobiens = jacobiens(triangles, points)
    liste_gradients = gradients(triangles, points)

    print('Génération des matrices de masse et de rigidité')
    A, B = sys_lin(triangles, labels, a, f, V1, V2,
                   liste_jacobiens, liste_gradients)

    print('Résolution du système linéaire')
    sol_0 = spsolve(A, B)

    print('Tracé de la solution')
    plt.figure(figsize=(8, 6))
    plt.triplot(points[:, 0], points[:, 1],
                triangles, lw=0.5, ms=100)
    plt.plot(points[labels == 0, 0], points[labels == 0, 1], 'bo', ms=1)
    plt.plot(points[labels == 1, 0], points[labels == 1, 1], 'ro', ms=2)
    plt.plot(points[labels == 2, 0], points[labels == 2, 1], 'go', ms=2)
    plt.title('Maillage Python')
    plt.axis('equal')

    inconnus = np.where(labels == 0)[0]
    cond_dir_1 = np.where(labels == 1)[0]
    cond_dir_2 = np.where(labels == 2)[0]
    sol_1 = V1*np.ones(len(cond_dir_1))
    sol_2 = V2*np.ones(len(cond_dir_2))
    sol = np.concatenate((sol_0, sol_1, sol_2))
    trace_points = np.concatenate((inconnus, cond_dir_1, cond_dir_2))

    ax = plt.figure(figsize=(8, 6)).add_subplot(projection='3d')
    ax.scatter(points[trace_points, 0], points[trace_points, 1], sol, s=1)
    ax.set_title(
        f'Solution Python  r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')

    plt.show()


def init_maillage_coaxial_hexa(r1=0.2, r2=1.0):
    h = 0.05
    xmin = -1
    xmax = 1+h
    ymin = -25*h*np.sqrt(3)/2
    ymax = 26*h*np.sqrt(3)/2
    xl = np.arange(xmin, xmax, h)
    yl = np.arange(ymin, ymax, h*np.sqrt(3)/2)
    a1cer = np.arange(0, 2*np.pi, 2*np.pi/32)
    a2cer = np.arange(0, 2*np.pi, 2*np.pi/100)

    points = np.array([[0, 0]])
    labels = np.array([-1], dtype=np.int64)
    for a1 in a1cer:
        points = np.append(points, [[r1*np.cos(a1), r1*np.sin(a1)]], axis=0)
        labels = np.append(labels, 1)
    for a2 in a2cer:
        points = np.append(points, [[r2*np.cos(a2), r2*np.sin(a2)]], axis=0)
        labels = np.append(labels, 2)

    ind = 0
    for y in yl:
        for x in xl:
            if ind % 2 == 0:
                rexp2 = x**2+y**2
                if rexp2 > r1**2 and rexp2 < r2**2:
                    points = np.append(points, [[x, y]], axis=0)
                    labels = np.append(labels, 0)
            else:
                rexp2 = (x+h/2)**2+y**2
                if rexp2 > r1**2 and rexp2 < r2**2:
                    points = np.append(points, [[x+h/2, y]], axis=0)
                    labels = np.append(labels, 0)
        ind = ind+1

    for noeud in range(points.shape[0]):
        if labels[noeud] == 0:
            liste_noeuds = np.arange(points.shape[0])
            liste_noeuds = np.delete(liste_noeuds, noeud)
            if any(np.sqrt((points[liste_noeuds, 0]-points[noeud, 0])**2+(points[liste_noeuds, 1]-points[noeud, 1])**2) < h/3):
                labels[noeud] = -2
    points = np.delete(points, np.where(labels == -2)[0], axis=0)
    labels = np.delete(labels, np.where(labels == -2)[0])

    # Perform Delaunay triangulation
    tri = Delaunay(points)
    triangles = tri.simplices
    noeud = 0
    liste_adj_triangles = adj_triangles(triangles, points)
    index_tri = liste_adj_triangles[noeud]
    triangles = np.delete(triangles, index_tri, axis=0)

    return triangles, points, labels


def coaxial_hexa():
    r1 = 0.2
    r2 = 1.0
    V1 = 12
    V2 = 10
    a = 10.0
    f = -10.0
    print('Paramètres')
    print(
        f'r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')
    print('Construction du maillage')
    triangles, points, labels = init_maillage_coaxial_hexa()
    print('Construction des listes')
    liste_jacobiens = jacobiens(triangles, points)
    liste_gradients = gradients(triangles, points)

    print('Génération des matrices de masse et de rigidité')
    A, B = sys_lin(triangles, labels, a, f, V1, V2,
                   liste_jacobiens, liste_gradients)

    print('Résolution du système linéaire')
    sol_0 = spsolve(A, B)

    print('Tracé de la solution')
    plt.figure(figsize=(8, 6))
    plt.triplot(points[:, 0], points[:, 1],
                triangles, lw=0.5, ms=100)
    plt.plot(points[labels == 0, 0], points[labels == 0, 1], 'bo', ms=1)
    plt.plot(points[labels == 1, 0], points[labels == 1, 1], 'ro', ms=2)
    plt.plot(points[labels == 2, 0], points[labels == 2, 1], 'go', ms=2)
    plt.title('Maillage Python')
    plt.axis('equal')

    inconnus = np.where(labels == 0)[0]
    cond_dir_1 = np.where(labels == 1)[0]
    cond_dir_2 = np.where(labels == 2)[0]
    sol_1 = V1*np.ones(len(cond_dir_1))
    sol_2 = V2*np.ones(len(cond_dir_2))
    sol = np.concatenate((sol_0, sol_1, sol_2))
    trace_points = np.concatenate((inconnus, cond_dir_1, cond_dir_2))

    ax = plt.figure(figsize=(8, 6)).add_subplot(projection='3d')
    ax.scatter(points[trace_points, 0], points[trace_points, 1], sol, s=1)
    ax.set_title(
        f'Solution Python  r1={r1:.2f}, r2={r2:.2f}, V1={V1:.2f}, V2={V2:.2f}, a={a:.2f}, f={f:.2f}')

    plt.show()
