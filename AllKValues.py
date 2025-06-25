import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg
np.set_printoptions(precision=250)

# A user can input their own matrix or they can use a predefined matrix
def inputMatrix():
    # while True:
    #     n = int(input("Enter the n size of your matrix: "))
    #     if n > 1:
    #         break
    #     else:
    #         print("\nInvalid input, please try again...")
    
    # matrixInput = input(f"Please put your {n * n} elements here separated by space: ")
    # elements = [eval(x) for x in matrixInput.split()]
    # A = np.matrix(elements).reshape(n, n)

    omega = (-1+1j*np.sqrt(3))/2
    A = np.matrix([
        [1, 0, 0, 0, 0],
        [0, omega, 0, 0, 0],
        [0, 0, omega**2, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0]])
    
    # A = np.matrix([
    #     [1, 0, 0, 0, 0],
    #     [0, 2, 0, 0, 0],
    #     [0, 0, 3j, 0, 0],
    #     [0, 0, 0, 4j, 0],
    #     [0, 0, 0, 0, 5j]])

    # ep1 = 1e-25 + 1 
    # ep2 = 1e-16 + 1
    # ep3 = 1e-20 + 1
    # A = np.matrix([
    #         [2, 0, 0, 0, 0],
    #         [0, ep1 + 1j, 0, 0, 0],
    #         [0, 0, ep2, 0, 0],
    #         [0, 0, 0, ep3 - 1j, 0],
    #         [0, 0, 0, 0, -1]])

    # A = np.matrix([
    #         [1j, 0, 0, 0, 0],
    #         [0, -1j, 0, 0, 0],
    #         [0, 0, -1, 0, 0],
    #         [0, 0, 0, -2, 1],
    #         [0, 0, 0, 0, 0]])

    # A = np.matrix([
    #         [1, 0, 0, 0],
    #         [0, 1j, 0, 0],
    #         [0, 0, -1, 0],
    #         [0, 0, 0, -1j]])
    
    # A = np.asmatrix(np.diag(np.array([1j,-1j,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])))

    # A = np.matrix([
    #         [0, 0, 0, 0, 0, 0, 0],
    #         [0, 1+1j, 0, 0, 0, 0, 0],
    #         [0, 0, 1+1j, 0, 0, 0, 0],
    #         [0, 0, 0, 1-1j, 0, 0, 0],
    #         [0, 0, 0, 0, 2, 0, 0],
    #         [0, 0, 0, 0, 0, 2, 0],
    #         [0, 0, 0, 0, 0, 0, 3+2j]])

    # A = np.asmatrix(np.diag(np.array([(1 + 1/(m+1))*np.exp(2 * np.pi * 1j * m/4) for m in np.arange(8)])))
    # A = np.matrix([
    #         [3j, 0, 0, 0, 0],
    #         [0, 1j, 0, 0, 0],
    #         [0, 0, -1j, 0, 0],
    #         [0, 0, 0, -3j, 0],
    #         [0, 0, 0, 0, -1]])

    # A = np.matrix([[1, (5-9j)/2], 
    #           [(5+3j)/2, 4]])

    # A = np.matrix([
    #               [0, 1], 
    #               [0, 0]])

    # A = np.matrix([[1, (5-9j)/2, 2+1j], 
    #           [(5+3j)/2, 4, -1j], 
    #           [3-2j, 1j, 2]])
    # print("Matrix A:\n", A)
    n = A.shape[0]
    return A, n

# Returns the top-k eigenvectors of a Hermitian matrix A
def vecs(A, k):
    eigval, eigvec = linalg.eigh(A)
    return eigvec[:, -k:]

# Computes the weighted sum x_1*Ax_1 + ... + x_k*Ax_k divided by k
def nrsum(A, k, V):
    return np.einsum('ji, jk, ki ->', np.conj(V), A, V) / k 

# Computes the boundary points of the k-numerical range for a given matrix A
def compute_numerical_range(A, k, n, step_size=0.01):
    boundary_points = []
    
    for theta in np.arange(0, 2 * np.pi, step_size):
        A_theta = np.exp(-1j * theta) * A
        H_theta = (A_theta + A_theta.H) / 2
        
        eigval, eigvec = linalg.eigh(H_theta)
        
        if eigval[n-k] != eigval[n-k-1]:
            V = vecs(H_theta, k)
            boundary_points.append(nrsum(A, k, V))
        else:
            alpha = eigval[-k]
            is_close = np.isclose(eigval, alpha)
            vals1 = eigval[(eigval > alpha) & (~is_close)]
            vecs1 = eigvec[:, (eigval > alpha) & (~is_close)] 
            vecs2 = eigvec[:, is_close]

            m = vals1.size
            k_prime = k - m
            
            B = vecs2
            K_theta = (A_theta - A_theta.H) / (2j)
            
            W = B * vecs(B.H * K_theta * B, k_prime)
            U = B * vecs(-B.H * K_theta * B, k_prime)
            
            boundary_points.append(nrsum(A, k, np.hstack([U, vecs1])))
            boundary_points.append(nrsum(A, k, np.hstack([W, vecs1]))) 
            
    boundary_points.append(boundary_points[0])
    return np.array(boundary_points)

# Plots the k-numerical range for all k = 1 to n
def plot_all_k(A, n):
    plt.figure(figsize=(8, 8))

    for k in range(1, n + 1):
        boundary_points = compute_numerical_range(A, k, n)
        plt.plot(boundary_points.real, boundary_points.imag, marker='o',
                 linestyle='-', markersize=4, label=f'k = {k}')
    
    plt.scatter([0], [0], color='red', marker='x', label='Origin (0,0)')
    plt.title("Numerical Range Boundary Visualization for All k")
    plt.xlabel("Real Part")
    plt.ylabel("Imaginary Part")
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.show()

# Runs the script
A, n = inputMatrix()
plot_all_k(A, n)
