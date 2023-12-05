import numpy as np
from numba import njit

def exclude_unconnected_samples(distance_matrix): #ensures that every sample in the matrix is connected to every other sample at arbitrary degree
    if distance_matrix.shape[0] < 2:
        return
    unconnected_sample_idxs = get_unconnected_sample_idxs(distance_matrix)
    distance_matrix[unconnected_sample_idxs, :] = np.inf
    distance_matrix[:, unconnected_sample_idxs] = np.inf

def get_unconnected_sample_idxs(lower_matrix):
    full_matrix = convert_lower_to_full_matrix(lower_matrix)
    sums = np.sum(np.isfinite(full_matrix), axis=1)
    starting_sample = np.argmax(sums) #one of the traces with the most neighbors
    num_samples = full_matrix.shape[0]

    connected_samples = np.zeros(num_samples, dtype=np.bool_)
    check_connected_traces(full_matrix, starting_sample, connected_samples) #fills the set
    unconnected_samples = np.where(~connected_samples)[0]
    return unconnected_samples


def convert_lower_to_full_matrix(lower_matrix):
    full_matrix = np.copy(lower_matrix)
    rows, cols = np.where((lower_matrix != np.inf))
    full_matrix[cols, rows] = lower_matrix[rows, cols]
    return full_matrix

@njit
def check_connected_traces(matrix, trace_idx, visited):
    neighbors = np.where(matrix[trace_idx] != np.inf)[0]
    for neighbor in neighbors:
        if not visited[neighbor]:
            visited[neighbor] = True
            check_connected_traces(matrix, neighbor, visited)


