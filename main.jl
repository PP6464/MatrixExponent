using LinearAlgebra

function matrix_exp(M::Matrix{Number})
  # Throw error if matrix is not square
  if size(M)[1] != size(M)[2]
    throw(ArgumentError("The matrix must be square"))
  end
  # First check if the matrix is diagonal
  is_diagonal = true
  for i ∈ 1:size(M)[1]
    if !is_diagonal
      break
    end
    for j ∈ 1:size(M)[2]
      if i != j && M[i, j] != 0
        is_diagonal = true
        break
      end
    end
  end
  # If the matrix is diagonal simply return elementwise exponentiated M (only along lead diagonal)
  if is_diagonal
    res_mat = zeros(Number, size(M)[1], size(M)[2])
    for i ∈ 1:size(M)[1]
      res_mat[i, i] = exp(M[i, i])
    end
    return res_mat
  end
  # Check if matrix has sufficient eigenvalues
  eigen_info = eigen(M)
  if length(unique(eigen_info.values)) == size(M)[1]
    # There are sufficient eigenvalues
    eigenvalues = zeros(Number, size(M)[1], size(M)[1])
    for i ∈ 1:eigen_info.values
      eigenvalues[i, i] = eigen_info.values[i]
    end
    return eigen_info.vectors * matrix_exp(eigenvalues) * inv(eigen_info.vectors)
  end
  # Check if matrix is nilpotent
  k = 1
  current_mat = M
  while true
    if M * current_mat == zeros(Number, size(M)[1], size(M)[2])
      # The matrix is nilpotent with nilpotency degree k
      res = zero(Number)
      for i ∈ 1:k - 1
        res += M^k / factorial(i)
      end
    else if k <= size(M)[1]
      k += 1
      current_mat = M * current_mat
    else
      break
    end
  end
  # Return default matrix exponent function (uses approximation)
  return exp(M)
end
