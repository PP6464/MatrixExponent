using LinearAlgebra

function matrix_exp(M::Matrix{T}) where T <: Number
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
        is_diagonal = false
        break
      end
    end
  end
  # If the matrix is diagonal simply return elementwise exponentiated M (only along lead diagonal)
  if is_diagonal
    res_mat = zeros(Number, size(M)[1], size(M)[2])
    for i ∈ 1:size(M)[1]
      res_mat[i, i] = ℯ^M[i, i]
    end
    return res_mat
  end
  # Check if matrix has sufficient eigenvalues
  eigen_info = eigen(M)
  if length(unique(eigen_info.values)) == size(M)[1]
    # There are sufficient eigenvalues
    eigenvalues = zeros(Number, size(M)[1], size(M)[2])
    for (i, v) ∈ enumerate(eigen_info.values)
      eigenvalues[i, i] = v
    end
    return eigen_info.vectors * matrix_exp(eigenvalues) * inv(eigen_info.vectors)
  end
  # Check if matrix is nilpotent
  k = 0
  M_powers = [Matrix{Int64}(I, size(M)[1], size(M)[2])]
  while k < size(M)[1]
    push!(M_powers, M_powers[end] * M)
    k += 1
    if M_powers[end] == zeros(ComplexF64, size(M)[1], size(M)[2])
      println("nilpotent")
      # M is nilpotent with degree k
      res_mat = zeros(ComplexF64, size(M)[1], size(M)[2])
      for i ∈ 1:k+1
        res_mat += M_powers[i] / factorial(i - 1)
      end
    end
  end
  # Return default matrix exponent function (uses approximation)
  return exp(M)
end
