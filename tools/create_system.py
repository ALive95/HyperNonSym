from sympy import Matrix, symbols, zeros, sympify
import numpy as np
import random # For generating random symbols or values if desired

class Create_System:
    def __init__(self, size):
        self.size = size

    def is_symmetric(self, matrix):
        """Checks if a SymPy matrix is symmetric."""
        return matrix == matrix.T

    def splitsym(self, matrix):
        """Splits a SymPy matrix into its symmetric and antisymmetric parts."""
        symmetric = (matrix + matrix.T) / 2
        antisymmetric = (matrix - matrix.T) / 2
        return symmetric, antisymmetric

    def input_symbolic_matrix(self, matrix_name="matrix"):
        """
        Prompts the user to input elements for a symbolic matrix.
        Each element can be a number or a symbolic expression.
        """
        matrix_elements = []
        print(f"Enter {self.size}x{self.size} {matrix_name} (elements can be numbers or symbols, separated by spaces):")
        for i in range(self.size):
            while True:
                row_str = input(f"Row {i+1}: ").strip()
                try:
                    # Evaluate each element as a SymPy expression
                    row_elements = [sympify(expr) for expr in row_str.split()]
                    if len(row_elements) != self.size:
                        print(f"Invalid input! Please enter exactly {self.size} elements for the row.")
                    else:
                        matrix_elements.append(row_elements)
                        break
                except (ValueError) as e: # Catch ValueError specifically for input parsing
                    print(f"Invalid input! Please ensure elements are valid numbers or symbolic expressions. Error: {e}")
                except Exception as e: # Catch other potential sympify errors
                    print(f"Error processing input: {e}")
        return Matrix(matrix_elements)

    def split_matrices(self):
        """
        Gets matrix A and B from user input, ensuring A is symmetric,
        then splits B into its symmetric (Bs) and antisymmetric (Ba) parts.
        All matrices are SymPy symbolic matrices.
        """
        print("\n--- Input Matrix A ---")
        A = self.input_symbolic_matrix("Matrix A")
        if A is None:
            return None, None, None

        # Check if matrix A is symmetric
        if not self.is_symmetric(A):
            print("Warning: Matrix A is not symmetric. Proceeding anyway, but be aware.")
            # A = (A + A.T) / 2 # Uncomment to force A to be symmetric
            # return None, None, None # Uncomment to exit if A must be symmetric

        print("\n--- Input Matrix B ---")
        B = self.input_symbolic_matrix("Matrix B")
        if B is None:
            return None, None, None

        # Splitting matrix B into symmetric and antisymmetric parts
        B_symmetric, B_antisymmetric = self.splitsym(B)
        return A, B_antisymmetric, B_symmetric

    def generate_random_matrices(self):
        """
        Generates random symbolic matrices A (symmetric), Ba (antisymmetric),
        and Bs (symmetric) with a user-specified rank for Bs.
        """
        # Generate random integer elements for A using NumPy, then convert to SymPy Matrix
        A_np = np.random.randint(-2, 2, size=(self.size, self.size))
        A = Matrix((A_np + A_np.T)) # Ensure A is symmetric from random gen

        # --- Generate Bs with specific rank ---
        while True:
            try:
                desired_rank_Bs_str = input(f"Enter the desired rank for Bs (1 to {self.size}): ")
                desired_rank_Bs = int(desired_rank_Bs_str)
                if 1 <= desired_rank_Bs <= self.size:
                    break
                else:
                    print(f"Rank must be between 1 and {self.size}. Please try again.")
            except ValueError:
                print("Invalid input. Please enter an integer.")

        # Construct Bs with the desired rank
        U_elements = []
        for i in range(self.size):
            row = []
            for j in range(desired_rank_Bs):
                row.append(random.randint(-2, 2)) # Use random integers for simplicity
            U_elements.append(row)
        U = Matrix(U_elements)

        current_rank_U = U.rank()
        attempts = 0
        max_attempts = 100
        while current_rank_U < desired_rank_Bs and attempts < max_attempts:
            U_elements = []
            for i in range(self.size):
                row = []
                for j in range(desired_rank_Bs):
                    row.append(random.randint(-2, 2))
                U_elements.append(row)
            U = Matrix(U_elements)
            current_rank_U = U.rank()
            attempts += 1

        if current_rank_U < desired_rank_Bs:
            print(f"Warning: Could not reliably generate U with rank {desired_rank_Bs} after {max_attempts} attempts. Generated U has rank {current_rank_U}.")

        Bs = U * U.T

        B_temp_np = np.random.randint(-2, 2, size=(self.size, self.size))
        B_temp = Matrix(B_temp_np)
        Ba = (B_temp - B_temp.T)

        return A, Ba, Bs