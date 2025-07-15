import sys
from tools.tree import *
from tools.latex import *
from tools.matrix import *


def main():
    """Main function to run binary tree exploration and generate LaTeX output."""
    output_file_name = "output.txt"

    # Get matrices - user input happens here before redirection
    A, Ba, Bs, size = get_matrices()

    if A is None:
        print("Matrix acquisition failed or was cancelled. Exiting.")
        return

    # Save reference to console output
    original_stdout = sys.stdout

    try:
        # Redirect output to file
        with open(output_file_name, 'w', encoding='utf-8') as f:
            sys.stdout = f

            # Display matrices
            print_matrix(A, "Matrix A")
            print_matrix(Bs, "Matrix Bs (Symmetric Part)")
            print_matrix(Ba, "Matrix Ba (Antisymmetric Part)")

            # Explore the binary tree
            root, final_rank = explore_tree(A, Ba, Bs, size)

            # Print results
            print(f"\n{'=' * 60}")
            print(f"FINAL TREE STRUCTURE")
            print(f"{'=' * 60}")
            print_custom_tree(root)

            print(f"\nFinal rank achieved: {final_rank}")
            print(f"Target rank: {size}")
            print(f"Exploration {'completed successfully' if final_rank >= size else 'incomplete'}")

            # Build Lyapunov functional
            print(f"\n{'=' * 60}")
            print(f"BUILDING LYAPUNOV FUNCTIONAL")
            print(f"{'=' * 60}")

            u_symbols = symbols(' '.join([f'u_{i + 1}' for i in range(size)]))
            U = Matrix(u_symbols)
            m = symbols('m', real=True)

            functional = build_lyapunov(root, U, A, Ba, size, m)

            # Output LaTeX
            print(f"\n{'=' * 60}")
            print(f"LATEX OUTPUT")
            print(f"{'=' * 60}")

            latex_output = functional_to_latex(A, Ba, U, root)
            print(latex_output)

    finally:
        # Restore console output
        sys.stdout = original_stdout

    print(f"\nAll program output has been saved to '{output_file_name}'")


if __name__ == "__main__":
    main()