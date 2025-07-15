from tools.create_system import Create_System
from sympy import *


def compute_rank(matrix):
    """Compute the column rank of a symbolic matrix."""
    try:
        rank = matrix.rank()
        print(f"Matrix dimensions: {matrix.rows} x {matrix.cols}")
        print(f"Computed column rank: {rank}")
        return rank
    except Exception as e:
        print(f"Error computing rank: {e}")
        return None


def print_matrix(matrix, title="Matrix"):
    """Print a matrix with aligned rows."""
    print(f"\n{title}:")
    for i in range(matrix.rows):
        row_str = "["
        for j in range(matrix.cols):
            row_str += str(matrix[i, j])
            if j < matrix.cols - 1:
                row_str += ", "
        row_str += "]"
        print(row_str)


def check_rank_condition(M, X, r, A, Ba):
    """
    Check rank conditions for matrices.

    Returns:
        0: if rank([M; X*A; X*Ba]) > rank([M; X*A])
        1: if rank([M; X*A; X*Ba]) == rank([M; X*A])
        -1: if rank([M; X*Ba]) > r
        None: if rank([M; X*Ba]) == r
    """
    XA = X * A
    M_XA = M.col_join(XA)
    rank_M_XA = M_XA.rank()

    print(f"Rank of [M; X*A]: {rank_M_XA}")

    if rank_M_XA > r:
        print(f"Rank of [M; X*A] ({rank_M_XA}) > r ({r})")

        XBa = X * Ba
        M_XA_XBa = M_XA.col_join(XBa)
        rank_M_XA_XBa = M_XA_XBa.rank()

        print(f"Rank of [M; X*A; X*Ba]: {rank_M_XA_XBa}")

        if rank_M_XA_XBa > rank_M_XA:
            print(f"Rank of [M; X*A; X*Ba] ({rank_M_XA_XBa}) > rank of [M; X*A] ({rank_M_XA})")
            return 0
        else:
            print(f"Rank of [M; X*A; X*Ba] ({rank_M_XA_XBa}) == rank of [M; X*A] ({rank_M_XA})")
            return 1
    else:
        print(f"Rank of [M; X*A] ({rank_M_XA}) == r ({r})")

        XBa = X * Ba
        M_XBa = M.col_join(XBa)
        rank_M_XBa = M_XBa.rank()

        print(f"Rank of [M; X*Ba]: {rank_M_XBa}")

        if rank_M_XBa > r:
            print(f"Rank of [M; X*Ba] ({rank_M_XBa}) > r ({r})")
            return -1
        else:
            print(f"Rank of [M; X*Ba] ({rank_M_XBa}) == r ({r})")
            return None


def get_preset_matrices(preset_num):
    """Get specific preset matrices based on selection."""
    # Define symbolic parameters
    a, b, c, d, k = symbols('a b c d k', real=True, nonzero=True)

    presets = {
        1: {
            'A': Matrix([
                [0, a, 0],
                [a, 0, 0],
                [0, 0, b]
            ]),
            'B': Matrix([
                [1, 0, -1],
                [0, 0, 0],
                [1, 0, 0]
            ]),
            'size': 3,
            'description': "3x3 system with cancellation",
            'parameters': 'a, b'
        },

        2: {
            'A': Matrix([
                [0, -1, 0, 0],
                [-1, 0, 0, 0],
                [0, 0, 0, a],
                [0, 0, a, 0]
            ]),
            'B': Matrix([
                [0, 0, 0, -1],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [1, 0, 0, b]
            ]),
            'size': 4,
            'description': "Timoshenko",
            'parameters': 'a, b'
        },

        3: {
            'A': Matrix([
                [0, 0, -1, 0, 0],
                [0, 0, 0, -c, d],
                [-1, 0, 0, 0, 0],
                [0, -c, 0, 0, 0],
                [0, d, 0, 0, 0]
            ]),
            'B': Matrix([
                [0, -1, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1]
            ]),
            'size': 5,
            'description': "Timoshenko with Memory",
            'parameters': 'c, d'
        },

        4: {
            'A': Matrix([
                [0, 0, -1, 0, 0, 0],
                [0, 0, 0, -c, d, 0],
                [-1, 0, 0, 0, 0, 0],
                [0, -c, 0, 0, 0, 0],
                [0, d, 0, 0, 0, k],
                [0, 0, 0, 0, k, 0]
            ]),
            'B': Matrix([
                [0, -1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1]
            ]),
            'size': 6,
            'description': "Timoshenko-Cattaneo",
            'parameters': 'c, d, k'
        }
    }

    return presets.get(preset_num, None)


def get_matrices():
    """Get matrices from user input or presets."""
    print("Choose an option:")
    print("1. Input matrices manually")
    print("2. Generate random matrices (HIGHLY UNSTABLE FOR NOW! Not recommended)")
    print("3. Use preset matrices")

    choice = input("Enter your choice (1/2/3): ").strip()

    if choice == '1':
        size = int(input("Enter the size of the matrices: "))
        splitter = Create_System(size)
        A, Ba, Bs = splitter.split_matrices()
        return A, Ba, Bs, size

    elif choice == '2':
        size = int(input("Enter the size of the matrices: "))
        splitter = Create_System(size)
        A, Ba, Bs = splitter.generate_random_matrices()
        return A, Ba, Bs, size

    elif choice == '3':
        print("\nAvailable presets:")
        print("1. 3x3 system with cancellation")
        print("2. Timoshenko")
        print("3. Timoshenko with Memory")
        print("4. Timoshenko-Cattaneo")

        preset_choice = input("Enter preset number (1-4): ").strip()

        try:
            preset_num = int(preset_choice)
            preset_data = get_preset_matrices(preset_num)

            if preset_data is None:
                print("Invalid preset number. Exiting.")
                return None, None, None, None

            A_preset = preset_data['A']
            B_preset = preset_data['B']
            size = preset_data['size']
            description = preset_data['description']
            parameters = preset_data['parameters']

            # Calculate Ba and Bs
            Bs = (B_preset + B_preset.T) / 2
            Ba = (B_preset - B_preset.T) / 2
            A = A_preset

            print(f"\nUsing preset {preset_num}: {description}")
            print(f"Size: {size}x{size}")
            print(f"Parameters: {parameters} (all nonzero)")

            return A, Ba, Bs, size

        except ValueError:
            print("Invalid input. Please enter a number between 1 and 4.")
            return None, None, None, None

    else:
        print("Invalid choice. Exiting.")
        return None, None, None, None