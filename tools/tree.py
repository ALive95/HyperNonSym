from sympy import *
from collections import deque
from tools.latex import *


class TreeNode:
    """Node in the binary tree exploration."""

    def __init__(self, matrix, name, level=0, parent=None):
        self.matrix = matrix
        self.name = name
        self.level = level
        self.parent = parent
        self.children = []
        self.processed = False
        self.direction = None
        self.discrepancy = self._calc_discrepancy()
        self.number = self._calc_number()

    def _calc_discrepancy(self):
        """Calculate discrepancy: 0 for root or nodes ending in A, 1 for nodes ending in Ba."""
        if self.name == "Bs":
            return 0
        elif self.name.endswith("A"):
            return 0
        elif self.name.endswith("Ba"):
            return 1

    def _calc_number(self):
        """Calculate number based on parent's discrepancy."""
        if self.name == "Bs":
            return 0
        elif self.parent is not None:
            parent_number = self.parent.number
            parent_discrepancy = self.parent.discrepancy

            if parent_discrepancy == 0:
                return parent_number
            elif parent_discrepancy == 1:
                return parent_number + 1
            else:
                return parent_number
        else:
            return 0

    def add_child(self, child):
        """Add a child node."""
        self.children.append(child)
        child.parent = self
        child.level = self.level + 1


def explore_tree(A, Ba, Bs, size, max_iterations=10):
    """Explore the binary tree based on rank conditions."""
    root = TreeNode(Bs, "Bs", level=0)
    M = Bs.copy()
    current_rank = compute_rank(M)

    print(f"\n{'=' * 60}")
    print(f"BINARY TREE EXPLORATION")
    print(f"{'=' * 60}")
    print(f"Initial rank: {current_rank}")
    print(f"Target rank: {size}")

    queue = deque([(root, M, [root])])
    iteration = 0

    while queue and iteration < max_iterations:
        iteration += 1
        print(f"\n{'=' * 40}")
        print(f"ITERATION {iteration}")
        print(f"{'=' * 40}")

        current_node, M, current_leaves = queue.popleft()

        if current_rank >= size:
            print(f"Target rank {size} reached! Stopping exploration.")
            break

        new_leaves = []
        new_matrices = []

        for leaf in current_leaves:
            if leaf.processed:
                continue

            print(f"\nProcessing leaf: {leaf.name}")
            print(f"Current M rank: {current_rank}")

            result = check_rank_condition(M, leaf.matrix, current_rank, A, Ba)
            leaf.direction = result

            print(f"Rank condition result for {leaf.name}: {result}")

            # Add children based on result
            if result == 1:
                child = TreeNode(leaf.matrix * A, f"{leaf.name} A", parent=leaf)
                leaf.add_child(child)
                new_leaves.append(child)
                new_matrices.append(child.matrix)
                print(f"Added child: {child.name}")

            elif result == -1:
                child = TreeNode(leaf.matrix * Ba, f"{leaf.name} Ba", parent=leaf)
                leaf.add_child(child)
                new_leaves.append(child)
                new_matrices.append(child.matrix)
                print(f"Added child: {child.name}")

            elif result == 0:
                # Add both children
                child_A = TreeNode(leaf.matrix * A, f"{leaf.name} A", parent=leaf)
                child_Ba = TreeNode(leaf.matrix * Ba, f"{leaf.name} Ba", parent=leaf)
                leaf.add_child(child_A)
                leaf.add_child(child_Ba)
                new_leaves.extend([child_A, child_Ba])
                new_matrices.extend([child_A.matrix, child_Ba.matrix])
                print(f"Added children: {child_A.name}, {child_Ba.name}")

            else:
                print(f"No children added for {leaf.name}")

            leaf.processed = True

        # Update M with new matrices
        if new_matrices:
            new_M = M.copy()
            for matrix in new_matrices:
                new_M = new_M.col_join(matrix)

            current_rank = compute_rank(new_M)
            print(f"\nUpdated M with new leaves. New rank: {current_rank}")

            if current_rank < size and new_leaves:
                queue.append((current_node, new_M, new_leaves))

        if current_rank >= size:
            print(f"\nTarget rank {size} reached!")
            break

    if iteration >= max_iterations:
        print(f"\nMaximum iterations ({max_iterations}) reached. Stopping exploration.")

    return root, current_rank


def build_lyapunov(root, U, A, Ba, size, m=1, dx=None, xi=None):
    """Build the Lyapunov functional based on tree exploration results."""
    from sympy import symbols, Rational

    if dx is None:
        dx = symbols('dx')
    if xi is None:
        xi = symbols('xi')

    # Start with base term: (1/2)||U||^2
    functional = Rational(1, 2) * U.T * U
    functional = functional[0]

    print(f"\nBuilding Lyapunov functional...")
    print(f"Initial term: (1/2)||U||^2")

    def add_terms(node, level=1):
        nonlocal functional

        if node.parent and node.parent.direction is not None:
            indent = "  " * level
            print(
                f"\n{indent}Processing node {node.name}, node number {node.number}, Direction of parent {node.parent.direction}")

            node_U = node.parent.matrix * U

            if node.parent.direction == 1:
                # Add 1/ξ^(2*(1+node.number)) * ⟨parent*U, node*∂_x*U⟩
                node_dx_U = node.matrix * dx * U
                term = (node_dx_U.T * node_U)[0] / (xi ** (2 * (1 + node.number)))
                functional += term
                print(f"{indent}Added term: (1/ξ^{2 * (1 + node.number)}) ⟨{node.parent.name} U, {node.name} ∂_x U⟩")
                print(f"{indent}Term: {term}")

            elif node.parent.direction == -1:
                # Add 1/ξ^(2*(1+node.number)) * ⟨parent*U, node*U⟩
                node_U_term = node.matrix * U
                term = (node_U_term.T * node_U)[0] / (xi ** (2 * (1 + node.number)))
                functional += term
                print(f"{indent}Added term: (1/ξ^{2 * (1 + node.number)}) ⟨{node.parent.name} U, {node.name} U⟩")
                print(f"{indent}Term: {term}")

            elif node.parent.direction == 0:
                if node.name.endswith('A'):
                    # Mixed terms for A direction
                    node_dx_U = node.matrix * dx * U
                    node_A_U = node.parent.matrix * A * U
                    node_Ba_A_U = node.parent.matrix * Ba * A * U

                    term1 = (node_dx_U.T * node_U)[0]
                    term2 = 2 * m * (node_Ba_A_U.T * node_A_U)[0]

                    total_term = (term1 + term2) / (xi ** (2 * (1 + node.number)))
                    functional += total_term
                    print(
                        f"{indent}Added mixed A terms: (1/ξ^{2 * (1 + node.number)}) [⟨{node.parent.name} U, {node.name} ∂_x U⟩ + m⟨{node.parent.name} A U, {node.parent.name} Ba A U⟩]")
                    print(f"{indent}Term: {total_term}")

                else:  # Ends with Ba
                    # Mixed terms for Ba direction
                    node_Ba_U = node.matrix * U
                    node_A_U = node.parent.matrix * A * U
                    node_Ba_A_U = node.parent.matrix * Ba * A * U

                    term1 = (node_Ba_U.T * node_U)[0]
                    term2 = 2 * m * (node_Ba_A_U.T * node_A_U)[0]

                    total_term = (term1 + term2) / (xi ** (2 * (1 + node.number)))
                    functional += total_term
                    print(
                        f"{indent}Added mixed Ba terms: (1/ξ^{2 * (1 + node.number)}) [⟨{node.parent.name} U, {node.name} U⟩ + m⟨{node.parent.name} A U, {node.parent.name} Ba A U⟩]")
                    print(f"{indent}Term: {total_term}")

        # Process children
        for child in node.children:
            add_terms(child, level + 1)

    # Start from root's children
    for child in root.children:
        add_terms(child, level=0)

    print(f"\nFinal Lyapunov functional (use this to check the LaTeX one):\n"
          f" {functional}")

    return functional


def print_custom_tree(node, indent="", is_last=True):
    """Print the tree structure."""
    if node is None:
        return

    print(indent + ("└── " if is_last else "├── ") + node.name)

    for i, child in enumerate(node.children):
        is_last_child = (i == len(node.children) - 1)
        print_custom_tree(child, indent + ("    " if is_last else "│   "), is_last_child)

