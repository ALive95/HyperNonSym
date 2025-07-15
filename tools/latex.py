from sympy import *
from tools.matrix import *


def generate_l2_latex(vec1, vec2):
    """Generate L2 scalar product LaTeX notation."""
    from sympy import S, latex, Mul, Add

    if not (isinstance(vec1, Matrix) and isinstance(vec2, Matrix) and
            vec1.cols == 1 and vec2.cols == 1 and vec1.rows == vec2.rows):
        raise ValueError("Inputs must be SymPy column matrices of the same size.")

    l2_terms = []
    for i in range(vec1.rows):
        component1 = vec1[i, 0]
        component2 = vec2[i, 0]

        product = (component1 * component2).simplify()
        if product == 0:
            continue

        def extract_coeff_and_var(expr):
            """Extract coefficient and variable parts from expression."""
            if expr == 0:
                return S.Zero, S.One

            coeff_part = S.One
            var_part = S.One

            if isinstance(expr, Add):
                return expr, S.One

            for factor in Mul.make_args(expr):
                is_var = False

                if hasattr(factor, 'free_symbols') and factor.free_symbols:
                    # Check if it's a variable or derivative
                    factor_str = str(factor)
                    if any(sym.name.startswith('u_') for sym in factor.free_symbols):
                        is_var = True
                    elif factor_str.startswith('dx(') and ')' in factor_str:
                        is_var = True
                    elif str(factor) == 'dx':
                        is_var = True

                if is_var:
                    var_part *= factor
                else:
                    coeff_part *= factor

            if var_part == 1 and expr != 1:
                return expr, S.One

            return var_part.simplify(), coeff_part.simplify()

        var1, coeff1 = extract_coeff_and_var(component1)
        var2, coeff2 = extract_coeff_and_var(component2)

        overall_coeff = (coeff1 * coeff2).simplify()

        var1_latex = latex(var1).replace(r"dx", r"\partial_x")
        var2_latex = latex(var2).replace(r"dx", r"\partial_x")

        if overall_coeff == 1:
            l2_terms.append(f"\\langle {var1_latex}, {var2_latex} \\rangle")
        elif overall_coeff == -1:
            l2_terms.append(f"- \\langle {var1_latex}, {var2_latex} \\rangle")
        else:
            l2_terms.append(f"{latex(overall_coeff)}\\langle {var1_latex}, {var2_latex} \\rangle")

    if not l2_terms:
        return "0"

    # Join terms with proper signs
    result = ""
    for i, term in enumerate(l2_terms):
        if i > 0 and not term.strip().startswith('-'):
            result += " + "
        elif term.strip().startswith('-') and i > 0:
            result += " "
        result += term

    return result


def compute_m_symbolic(X, A, Ba, U):
    """
    Compute the symbolic value of m that satisfies:
    ⟨XAU, -XBaU + mXBaA²U⟩ = 0

    This gives us: m = ⟨XAU, XBaU⟩ / ⟨XAU, XBaA²U⟩

    Args:
        X: Matrix representing the node.parent.matrix
        A: Matrix A
        Ba: Matrix Ba (antisymmetric part of B)
        U: Vector U

    Returns:
        Symbolic expression for m, or None if denominator is zero
    """
    from sympy import simplify, Matrix

    # Compute the vectors for the scalar products
    XAU = X * A * U
    XBaU = X * Ba * U
    XBaA2U = X * Ba * A * A * U

    # Compute the numerator: ⟨XAU, XBaU⟩
    numerator = 0
    for i in range(XAU.rows):
        numerator += XAU[i, 0] * XBaU[i, 0]
    numerator = simplify(numerator)

    # Compute the denominator: ⟨XAU, XBaA²U⟩
    denominator = 0
    for i in range(XAU.rows):
        denominator += XAU[i, 0] * XBaA2U[i, 0]
    denominator = simplify(denominator)

    # Check if denominator is zero
    if denominator == 0:
        print("Warning: Denominator is zero, cannot compute m")
        return None

    # Compute m
    m = simplify(numerator / denominator)
    return m


def check_cancellations(X, A, Ba, U, m_value):
    """
    Check for parameter cancellations by solving:
    ⟨XU - mXA²U, XBaAU⟩ = 0

    Args:
        X: Matrix representing the node.parent.matrix
        A: Matrix A
        Ba: Matrix Ba (antisymmetric part of B)
        U: Vector U
        m_value: The previously computed value of m

    Returns:
        dict with:
        - 'solutions': List of parameter constraints that make the expression zero
        - 'expression': The simplified scalar product expression
        - 'status': 'solved', 'underdetermined', 'no_solution', or 'always_zero'
    """
    from sympy import solve, simplify, Eq, symbols

    # Compute the vectors for the scalar product
    XU = X * U
    XA2U = X * A * A * U
    XBaAU = X * Ba * A * U

    # Compute the first vector: XU - m*XA²U
    vec1 = XU - m_value * XA2U
    vec2 = XBaAU

    # Compute the scalar product ⟨XU - mXA²U, XBaAU⟩
    scalar_product = 0
    for i in range(vec1.rows):
        scalar_product += vec1[i, 0] * vec2[i, 0]

    scalar_product = simplify(scalar_product)

    # If already zero, no constraints needed
    if scalar_product == 0:
        return {
            'solutions': [],
            'expression': scalar_product,
            'status': 'always_zero'
        }

    # Get all free symbols (parameters) in the expression
    free_params = scalar_product.free_symbols

    # Remove symbols that represent derivatives or functions (like dx, u_0, etc.)
    actual_params = []
    for param in free_params:
        param_str = str(param)
        # Keep only simple parameter symbols (like a, b, c, d, k)
        if (len(param_str) == 1 and param_str.isalpha()) or param_str in ['xi']:
            actual_params.append(param)

    if not actual_params:
        # No parameters to solve for
        if scalar_product == 0:
            return {
                'solutions': [],
                'expression': scalar_product,
                'status': 'always_zero',
                'parameters': '--'
            }
        else:
            return {
                'solutions': [],
                'expression': scalar_product,
                'status': 'no_solution',
                'parameters': '--'
            }

    try:
        # Try to solve for parameter constraints
        solutions = solve(Eq(scalar_product, 0), actual_params)

        if not solutions:
            return {
                'solutions': [],
                'expression': scalar_product,
                'status': 'no_solution',
                'parameters': '--'
            }

        # Handle different types of solutions
        if isinstance(solutions, dict):
            # Single solution as dictionary
            solution_list = [solutions]
        elif isinstance(solutions, list):
            # Multiple solutions
            solution_list = solutions
        else:
            # Single solution, convert to list
            solution_list = [solutions]

        # Determine if the system is underdetermined
        num_params = len(actual_params)
        num_constraints = 1  # We have one equation

        if num_params > num_constraints:
            status = 'underdetermined'
        else:
            status = 'solved'

        return {
            'solutions': solution_list,
            'expression': scalar_product,
            'status': status,
            'parameters': actual_params
        }

    except Exception as e:
        # print(f"Error solving for cancellations: {e}") # Removed for cleaner output
        return {
            'solutions': [],
            'expression': scalar_product,
            'status': 'no_solution',
            'parameters': '--'
        }


def analyze_cancellations(X, A, Ba, U):
    """
    Complete analysis: compute m and check for cancellations.

    Returns:
        dict with both m computation and cancellation analysis
    """
    # First compute m
    m_value = compute_m_symbolic(X, A, Ba, U)

    if m_value is None:
        return {
            'm_value': None,
            'cancellation_analysis': None,
            'status': 'failed_m_computation'
        }

    # Then check for cancellations
    cancellation_result = check_cancellations(X, A, Ba, U, m_value)

    return {
        'm_value': m_value,
        'cancellation_analysis': cancellation_result,
        'status': 'completed'
    }


def functional_to_latex(A, Ba, U, root):
    """Convert the Lyapunov functional to LaTeX format with symbolic m computation."""
    from sympy import symbols

    dx = symbols('dx')

    base_term = r"\frac{1}{2}\|\mathbf{u}\|^2"
    terms_by_level = {0: [base_term]}
    cancellation_summaries = []  # List to store cancellation analysis results

    def collect_terms(node, level=1):

        if node.parent and node.parent.direction is not None:
            if level not in terms_by_level:
                terms_by_level[level] = []

            xi_exp = 2 * (1 + node.number)
            xi_factor = f"\\frac{{1}}{{\\xi^{{{xi_exp}}}}}"

            if node.parent.direction == 1:
                vec1 = node.parent.matrix * U
                vec2 = node.matrix * dx * U
                scalar_product = generate_l2_latex(vec1, vec2)
                if scalar_product != "0":
                    term = f"{xi_factor}\\left({scalar_product}\\right)"
                    terms_by_level[level].append(term)

            elif node.parent.direction == -1:
                vec1 = node.parent.matrix * U
                vec2 = node.matrix * U
                scalar_product = generate_l2_latex(vec1, vec2)
                if scalar_product != "0":
                    term = f"{xi_factor}\\left({scalar_product}\\right)"
                    terms_by_level[level].append(term)

            elif node.parent.direction == 0:
                # Compute m symbolically and check for cancellations
                X = node.parent.matrix
                analysis = analyze_cancellations(X, A, Ba, U)
                cancellation_summaries.append({
                    'node_info': f"Node at level {level}, name {node.name}",
                    'analysis_result': analysis
                })

                if analysis['status'] == 'completed':
                    m_symbolic = analysis['m_value']

                    # Term 1: main term
                    vec1_1 = node.parent.matrix * U
                    vec2_1 = node.matrix * dx * U
                    term1_latex = generate_l2_latex(vec1_1, vec2_1)

                    # Term 2: mixed term with symbolic m
                    vec1_2 = node.parent.matrix * A * U
                    vec2_2 = node.parent.matrix * Ba * A * U
                    term2_latex = generate_l2_latex(vec1_2, vec2_2)

                    if term1_latex != "0":
                        term1_line = f"{xi_factor}\\left({term1_latex}\\right)"
                        terms_by_level[level].append(term1_line)

                    if term2_latex != "0":
                        # Convert symbolic m to LaTeX
                        from sympy import latex
                        m_latex = latex(m_symbolic)
                        if m_latex == '1':
                            term2_line = f"{xi_factor}\\left({term2_latex}\\right)"
                        else:
                            term2_line = f"{xi_factor}{m_latex}\\left({term2_latex}\\right)"
                        terms_by_level[level].append(term2_line)

                else:
                    # Fall back to original behavior
                    vec1_1 = node.parent.matrix * U
                    vec2_1 = node.matrix * dx * U
                    term1_latex = generate_l2_latex(vec1_1, vec2_1)

                    vec1_2 = node.parent.matrix * A * U
                    vec2_2 = node.parent.matrix * Ba * A * U
                    term2_latex = generate_l2_latex(vec1_2, vec2_2)

                    if term1_latex != "0":
                        term1_line = f"{xi_factor}\\left({term1_latex}\\right)"
                        terms_by_level[level].append(term1_line)

                    if term2_latex != "0":
                        term2_line = f"{xi_factor}m\\left({term2_latex}\\right)"
                        terms_by_level[level].append(term2_line)

        for child in node.children:
            collect_terms(child, level + 1)

    # Collect all terms
    for child in root.children:
        collect_terms(child, level=1)


    # Build LaTeX output
    latex_output = r"\begin{align*}" + "\n"
    latex_output += r"\mathcal{L} &= " + terms_by_level[0][0]

    for level in sorted(terms_by_level.keys()):
        if level == 0:
            continue
        if terms_by_level[level]:
            filtered_terms = [t for t in terms_by_level[level]
                              if t != "0" and t.strip() != r"\left(0\right)"]
            if filtered_terms:
                latex_output += r" \\" + "\n" + r"&\quad+ " + (r" \\" + "\n" + r"&\quad+ ").join(filtered_terms)

    latex_output += r"\end{align*}"

    # --- Print cancellation summary outside the recursion ---
    print("\n" + "="*30)
    print("Cancellation Analysis Summary")
    print("="*30)

    if not cancellation_summaries:
        print("No mixed terms (direction 0) found for cancellation analysis.")
    else:
        for summary in cancellation_summaries:
            node_info = summary['node_info']
            analysis = summary['analysis_result']

            print(f"\n--- {node_info} ---")
            print(f"Computed m = {analysis['m_value']}")

            if analysis['status'] == 'completed':
                cancellation_info = analysis['cancellation_analysis']
                print(f"Cancellation status: {cancellation_info['status']}")
                print(f"Cancellation parameters: {cancellation_info['parameters']}; "
                      f"Expression: {cancellation_info['expression']}")

                if cancellation_info['status'] == 'solved':
                    print(f"Parameter constraints for cancellation: {cancellation_info['solutions']}")
                elif cancellation_info['status'] == 'underdetermined':
                    print(f"Underdetermined system. Possible constraints: {cancellation_info['solutions']}")
                elif cancellation_info['status'] == 'always_zero':
                    print("Cancellation expression is always zero - no constraints needed.")
                else:
                    print("No parameter constraints found for cancellation.")
            else:
                print(f"Analysis status: {analysis['status']}. Could not complete m computation or cancellation check.")

    print("\n" + "="*30 + "\n")
    # --- End of summary printing ---

    return latex_output