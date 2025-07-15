from sympy import *


class Operator:
    def __init__(self, A, Ba, Bs, X, u, d_x_u, direction, size):
        self.A = A
        self.Ba = Ba
        self.Bs = Bs
        self.X = X
        self.u = u
        self.d_x_u = d_x_u
        self.direction = direction
        self.size = size

    def compute_product_Ba(self):
        # Compute the products
        term1 = (self.X * self.Ba * self.u).dot(self.X * self.Ba * self.u)
        term2 = (self.X * self.Bs * self.u).dot(self.X * self.Ba * self.u)
        term3 = (self.X * self.A * self.d_x_u).dot(self.X * self.Ba * self.u)
        term4 = (self.X * self.d_x_u).dot(self.X * self.Ba * self.A * self.u)
        term5 = (self.X * self.u).dot(self.X * self.Ba * self.Bs * self.u)
        term6 = (self.X * self.u).dot(self.X * self.Ba * self.Ba * self.u)

        result = - term1 - term2 - term3 + term4 - term5 - term6
        result = expand(result)

        # Remove terms of the form d_x_u{i}*u{i}
        for i in range(1, self.size + 1):
            result = result.subs(self.d_x_u[i - 1] * self.u[i - 1], 0)

        # Put the derivative on the smallest index
        for i in range(1, self.size + 1):
            for j in range(i + 1, self.size + 1):
                result = result.subs(self.d_x_u[j - 1] * self.u[i - 1],
                                         self.d_x_u[i - 1] * self.u[j - 1])

        # Collect and simplify the expression
        result = collect(result, self.d_x_u)

        # Return the sum of terms
        return result

    def compute_product_A(self):
        # Compute the products
        term1 = (self.X * self.A * self.d_x_u).norm() ** 2
        term2 = simplify((self.X * self.Bs * self.d_x_u).dot(self.X * self.A * self.u))
        term3 = simplify((self.X * self.Ba * self.d_x_u).dot(self.X * self.A * self.u))
        term4 = simplify((self.X * self.d_x_u).dot(self.X * self.A ** 2 * self.d_x_u))
        term5 = simplify((self.X * self.d_x_u).dot(self.X * self.A * self.Bs * self.u))
        term6 = simplify((self.X * self.d_x_u).dot(self.X * self.A * self.Ba * self.u))

        result = - term1 + term2 + term3 + term4 + term5 + term6

        # Remove terms of the form d_x_u{i}*u{i}
        for i in range(1, self.size + 1):
            result = result.subs(self.d_x_u[i - 1] * self.u[i - 1], 0)

        # Put the derivative on the smallest index
        for i in range(1, self.size + 1):
            for j in range(i + 1, self.size + 1):
                result = result.subs(self.d_x_u[j - 1] * self.u[i - 1],
                                         self.d_x_u[i - 1] * self.u[j - 1])

        # Collect and simplify the expression
        result = collect(result, self.d_x_u)

        # Return the sum of terms
        return result

    def compute_product_mix(self):

        #Compute the products
        term1 = simplify((self.X * self.A * self.d_x_u).dot((self.X * self.Ba * self.A * self.u)))
        term2 = simplify((self.X * self.Ba * self.A * self.d_x_u).dot(self.X * self.A ** 2 * self.u))
        term3 = simplify((self.X * self.A * self.Ba * self.u).dot(self.X * self.Ba * self.A * self.u))
        term4 = simplify((self.X * self.A * self.Bs * self.u).dot(self.X * self.Ba * self.A * self.u))
        term5 = simplify((self.X * self.A * self.u).dot(self.X * self.Ba * self.A * self.Ba * self.u))
        term6 = simplify((self.X * self.A * self.u).dot(self.X * self.Ba * self.A * self.Bs * self.u))

        result = - term1 + term2 - term3 + term4 - term5 - term6

        # Remove terms of the form d_x_u{i}*u{i}
        for i in range(1, self.size + 1):
            result = result.subs(self.d_x_u[i - 1] * self.u[i - 1], 0)

        # Put the derivative on the smallest index
        for i in range(1, self.size + 1):
            for j in range(i + 1, self.size + 1):
                result = result.subs(self.d_x_u[j - 1] * self.u[i - 1],
                                 self.d_x_u[i - 1] * self.u[j - 1])

        # Collect and simplify the expression
        result = collect(result, self.d_x_u)
        result = result.simplify()

        # Return the sum of terms
        return result
    def compute_selected_product(self):
        if self.direction == 1:
            return self.compute_product_A()
        elif self.direction == -1:
            return self.compute_product_Ba()
        elif self.direction == 0:
            return self.compute_product_mix()
        else:
            raise ValueError("Invalid direction. Direction must be 1 or -1.")