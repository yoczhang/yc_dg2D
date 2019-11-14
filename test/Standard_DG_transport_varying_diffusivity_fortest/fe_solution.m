function r=fe_solution(x,y,solution,der_order_x,der_order_y,basis_type)


r=triangular_local_basis(x,y,basis_type,der_order_x,der_order_y)*solution;